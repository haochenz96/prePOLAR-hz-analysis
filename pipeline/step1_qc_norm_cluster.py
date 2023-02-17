# data
import scanpy as sc
import numpy as np
import scrublet as scr
from anndata import AnnData
from utils import load_h5_data
from glob import glob
from pathlib import Path
import logging
logger = logging.getLogger(__name__)


# plotting
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import seaborn.objects as so
from scanpy.plotting.palettes import *
my_palette = godsnot_102

# reference
from fetch_references import load_karthik_genemarkers, get_broad_celltypes

from IPython import embed

# # R
# # Might not be necessary but about setting R_HOME in conda environment: https://stackoverflow.com/questions/38194876/how-to-set-the-r-home-environment-variable-to-the-r-home-directory
# import rpy2
# import anndata2ri
# import rpy2.robjects as ro
# import rpy2.rinterface_lib.callbacks as rcb
# from rpy2.robjects.packages import importr
# base = importr('base')
# rcb.logger.setLevel(logging.ERROR)
# ro.pandas2ri.activate()
# anndata2ri.activate()

# ----- calculate QC metrics -----

def calculate_qc_metrics(adata):
    # 1. label MT, RB, HB genes
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

    # 2. compute the doublet scores
    scrub = scr.Scrublet(adata.X)
    doublet_scores, _ = scrub.scrub_doublets()
    adata.obs['scrublet_scores'] = doublet_scores
    
    #adata = adata[adata.obs['scrublet_scores'] <= 0.2]
    
    # 3. calculate the saturation counts
    sc.pp.calculate_qc_metrics(adata, percent_top=(5, 10, 20, 50, 100, 200, 500), inplace=True)
    adata.obs['n_counts_sat'] = [1000]*adata.shape[0]
    adata.obs['n_counts_sat'] = np.min(adata.obs[['total_counts', 'n_counts_sat']], axis=1)
    adata.obs['n_genes_sat'] = [1000]*adata.shape[0]
    adata.obs['n_genes_sat'] = np.min(adata.obs[['n_genes_by_counts', 'n_genes_sat']], axis=1)
    mito_genes = adata.var_names.str.startswith('MT-')
    adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    
    return adata

# ----- basic QC plots -----
def basic_plots(adata, out_pdf, show_outlier = None):
    if not all([ x in adata.obs.columns for x in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'] ]):
        raise ValueError("adata.obs must contain columns 'n_genes_by_counts', 'total_counts', 'pct_counts_mt'")
    
    # 1. Top 20 highest expressed genes
    fig1, ax = plt.subplots(figsize=(8,12))
    fig1.tight_layout(pad=5)
    p_box = sc.pl.highest_expr_genes(adata, n_top=20, show=False, ax=ax)
    
    # 2. scatter QC metrics

    # https://github.com/scverse/scanpy/issues/2135
    p_scatter = sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color="pct_counts_mt",ax=None, show=False)
    fig2 = plt.gcf()
    fig2.set_size_inches(8, 8)

    # 3. violin QC metrics
    fig3, axes = plt.subplots(1,3, figsize = (12, 8))
    fig3.tight_layout(pad=3)
    p_violin = sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4, show=False, ax=axes[0])
    
    p_violin = sc.pl.violin(adata, ['total_counts'], jitter=0.4, show=False, ax=axes[1])
    p_violin = sc.pl.violin(adata, ['pct_counts_mt'], jitter=0.4, show=False, ax=axes[2])

    if show_outlier is not None:
        if type(show_outlier) is not int:
            raise ValueError("show_outlier must be an integer")
        # show outliers on Violin plots
        for i, prop in zip(range(3), ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']):
            outlier_threshs = get_outlier_thresholds(adata.obs[prop], show_outlier)

            for y in outlier_threshs:
                if y is not None:
                    axes[i].axhline(y=y, xmin = 0.1, xmax = 0.9, color='r', linestyle='--', linewidth = 1, label = f"outlier ({show_outlier} MADs) threshold")

                axes[i].legend(bbox_to_anchor = (1.0, 1), loc = 'upper right', prop={'size':5})

    # save to pdf
    with PdfPages('foo.pdf') as out_pdf:
        out_pdf.savefig(fig1)
        out_pdf.savefig(fig2)
        out_pdf.savefig(fig3)

# ----- filters -----

# ====== from Hwang et al., Nature Genetics (2021) =====

def karthik_basic_filter(adata):
    #sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=0)
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    #pdac_thresh_dict = load_umithreshold()
    #min_counts = pdac_thresh_dict.get(pid, 200)
    #adata = adata[adata.obs['n_counts'] >= min_counts]
    #adata = adata[adata.obs['n_counts'] <= 2500]
    return adata

# ====== from sc-best-practices.org =====

# automatic thresholding via MAD (median absolute deviations)
def get_outlier_thresholds(M: np.array, nmads: int):
    '''
    M: numpy.array
        1D array of values
    nmads: int
        number of median absolute deviations (MADs) to use as threshold

    returns: Tuple
        (lower threshold, upper threshold)
    '''

    threshs = (np.median(M) - nmads * M.mad(), np.median(M) + nmads * M.mad())
    if threshs[0] < M.min():
        threshs[0] = None
    if threshs[1] > M.max():
        threshs[1] = None
    return threshs

# ===== from Yubin =====
def load_yubin_filters():
    yubin_filters = {
        'min_genes_per_cell': 127, 
        'min_counts_per_cell': 255, 
        'min_cells_per_gene': 10, 
        'max_counts_per_gene': None, 
        'remove_cells_with_high_mt': 5, 
        'remove_mt_genes': True,
        'remove_rb_genes': True,
    }

    return yubin_filters

def hard_filter_sample(adata, filters: dict):
    '''

    '''
    orig_data_shape = adata.shape

    if (min_counts_per_cell := filters.get('min_counts_per_cell')):
        sc.pp.filter_cells(adata, min_counts=min_counts_per_cell)
        # logger.info(f"Filtered {np.sum(adata.obs['total_counts'] > min_counts_per_cell)} / {adata.obs.shape[0]} cells with less than {min_genes_per_cell} genes expressed")

    if (min_genes_per_cell := filters.get('min_genes_per_cell')):
        sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)
        # logger.info(f'Filtered {} / {} cells with less than {min_genes_per_cell} genes expressed')

    if (max_counts_per_gene := filters.get('max_counts_per_gene')):
        sc.pp.filter_genes(adata, max_counts=max_counts_per_gene)

    if (min_cells_per_gene := filters.get('min_cells_per_gene')):
        sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    
    if (remove_cells_with_high_mt := filters.get('remove_cells_with_high_mt')):
        if not 'pct_counts_mt' in adata.obs.columns:
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
        adata = adata[adata.obs['pct_counts_mt'] < remove_cells_with_high_mt, :]
    if (remove_mt_genes := filters.get('remove_mt_genes')):
        adata = adata[:, ~adata.var_names.str.startswith('mt-')]
    if (remove_rb_genes := filters.get('remove_rb_genes')):
        adata = adata[:, ~adata.var_names.str.startswith('RPS') | ~adata.var_names.str.startswith('RPL')]
    if (remove_hb_genes := filters.get('remove_hb_genes')):
        adata = adata[:, ~adata.var_names.str.contains(("^HB[^(P)]"))]

    filtered_data_shape = adata.shape
    logging.info(f"Filtered {orig_data_shape[0] - filtered_data_shape[0]} cells and {orig_data_shape[1] - filtered_data_shape[1]} genes")

    return adata


# ----- normalization -----

# (1)
# Hwang et al log1p; scanpy function `normalize_per_cell` is deprecated though
# Karthik method
def karthik_normalize_data(adata):
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    return adata

def log1p_normalize_data(adata):
    # proportional fitting to mean of cell depth
    proportional_fitting = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    # log1p transform
    adata.layers["log1pPF_normalization"] = sc.pp.log1p(proportional_fitting["X"])
    # proportional fitting
    adata.layers["PFlog1pPF_normalization"] = sc.pp.normalize_total(
        adata, target_sum=None, layer="log1pPF_normalization", inplace=False
    )["X"]

# # (2) scran
# def scran_normalize_data(adata):
#     import scipy
#     from scipy.sparse import csr_matrix, issparse

#     scran = importr("scran")
#     BiocParallel = importr("BiocParallel")

#     # Preliminary clustering for differentiated normalisation
#     adata_pp = adata.copy()
#     sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
#     sc.pp.log1p(adata_pp)
#     sc.pp.pca(adata_pp, n_comps=15)
#     sc.pp.neighbors(adata_pp)
#     sc.tl.leiden(adata_pp, key_added="groups")

#     data_mat = adata_pp.X.T
#     # convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
#     if scipy.sparse.issparse(data_mat):
#         if data_mat.nnz > 2**31 - 1:
#             data_mat = data_mat.tocoo()
#         else:
#             data_mat = data_mat.tocsc()
#     ro.globalenv["data_mat"] = data_mat
#     ro.globalenv["input_groups"] = adata_pp.obs["groups"]

#     size_factors = scran.sizeFactors(
#         scran.computeSumFactors(
#             scran.SingleCellExperiment(
#                 list(counts=data_mat)), 
#                 clusters = input_groups,
#                 min.mean = 0.1,
#                 BPPARAM = BiocParallel.MulticoreParam()
#         )
#     )
#     del adata_pp

    



# ----- plot highly variable genes data -----
def plot_hvg(adata, sample_name = '', write_figure=True):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=1000)
    if write_figure:
        out_f = f'{sample_name}-HVG_plot.png'
        sc.pl.highly_variable_genes(adata, show=False, save=out_f)

# ----- cluster data (pca, neighbors, louvain, umap) -----
def cluster_data(adata, sample_name = '', output_dir = None):
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.leiden(adata)
    sc.tl.umap(adata)

    if output_dir is not None:
        out_f_prefix = str(Path(output_dir) / f'{sample_name}-')
    else:
        out_f_prefix = f'{sample_name}-'

    sc.pl.pca(adata, color='CST3', show=False, save=f'{out_f_prefix}pca.png')
    sc.pl.pca_variance_ratio(adata, log=True, show=False, save=f'{out_f_prefix}pcavariance.png')


# ----- 

# ===== Karthik wrappers =====
def karthik_annotate_data(adata):
    print(adata)
    mito_genes = adata.var_names.str.startswith('MT-')
    adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    
    # standard normalization
    adata = adata[adata.obs['n_genes'] < 2500, :]
    
    # rawcountdata = AnnData(X=adata.X)
    # rawcountdata.obs_names = adata.obs_names
    # rawcountdata.var_names = adata.var_names
    # adata.uns['rawcount'] = rawcountdata # <---- @HZ: problematic
    
    #adata = adata[adata.obs['percent_mito'] < 0.05, :]
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    return adata

def karthik_plot_QC(adata, sample_name):

    sc.pl.highest_expr_genes(adata, n_top=20, showfliers=False, save=f'{sample_name}-QC-HEG.pdf')
    sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save=f'{sample_name}-QC-violins.pdf')
    sc.pl.scatter(adata, x='n_counts', y='percent_mito', save=f'{sample_name}-QC-scatter_mito.pdf')
    sc.pl.scatter(adata, x='n_counts', y='n_genes', save=f'{sample_name}-QC-scatter_ngenes.pdf')
    #sc.pl.highly_variable_genes(adata)
    return adata

def karthik_analyze_data(adata, resolution=None):

    genemarkers = load_karthik_genemarkers('/home/zhangh5/work/prePOLAR/analysis/reference/karthik-gene_markers.use.txt')
    broad_celltypes = get_broad_celltypes()
    
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    if resolution:
        sc.tl.leiden(adata, resolution=resolution)
    else:
        sc.tl.leiden(adata, resolution)
        sc.tl.umap(adata)
        for broad_celltype, specific_celltypes in broad_celltypes.items():
            allmarkers = set()
            for celltype in specific_celltypes:
                markers = genemarkers[celltype]
                allmarkers = allmarkers.union(set(markers))
                sc.tl.score_genes(adata, markers, score_name=celltype, use_raw=True)
            sc.tl.score_genes(adata, allmarkers, score_name=broad_celltype, use_raw=True)
        sc.tl.rank_genes_groups(adata, 'leiden', method='t-test', use_raw=True)
    return adata

def karthik_more_plots(adata, sample_name):
    COLORS = my_palette


    broad_celltypes = get_broad_celltypes()
    sc.pl.pca(adata, color=['percent_mito', 'n_genes', 'n_counts'], cmap="jet", save=f'{sample_name}-PCA-QC.pdf')
    sc.pl.pca_variance_ratio(adata, log=True, save=f'{sample_name}-PCA-var_ratio.pdf')
    sc.pl.umap(adata, color=['n_genes', 'scrublet_scores', 'batch', 'leiden'], palette=COLORS, save=f'{sample_name}-UMAP-QC.pdf')
    sc.pl.umap(adata, color=['CFTR', 'ANXA2', 'IL7R', 'KRT19', 'PTPRC', 'PECAM1', 'ACTA2', 'CD96', 'MRC1', 'CD163', 'KRAS'],  cmap="jet", save=f'{sample_name}-UMAP-marker_genes.pdf')

    print ("Broad Cell types")
    sc.pl.umap(adata, color=broad_celltypes.keys(), cmap="jet")
    
    for broad_celltype, specific_celltypes in broad_celltypes.items():
        if len(specific_celltypes) > 1:
            print(broad_celltype)
            sc.pl.umap(adata, color=specific_celltypes, cmap="jet", save=f'{sample_name}-UMAP-{broad_celltype}.pdf')

    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    return adata

def karthik_process_pdac_data(h5_in_dir, data_type, h5_out_dir):
    #filenames = glob("/ahg/regevdata/projects/Pancreas/alignreads/*_10x/*-NSTnPo*/raw_*.h5")
    h5_fs = glob(f"{h5_in_dir}/*.h5")
    combined_data = []
    for filename in h5_fs:
        naivedata = load_h5_data(filename)
        combined_data.append(naivedata)

    # direct concatentation of all data
    naivedata = combined_data[0].concatenate(combined_data[1:])

    naivedata = karthik_basic_filter(naivedata)
    # naivedata.write(str(Path(h5_out_dir) / f'{data_type}_adata.h5ad'))
    naivedata = karthik_annotate_data(naivedata)

    rawcount_adata = adata.layers['counts']
    rawcount_adata.write(str(Path(h5_out_dir) / f'{data_type}_raw_adata.h5ad'))

    highly_variable = []
    for sample_name in naivedata.obs['batch'].unique():

        result = sc.pp.highly_variable_genes(naivedata[naivedata.obs["batch"]==sample_name, :], min_mean=0.0125, max_mean=3, min_disp=0.5, inplace=False)
        highly_variable.append(result['highly_variable'])

    highly_variable = np.array(highly_variable)
    naivedata.var['highly_variable'] = np.sum(highly_variable, axis=0) > 6
    
    topgenes = list(naivedata.var_names[np.argsort(np.sum(naivedata.X, axis=0)/np.sum(naivedata.X))])[-50:]
    
    sc.tl.rank_genes_groups(naivedata, "batch", method='t-test', n_genes=500, use_raw=True)
    
    patientspecificgenes = set()
    for genes in naivedata.uns['rank_genes_groups']['names']:
        for gene in genes:
            patientspecificgenes.add(gene)

    highest_expressed = []
    embed()
    for gene in naivedata.var_names:
        if gene in topgenes: # or gene in patientspecificgenes:
            highest_expressed.append(True)
        else:
            highest_expressed.append(False)
    highest_expressed = np.array(highest_expressed)
    
    naivedata.var['highly_variable'] = ~(naivedata.var_names.str.startswith("MT-")) & ~(naivedata.var_names.str.startswith("RP")) & (naivedata.var['highly_variable']) & ~(highest_expressed)
    
    print(len(patientspecificgenes), sum(naivedata.var['highly_variable']))

    naivedata = karthik_plot_QC(naivedata)
    naivedata = karthik_analyze_data(naivedata)

    sc.settings.figdir = sc.settings.figdir + '/UMAPs'
    Path(sc.settings.figdir).mkdir(parents=True, exist_ok=True)
    naivedata = karthik_more_plots(naivedata)
        
    del naivedata.uns['rawcount']
    
    naivedata.write(str(Path(h5_out_dir) / f'{data_type}_adata.h5ad'))
    return naivedata, rawcount_adata