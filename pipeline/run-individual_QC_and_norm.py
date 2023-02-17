import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import scrublet as scr
import seaborn as sns
import seaborn.objects as so
from utils import load_h5_data
from pathlib import Path
from glob import glob

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) # suppress FutureWarning

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    # color_map="YlGnBu",
    frameon=False,
)
from scanpy.plotting.palettes import *
my_palette = godsnot_102
from matplotlib.pyplot import rc_context
from matplotlib.backends.backend_pdf import PdfPages

from IPython import embed

# ===== functions =====
# --> @QC
from QC_and_norm import calculate_qc_metrics, load_yubin_filters, hard_filter_sample
# --> @normalize
# Karthik method
def karthik_normalize_data(adata):
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    return adata

def plot_hvg(adata, sample_name = '', write_figure=True):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=1000)
    if write_figure:
        out_f = f'{sample_name}-HVG_plot.pdf'
        sc.pl.highly_variable_genes(adata, show=False, save=out_f)
# --> @cluster
def cluster_data(adata, sample_name = '', write_figure=True):
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.leiden(adata)
    sc.tl.umap(adata)
    
    if write_figure:
        out_f_prefix = f'{sample_name}-'
        sc.pl.pca(adata, color='CST3', show=False, save=f'{out_f_prefix}pca.pdf') # gene marker is hard coded to `CST3` which is a marker of Schwann cells and Langherans cells
        sc.pl.pca_variance_ratio(adata, log=True, show=False, save=f'{out_f_prefix}pca_variance.pdf')
# --> @celltype
from fetch_references import load_karthik_genemarkers, get_broad_celltypes
karthik_genemarkers = load_karthik_genemarkers('/home/zhangh5/work/prePOLAR/analysis/reference/karthik-gene_markers.use.txt')
broad_celltypes = get_broad_celltypes()


# ===== IO =====
h5_in_dir = '/home/zhangh5/work/prePOLAR/analysis/outputs/CellBender/processed_h5s'
h5_fs = glob(f"{h5_in_dir}/*.h5")

out_dir = Path('/home/zhangh5/work/prePOLAR/analysis/outputs/STEP1-individual_QC_normal_basic_plots')
# ==============

for h5_i in h5_fs:
    sample_name_i = Path(h5_i).name.split('.')[0]
    print(f"processing {sample_name_i} ...")

    adata = load_h5_data(h5_i)
    print(f'Before filtering, adata is of shape: {adata.shape}')

    
    sc.settings.figdir = out_dir / sample_name_i
    sc.settings.figdir.mkdir(parents = True, exist_ok = True)

    # ----- calculate QC metrics -----
    calculate_qc_metrics(adata)
    yubin_filters = load_yubin_filters()

    adata_filtered = hard_filter_sample(adata, yubin_filters)
    print(f'After filtering, adata is of shape: {adata_filtered.shape}')
    if adata_filtered.shape[0] < 200:
        print(f'[WARNING] {sample_name_i} has less than 200 cells after filtering, skip it.')
        continue

    # ----- log1p normalize -----
    karthik_normalize_data(adata_filtered)
    plot_hvg(adata_filtered, sample_name = sample_name_i,)

    # ----- cluster -----
    cluster_data(adata_filtered, sample_name = sample_name_i,)
    # embed()
    # score cells based on gene signatures
    for broad_celltype, specific_celltypes in broad_celltypes.items():
        allmarkers = set()
        for celltype in specific_celltypes:
            markers = karthik_genemarkers[celltype]
            allmarkers = allmarkers.union(set(markers))
            print(f"[INFO] ----> markers: {markers}")
            sc.tl.score_genes(adata_filtered, markers, score_name=celltype, use_raw=True)
        sc.tl.score_genes(adata_filtered, allmarkers, score_name=broad_celltype, use_raw=True)

    # ----- UMAPs -----
    sc.pl.umap(
        adata_filtered, 
        color=['leiden', 'n_genes', 'n_counts', 'scrublet_scores', 'n_counts_sat', 'percent_mito', 'pct_counts_in_top_5_genes', 'pct_counts_in_top_10_genes', 'pct_counts_in_top_20_genes', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes'], 
        palette=my_palette, 
        color_map='Reds',
        save = f"{sample_name_i}-UMAP-QC.pdf"
        )

    sc.pl.umap(adata_filtered, color=list(broad_celltypes.keys())+['scrublet_scores', 'KRT19', 'PTPRC', 'CFTR', 'COL1A1', 'PECAM1', 'leiden'], color_map='Reds', palette=my_palette, save = f"{sample_name_i}-UMAP-broad_celltype.pdf")

    for broad_celltype, specific_celltypes in broad_celltypes.items():
        if len(specific_celltypes) > 1:
            print (broad_celltype)
            sc.pl.umap(adata_filtered, color=specific_celltypes, palette=my_palette, show=False, save=f"{sample_name_i}-UMAP-broad_specific-{specific_celltypes}.pdf")
