import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from anndata import AnnData

from IPython import embed

# set Scanpy plotting parameters
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    # color_map="YlGnBu",
    frameon=False,
)
sc.settings.figdir = '/home/zhangh5/work/PREPOLAR/analysis/outputs/Step1-combined_QC_norm_basic_plots'
h5_out_dir = '/home/zhangh5/work/prePOLAR/analysis/outputs/Step1-combined_H5'
data_type = 'post_cellbender_all_samples-COMBINED'
from step1_qc_norm_cluster import karthik_process_pdac_data, karthik_annotate_data, karthik_plot_QC, karthik_analyze_data, karthik_more_plots, calculate_qc_metrics


# ----- 0. Preprocess all samples' data individually and combine -----
karthik_process_pdac_data(
    h5_in_dir = '/home/zhangh5/work/prePOLAR/analysis/outputs/CellBender/processed_h5s', 
    data_type = data_type, 
    h5_out_dir = '/home/zhangh5/work/prePOLAR/analysis/outputs/Step1-combined_H5',
)

# ----- 1. Load the combined data and perform joint analysis -----
naive_adata = sc.read_h5ad(f'/home/zhangh5/work/prePOLAR/analysis/outputs/Step1-combined_H5/{data_type}_adata.h5ad')

naive_adata = calculate_qc_metrics(naive_adata)
naive_adata = karthik_annotate_data(naive_adata)

#sc.pp.highly_variable_genes(naive_adata, batch_key='batch')

highly_variable = []

for sample_name in naive_adata.obs['batch'].unique():

    result = sc.pp.highly_variable_genes(naive_adata[naive_adata.obs["batch"]==sample_name, :], min_mean=0.0125, max_mean=3, min_disp=0.5, inplace=False)
    highly_variable.append(result['highly_variable'])

highly_variable = np.array(highly_variable)
naive_adata.var['highly_variable'] = np.sum(highly_variable, axis=0) > 6

topgenes = naive_adata.var_names[np.argsort(np.sum(naive_adata.X, axis=0)/np.sum(naive_adata.X))].flatten()[-50:] # @HZ: fixed bug...

sc.tl.rank_genes_groups(naive_adata, "batch", method='t-test', n_genes=500, use_raw=True)

patientspecificgenes = set()
for genes in naive_adata.uns['rank_genes_groups']['names']:
    for gene in genes:
        patientspecificgenes.add(gene)

#embed()
highest_expressed = []
for gene in naive_adata.var_names:
    if gene in topgenes: # or gene in patientspecificgenes:
        highest_expressed.append(True)
    else:
        highest_expressed.append(False)
highest_expressed = np.array(highest_expressed)

naive_adata.var['highly_variable'] = ~(naive_adata.var_names.str.startswith("MT-")) & ~(naive_adata.var_names.str.startswith("RP")) & (naive_adata.var['highly_variable']) & ~(highest_expressed)

print(len(patientspecificgenes), sum(naive_adata.var['highly_variable']))

# --- save the raw count
rawcount_adata = AnnData(X=naive_adata.layers['counts'])
rawcount_adata.obs_names = naive_adata.obs_names
rawcount_adata.var_names = naive_adata.var_names
#rawcount_adata = rawcount_adata[naive_adata.obs['cells_tokeep']==True, :]
rawcount_adata.write(str(Path(h5_out_dir) / f'{data_type}_raw_adata.h5ad'))

#naive_adata = naive_adata[naive_adata.obs['cells_tokeep']==True, :]
# naive_adata = karthik_plot_QC(naive_adata, sample_name = 'all_PP_combined')
naive_adata = karthik_analyze_data(naive_adata)

sc.settings.figdir = sc.settings.figdir / 'UMAPs'
Path(sc.settings.figdir).mkdir(parents=True, exist_ok=True)
naive_adata = karthik_more_plots(naive_adata, sample_name = 'all_PP_combined')
    
del naive_adata.uns['rawcount']

naive_adata.write(str(Path(h5_out_dir) / f'{data_type}_adata.h5ad'))