# This is based on the `process_pdac_data` function (L141) in pdac_utils.py of Hwang et al.'s repo

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

# ===== IO =====
h5ad_in_dir = '/home/zhangh5/work/prePOLAR/hz-analysis/outputs/Step1.0-individual_named_h5s'
h5_out_dir = '/home/zhangh5/work/prePOLAR/hz-analysis/outputs/Step1.2-combined_H5'
output_prefix = 'post_cellbender_all_samples-COMBINED'

sc.settings.figdir = '/home/zhangh5/work/prePOLAR/hz-analysis/outputs/Step1.2-combined_QC_norm_basic_plots'

from Step1_qc_norm_cluster import karthik_combine_and_process_pdac_data, karthik_annotate_data, karthik_plot_QC, karthik_analyze_data, karthik_more_plots, calculate_qc_metrics

# ----- 0. Preprocess all samples' data individually and combine -----
combined_adata, combined_adata_raw = karthik_combine_and_process_pdac_data(
    h5ad_in_dir = h5ad_in_dir, 
    output_prefix = output_prefix, 
    h5_out_dir = h5_out_dir,
)

# # ----- 1. Load the combined data and perform joint analysis -----
# # combined_adata = sc.read_h5ad(f'/home/zhangh5/work/prePOLAR/analysis/outputs/Step1-combined_H5/{data_type}_adata.h5ad')

# combined_adata = calculate_qc_metrics(combined_adata)
# combined_adata = karthik_annotate_data(combined_adata)

# #sc.pp.highly_variable_genes(combined_adata, batch_key='batch')

# highly_variable = []

# for sample_name in combined_adata.obs['batch'].unique():

#     result = sc.pp.highly_variable_genes(combined_adata[combined_adata.obs["batch"]==sample_name, :], min_mean=0.0125, max_mean=3, min_disp=0.5, inplace=False)
#     highly_variable.append(result['highly_variable'])

# highly_variable = np.array(highly_variable)
# combined_adata.var['highly_variable'] = np.sum(highly_variable, axis=0) > 6

# topgenes = combined_adata.var_names[np.argsort(np.sum(combined_adata.X, axis=0)/np.sum(combined_adata.X))].flatten()[-50:] # @HZ: fixed bug...

# sc.tl.rank_genes_groups(combined_adata, "batch", method='t-test', n_genes=500, use_raw=True)

# patientspecificgenes = set()
# for genes in combined_adata.uns['rank_genes_groups']['names']:
#     for gene in genes:
#         patientspecificgenes.add(gene)

# #embed()
# highest_expressed = []
# for gene in combined_adata.var_names:
#     if gene in topgenes: # or gene in patientspecificgenes:
#         highest_expressed.append(True)
#     else:
#         highest_expressed.append(False)
# highest_expressed = np.array(highest_expressed)

# combined_adata.var['highly_variable'] = ~(combined_adata.var_names.str.startswith("MT-")) & ~(combined_adata.var_names.str.startswith("RP")) & (combined_adata.var['highly_variable']) & ~(highest_expressed)

# print(len(patientspecificgenes), sum(combined_adata.var['highly_variable']))

# # --- save the raw count
# rawcount_adata = AnnData(X=combined_adata.layers['counts'])
# rawcount_adata.obs_names = combined_adata.obs_names
# rawcount_adata.var_names = combined_adata.var_names
# #rawcount_adata = rawcount_adata[combined_adata.obs['cells_tokeep']==True, :]
# rawcount_adata.write(str(Path(h5_out_dir) / f'{data_type}_raw_adata.h5ad'))

# #combined_adata = combined_adata[combined_adata.obs['cells_tokeep']==True, :]
# # combined_adata = karthik_plot_QC(combined_adata, sample_name = 'all_PP_combined')
# combined_adata = karthik_analyze_data(combined_adata)

# sc.settings.figdir = sc.settings.figdir / 'UMAPs'
# Path(sc.settings.figdir).mkdir(parents=True, exist_ok=True)
# combined_adata = karthik_more_plots(combined_adata, sample_name = 'all_PP_combined')
    
# del combined_adata.uns['rawcount']

# combined_adata.write(str(Path(h5_out_dir) / f'{data_type}_adata.h5ad'))