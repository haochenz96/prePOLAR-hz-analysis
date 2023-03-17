# This is based on the `process_individual` function (L477) in pdac_utils.py of Hwang et al.'s repo

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
from Step1_qc_norm_cluster import karthik_basic_filter, calculate_qc_metrics, load_yubin_filters, hard_filter_sample
# --> @normalize
from Step1_qc_norm_cluster import karthik_normalize_data, plot_hvg
# --> @cluster
from Step1_qc_norm_cluster import cluster_data
# --> @celltype
from fetch_references import load_karthik_genemarkers, get_broad_celltypes
karthik_genemarkers = load_karthik_genemarkers('/home/zhangh5/work/prePOLAR/hz-analysis/reference/karthik-gene_markers.use.txt')
broad_celltypes = get_broad_celltypes()


# ===== IO =====
h5_in_dir = '/home/zhangh5/work/prePOLAR/hz-analysis/outputs/Step0-CellBender-DONE_AT_LILAC/processed_h5s'
h5_fs = glob(f"{h5_in_dir}/*filtered_filtered.h5")

# # PP01 for testing
# h5_fs = [
#     '/home/zhangh5/work/prePOLAR/hz-analysis/outputs/Step0-CellBender-DONE_AT_LILAC/processed_h5s/PP01.cellbender_filtered_filtered.h5'
# ]

individual_h5_out_dir = Path('/home/zhangh5/work/prePOLAR/hz-analysis/outputs/Step1.0-individual_named_h5s')
individual_h5_out_dir.mkdir(parents = True, exist_ok = True)
fig_out_dir = Path('/home/zhangh5/work/prePOLAR/hz-analysis/outputs/Step1.1-individual_QC_normal_basic_plots')
fig_out_dir.mkdir(parents = True, exist_ok = True)
# ==============

for h5_i in h5_fs:

    sample_name_i = Path(h5_i).name.split('.')[0]
    if (individual_h5_out_dir / f'{sample_name_i}.h5ad').exists():
        print(f"[INFO] {sample_name_i} has been processed, skipping.")
        continue

    print(f"processing {sample_name_i} ...")

    adata = load_h5_data(h5_i)
    print(f'Before filtering, adata is of shape: {adata.shape}')

    # ----- save post-cellbender h5 immediately after sample_name_i is added -----
    adata.obs['sample_name'] = sample_name_i
    adata.write(f'{individual_h5_out_dir}/{sample_name_i}.h5ad')

    # ----- set figdir -----
    sc.settings.figdir = fig_out_dir / sample_name_i
    sc.settings.figdir.mkdir(parents = True, exist_ok = True)

    # ----- filter, calculate QC metrics -----
    adata_filtered = karthik_basic_filter(adata) # this adds the n_genes, n_counts layers
    adata_filtered = calculate_qc_metrics(adata) # this adds more metadata layers

    # # ----- Yubin filters ----- @HZ 03/07/2023: skip for now
    # yubin_filters = load_yubin_filters()
    # adata_filtered = hard_filter_sample(adata, yubin_filters)
    # print(f'After filtering, adata is of shape: {adata_filtered.shape}')
    # if adata_filtered.shape[0] < 200:
    #     print(f'[WARNING] {sample_name_i} has less than 200 cells after filtering, skip it.')
    #     continue

    # ----- log1p normalize -----
    adata_filtered = karthik_normalize_data(adata_filtered)

    # ----- plot highly variable genes -----
    plot_hvg(adata_filtered, sample_name = sample_name_i,)

    # ----- cluster and make the PCA plots -----
    cluster_data(adata_filtered, sample_name = sample_name_i,) # this creates the PCA and PCA variant plots
    # embed()
    # ----- score cells based on gene signatures -----
    for broad_celltype, specific_celltypes in broad_celltypes.items():
        allmarkers = set()
        for celltype in specific_celltypes:
            markers = karthik_genemarkers[celltype]
            allmarkers = allmarkers.union(set(markers))
            # print(f"[INFO] ----> markers: {markers}")
            sc.tl.score_genes(adata_filtered, markers, score_name=celltype, use_raw=True)
        sc.tl.score_genes(adata_filtered, allmarkers, score_name=broad_celltype, use_raw=True)

    # ----- plot UMAPs -----
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
            print(f"[INFO] making UMAP for cell type: {broad_celltype}")
            sc.pl.umap(adata_filtered, color=specific_celltypes, palette=my_palette, show=False, save=f"{sample_name_i}-UMAP-broad_specific-{broad_celltype}.pdf")
