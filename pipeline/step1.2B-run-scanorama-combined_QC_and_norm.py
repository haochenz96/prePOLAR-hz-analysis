import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scanorama
from step2_plot_annotated_UMAPs import plot_annotated_umaps

from pathlib import Path
import pickle as pkl

# ===== IO =====
combined_h5ad = '/home/zhangh5/work/prePOLAR/hz-analysis/outputs/Step1.2-combined_H5/post_cellbender_all_samples-COMBINED_adata.h5ad'

joint_ann_f = '/home/zhangh5/work/prePOLAR/hz-analysis/HZ-prePOLAR-joint_sample_sheet.csv'

output_integrated_h5ad = '/home/zhangh5/work/prePOLAR/hz-analysis/outputs/Step1.2-combined_H5/COMBINED_adata.scanorama_corrected.h5ad'
output_dir = Path(output_integrated_h5ad).parent
output_filename = Path(output_integrated_h5ad).stem
pkl_file = output_dir / f'{output_filename}.pkl'

figure_output_dir = '/home/zhangh5/work/prePOLAR/hz-analysis/outputs/Step1.2B-combined_scanorama_corrected_H5_annotated_UMAPs'

# ===== run params =====
route = 'correct'

if Path(output_integrated_h5ad).exists():
    print(f'[INFO] {output_integrated_h5ad} already exists! Skipping integration and loading the existing file.')

    adata_sc = sc.read_h5ad(output_integrated_h5ad)
    
    ## tsne and umap
    print(f'[INFO] {output_integrated_h5ad} loaded successfully.')
    if not 'neighbors' in adata_sc.uns.keys():
        print(f'[INFO] Computing neighbors')
        sc.pp.neighbors(adata_sc, n_pcs =40, use_rep = "X_scanorama")
        sc.tl.umap(adata_sc)
        sc.tl.tsne(adata_sc, n_pcs = 40, use_rep = "X_scanorama")
        print(f'[INFO] UMAP and tSNE performed successfully. Saving the file to {output_integrated_h5ad}')
        adata_sc.write(
            output_integrated_h5ad
        )
    else:
        print(f'[INFO] Neighbors already computed. Skipping.')

    # load clinical datasheet
    joint_ann_df = pd.read_csv(joint_ann_f, index_col=0)

    # ===== UMAPs =====
    # --- palette
    from scanpy.plotting.palettes import *
    my_palette = godsnot_102

    # (1) sample_name and leiden
    fig, axes = plt.subplots(1,2) # this creates the canvas for us to put the plots on
    sc.pl.umap(adata_sc, color="sample_name", show=False, palette = my_palette, ax=axes[0])
    sc.pl.umap(adata_sc, color="leiden", show=False, ax=axes[1])
    fig.set_size_inches(14, 6)
    fig.tight_layout(pad=1) # this adjusts the spacing between the subplots
    fig.savefig(f'{figure_output_dir}/UMAP_sample_name_and_leiden.png')

    # (2) general cell types

    features_of_interest = ['FIBROBLASTS', 'IMMUNE', 'PANCREATIC_SCHWANN_CELLS', 'ENDOTHELIAL', 'ENDOCRINE', 'CFTR', 'DUCTAL', 'scrublet_scores']
    fig = plot_annotated_umaps(adata_sc, features_of_interest, num_cols=4)
    fig.savefig(f'{figure_output_dir}/UMAP_general_cell_types.png')

    # (3) clinical features
    features_of_interest = ['Bx.Timing', 'metastatic_first_line_regimen', 'OS', 'LTS_V_STS', 'BRCA']
    for feature_i in features_of_interest:
        sample_to_feature_map = joint_ann_df[feature_i].to_dict() # this creates a dictionary mapping sample names to the feature of interest. i.e. {sample_name: feature_value}

        # now we add this dictionary to the adata object as an annotation
        adata_sc.obs[feature_i] = adata_sc.obs['sample_name'].map(sample_to_feature_map)
    fig = plot_annotated_umaps(adata_sc, features_of_interest, num_cols=4)
    fig.savefig(f'{figure_output_dir}/UMAP_clinical_features.png')


elif pkl_file.exists():
    print(f'[INFO] Pickled file {pkl_file} already exists, loading it and performing concatenation.')
    with open(pkl_file, "rb") as f:
        adatas = pkl.load(f)
    adata_integrated = adatas[0].concatenate(adatas[1:])

    ## tsne and umap
    sc.pp.neighbors(adata_integrated, n_pcs =40, use_rep = "X_scanorama")
    sc.tl.umap(adata_integrated)
    sc.tl.tsne(adata_integrated, n_pcs = 40, use_rep = "X_scanorama")

    adata_integrated.write(
        output_integrated_h5ad
    )

else:
    combined_adata = sc.read_h5ad(combined_h5ad)

    # Our batch variable is stored under 'sample_name'
    batch_id = 'sample_name'

    adatas = []
    for c in np.unique(combined_adata.obs[batch_id]):
        idx = combined_adata.obs[batch_id] == c
        adatas.append(combined_adata[idx, :])

    # # Run scanorama
    if route == 'correct':
        # -- A integrate + correct
        adatas_integrated = scanorama.correct_scanpy(
            adatas, 
            return_dimred=True,
            verbose=True, 
            dimred=40, # number of PC's to use 
            knn=40, # number of nearest neighbors to use
            )
    elif route == 'integrate':
        # -- B integrate only
        scanorama.integrate_scanpy(adatas) # HZ run1 params, dimred = 40, knn = 40)
        adatas_integrated = adatas.copy()
    else:
        raise ValueError(f'route {route} not supported')

    try:
        output_dir = Path(output_integrated_h5ad).parent
        output_filename = Path(output_integrated_h5ad).stem
        pkl_file = output_dir / f'{output_filename}.pkl'
        with open(pkl_file, "wb") as f: # "wb" because we want to write in binary mode
            pkl.dump(adatas_integrated, f)
        print(f'[SUCCESS] Pickled the scanorama integrated data to {pkl_file} successfully')
    except:
        print('[FAILURE] Failed to pickle the scanorama integrated data')
    adata_integrated = adatas_integrated[0].concatenate(adatas_integrated[1:])
    adata_integrated.write(
        output_integrated_h5ad
    )
