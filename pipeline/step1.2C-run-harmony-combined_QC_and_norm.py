# import 
from pathlib import Path
import scanpy as sc
from harmony import harmonize


# installed harmony via:
# `pip install harmony-pytorch``
# or can do 
# `mamba install -c bioconda harmony-pytorch`
# 
# Not sure why but got this error:
"""
"ImportError: /lib64/libstdc++.so.6: version `CXXABI_1.3.9' not found (required by /juno/work/iacobuzc/haochen/mambaforge/envs/scanpy/lib/python3.9/site-packages/pandas/_libs/window/aggregations.cpython-39-x86_64-linux-gnu.so)"
"""
# per: https://stackoverflow.com/questions/48831881/centos-7-libstdc-so-6-version-cxxabi-1-3-9-not-found
# this can be solved with:
# `export LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so`

# ===== IO =====
combined_h5ad = '/home/zhangh5/work/prePOLAR/hz-analysis/outputs/Step1.2-combined_H5/post_cellbender_all_samples-COMBINED_adata.h5ad'

joint_ann_f = '/home/zhangh5/work/prePOLAR/hz-analysis/HZ-prePOLAR-joint_sample_sheet.csv'

output_integrated_h5ad = '/home/zhangh5/work/prePOLAR/hz-analysis/outputs/Step1.2-combined_H5/COMBINED_adata.harmony_integrated.h5ad'
output_dir = Path(output_integrated_h5ad).parent
output_filename = Path(output_integrated_h5ad).stem
pkl_file = output_dir / f'{output_filename}.pkl'

figure_output_dir = '/home/zhangh5/work/prePOLAR/hz-analysis/outputs/Step1.2C-combined_harmony_integrated_H5_annotated_UMAPs'

if not Path(output_integrated_h5ad).exists():
    print(f'[INFO] {output_integrated_h5ad} does not exist. Running integration and saving the file to {output_integrated_h5ad}')
    # ----- load -----
    adata = sc.read_h5ad(combined_h5ad)
    assert 'X_pca' in adata.obsm.keys(), 'X_pca not found in adata.obsm. Please run PCA first.'

    # ----- harmonize -----
    Z = harmonize(adata.obsm['X_pca'], adata.obs, batch_key = 'sample_name')
    adata.obsm['X_harmony'] = Z
    adata.write_h5ad(output_integrated_h5ad)

else:  
    adata_ha = sc.read_h5ad(output_integrated_h5ad)
    ## tsne and umap
    print(f'[INFO] {output_integrated_h5ad} loaded successfully.')
    if not 'neighbors' in adata_ha.uns.keys():
        print(f'[INFO] Computing neighbors')
        sc.pp.neighbors(adata_ha, n_pcs =40, use_rep = "X_harmony")
        sc.tl.umap(adata_ha)
        sc.tl.tsne(adata_ha, n_pcs = 40, use_rep = "X_harmony")
        print(f'[INFO] UMAP and tSNE performed successfully. Saving the file to {output_integrated_h5ad}')
        adata_ha.write(
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
    sc.pl.umap(adata_ha, color="sample_name", show=False, palette = my_palette, ax=axes[0])
    sc.pl.umap(adata_ha, color="leiden", show=False, ax=axes[1])
    fig.set_size_inches(14, 6)
    fig.tight_layout(pad=1) # this adjusts the spacing between the subplots
    fig.savefig(f'{figure_output_dir}/UMAP_sample_name_and_leiden.png')

    # (2) general cell types

    features_of_interest = ['FIBROBLASTS', 'IMMUNE', 'PANCREATIC_SCHWANN_CELLS', 'ENDOTHELIAL', 'ENDOCRINE', 'CFTR', 'DUCTAL', 'scrublet_scores']
    fig = plot_annotated_umaps(adata_ha, features_of_interest, num_cols=4)
    fig.savefig(f'{figure_output_dir}/UMAP_general_cell_types.png')

    # (3) clinical features
    features_of_interest = ['Bx.Timing', 'metastatic_first_line_regimen', 'OS', 'LTS_V_STS', 'BRCA']
    for feature_i in features_of_interest:
        sample_to_feature_map = joint_ann_df[feature_i].to_dict() # this creates a dictionary mapping sample names to the feature of interest. i.e. {sample_name: feature_value}

        # now we add this dictionary to the adata object as an annotation
        adata_ha.obs[feature_i] = adata_ha.obs['sample_name'].map(sample_to_feature_map)
    fig = plot_annotated_umaps(adata_ha, features_of_interest, num_cols=4)
    fig.savefig(f'{figure_output_dir}/UMAP_clinical_features.png')
