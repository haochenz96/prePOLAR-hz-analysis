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

# ----- load -----
adata = sc.read_h5ad(combined_h5ad)
assert 'X_pca' in adata.obsm.keys(), 'X_pca not found in adata.obsm. Please run PCA first.'

# ----- harmonize -----
Z = harmonize(adata.obsm['X_pca'], adata.obs, batch_key = 'sample_name')
adata.obsm['X_harmony'] = Z
adata.write_h5ad(output_integrated_h5ad)
