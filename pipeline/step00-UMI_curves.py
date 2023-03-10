import glob
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import seaborn.objects as so

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    # color_map="YlGnBu",
    frameon=False,
)

def load_data(filename):
    adata = sc.read_10x_h5(filename)
    ncells, ngenes = adata.shape
    adata.var_names_make_unique()
    return adata

##### ----- IO ----- #####
data_dir = '/home/zhangh5/work/PREPOLAR/PREPOLAR_data/'
raw_h5_fs = sorted(glob.glob(data_dir + '**/raw*.h5', recursive = True))
filtered_h5_fs = sorted(glob.glob(data_dir + '**/filtered*.h5', recursive = True))

# sample_sheet_f = /home/zhangh5/work/PREPOLAR/analysis/Ronan-EK_prePOLAR_sample_sheet.csv
# sample_sheet_df = pd.read_csv(sample_sheet_f, index_col=None)

output_dir = '/home/zhangh5/work/PREPOLAR/analysis/STEP00-UMI_curves/'
Path(output_dir).mkdir(parents=True, exist_ok=True)

#############################
# make output log datasheet
sample_sheet_df = pd.DataFrame()
count = 0
for raw_count_f_i, filtered_h5_fs in zip(raw_h5_fs, filtered_h5_fs):
    
    # infer sample_name
    sample_name = Path(raw_count_f_i).parents[1].name
    sample_name_no_ek = sample_name.lstrip('EK-')

    if not sample_name == Path(filtered_h5_fs).parents[1].name:
        raise ValueError('sample_name does not match!')
        continue
    # if not sample_name_no_ek in sample_sheet_df['Sample ID'].values:
    #     raise ValueError('sample_name does not match!')
    #     exit(1)
    count += 1
    sample_id = 'PP' + str(count).zfill(2)
    print(f'[INFO] processing {sample_name} --> {sample_id}...')
    
    adatas = {}
    adatas['raw'] = load_data(raw_count_f_i)
    adatas['filtered'] = load_data(filtered_h5_fs)
    num_barcodes_over_threshold = {}
    # make output folder
    Path(output_dir + sample_id).mkdir(parents=True, exist_ok=True)


    for type in ['raw', 'filtered']:
        adata = adatas[type]
        # mitochondrial genes
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        # ribosomal genes
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        # hemoglobin genes.
        adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=(5, 10, 20, 50, 100, 200, 500), log1p=True
        )

        total_counts_sorted = adata.obs['total_counts'].sort_values(ascending=False).reset_index()
        num_barcodes_over_threshold[type] = total_counts_sorted.loc[total_counts_sorted['total_counts'] >= 2000].shape[0]

        # ----- plot UMI curve -----
        fig, ax = plt.subplots(figsize=(8, 8))
        p0 = sns.scatterplot(
            x = total_counts_sorted.loc[total_counts_sorted['total_counts'] > 18].index, 
            y = total_counts_sorted.loc[total_counts_sorted['total_counts'] > 18]['total_counts'],
            alpha = 0.8, linewidth=0
            )
        plt.xscale('log')
        plt.yscale('log')
        # plt.vlines(np.exp(kneedle.knee), plt.ylim()[0], plt.ylim()[1], linestyles='dashed', color = 'red')

        p0.set_xticks([500, 5000, 70000], [500, 5000, 70000])
        p0.set_yticks([20,200,2000,20000], [20,200,2000,20000])

        p0.set_xlabel('Droplet ID ranked by count')
        p0.set_ylabel('Total UMI counts')

        fig.savefig(output_dir + sample_id + '/' + f'{sample_id}-{type}-ranked_UMI_curve.png', dpi=200)
    
    # ----- log sample sheet -----
    
    sample_sheet_df.loc[sample_id, ['sample_name','sample_id','raw-num_barcodes_over_2000', 'filtered-num_barcodes_over_2000']] = [sample_name_no_ek, sample_id, num_barcodes_over_threshold['raw'], num_barcodes_over_threshold['filtered']]

sample_sheet_df.to_csv(output_dir + 'sample_sheet_with_UMI_data.csv', index=False)