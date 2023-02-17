import scanpy as sc

import seaborn as sns
import matpolib.pyplot as plt



def plot_UMI_curve(adata, sample_name, total_counts_threshold=None):
    sc.pp.calculate_qc_metrics(
        adata, inplace=True, log1p=True
    )
    fig, ax = plt.subplots(figsize=(8, 8))


    total_counts_sorted = adata.obs['total_counts'].sort_values(ascending=False).reset_index()

    if total_counts_threshold is not None:
        total_counts_sorted = total_counts_sorted.loc[total_counts_sorted['total_counts'] > total_counts_threshold]

    p0 = sns.scatterplot(
        x = total_counts_sorted.index, 
        y = total_counts_sorted['total_counts'],
        alpha = 0.8, linewidth=0
        )
    plt.xscale('log')
    plt.yscale('log')

    p0.set_xticks([500, 1500, 5000, 15000, 70000], [500, 1500, 5000, 15000, 70000])
    p0.set_yticks([20,200,2000,20000], [20,200,2000,20000])

    p0.set_xlabel('Droplet ID ranked by count')
    p0.set_ylabel('Total UMI counts')

    # plt.vlines(2000, plt.ylim()[0], plt.ylim()[1], linestyles='dashed', color = 'green', alpha = 0.8, label = 'expected_cell_number')

    # plt.vlines(15000, plt.ylim()[0], plt.ylim()[1], linestyles='dashed', color = 'red', alpha = 0.8, label = 'total_droplets_included')

    # plt.legend(bbox_to_anchor = (1.0, 1), loc = 'upper left')
    return fig
    # fig.savefig(f'{sample_name}-total_counts_sorted.png', dpi=200)