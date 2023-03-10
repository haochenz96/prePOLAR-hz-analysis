import scanpy as sc
from math import ceil
import matplotlib.pyplot as plt

def plot_annotated_umaps(adata, features_of_interest, num_cols=4):

    num_features = len(features_of_interest)
    if num_features > num_cols: # if there are more than 4 features, we'll need to make a grid with 2 rows for clarity
        num_rows = ceil(num_features/num_cols)
        fig, axes = plt.subplots(num_rows, num_cols) 
        
        for i, feature_i in enumerate(features_of_interest):
            sc.pl.umap(adata, color=feature_i, show=False, ax=axes[i//num_cols, i%num_cols], cmap="jet")
        fig.set_size_inches(6*num_features//2 + 1*(num_features//2 - 1), 12)

    else:
        fig, axes = plt.subplots(1,num_features)
        for i, feature_i in enumerate(features_of_interest):
            sc.pl.umap(adata, color=feature_i, show=False, ax=axes[i], cmap="jet")
        fig.set_size_inches(6*num_features + 1*(num_features - 1), 6)

    return fig