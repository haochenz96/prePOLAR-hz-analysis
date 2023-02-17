import scanpy as sc

def load_h5_data(filename):
    adata = sc.read_10x_h5(filename)
    ncells, ngenes = adata.shape
    adata.var_names_make_unique()
    return adata