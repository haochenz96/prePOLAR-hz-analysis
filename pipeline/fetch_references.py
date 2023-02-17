
# ----- Hwang et al., Nature Genetics (2021) -----
def load_karthik_genemarkers(filename):
    genemarkers = {}

    f = open(filename)
    for line in f:
        tokens = line.strip().split('\t')
        if len(tokens) > 1:
            markers = [gene.upper() for gene in tokens[1:]]
            genemarkers[tokens[0]] = markers
    return genemarkers

def get_broad_celltypes():
    broad_celltypes = {
        "MALIGNANT CELLS" : ['Moffitt_basal','Moffitt_classical','Bailey_squamous','Bailey_progenitor','Collison_QM','Collison_classical'],
        "ACINAR": ["ACINAR"],
        "ENDOCRINE": ["Alpha","Beta","Delta","Gamma","Episilon"],
        "ENDOTHELIAL": ["ENDOTHELIAL"],
        "IMMUNE": ['Pan_Immune','AntigenPresentingCells','Monocytes_1','Monocytes_2','Macrophage','cDC1','cDC2','DC_activated','pDC','Mast','Eosinophils','Neutrophils','M0','M1','M2','Mast_Resting','Mast_activated','CD8_Tcells','CD4_Tcells','NK','CD8_gammadelta','CD8_exhausted','CD4_naive','CD4_memory_resting','CD4_memory_activated','CD4_follicular_helper','CD4_regulatory','NK_resting','NK_activated','B_cell','Plasma','Bcell_naive','Bcell_memory'],
        "FIBROBLASTS": ['PanCAF','iCAF','myCAF','apCAF','CAF','Tuveson_iCAF','Tuveson_mCAF','Neuzillet_CAFa','Neuzillet_CAFb','Neuzillet_CAFc','Neuzillet_CAFd','Davidson_CAF1','Davidson_CAF2','Davidson_CAF3','Pan_Stellate','Quiescent_Stellate','Activated_Stellate','Immune_Stellate'],
        "PANCREATIC_SCHWANN_CELLS": ['PANCREATIC_SCHWANN_CELLS'],
        "DUCTAL":['ductal14', 'ductal2','ductal3','ductal4']
    }
    return broad_celltypes