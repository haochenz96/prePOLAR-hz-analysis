# # ----- it's better to start from mamba-forge instead of Anaconda -----
# cd /home/zhangh5/work/
# wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
# bash Mambaforge-$(uname)-$(uname -m).sh


################
#### SCANPY ####
################


# ----- scanpy -----
mamba create -n scanpy -c conda-forge scanpy python-igraph leidenalg r-essentials r-base=4.2.2
# # or we can directly install scanpy from Github
# cd /home/zhangh5/work/prePOLAR/analysis/tools
# git@github.com:scverse/scanpy.git
# cd scanpy
# pip install -e '.[dev,doc,test]'

# might need to update gcc for scrublet: 
# mamba install gcc=12.2.0
# pip install scrublet
# or if compilation does not work, use mamba
# # mamba install -c bioconda scrublet

# --- batch effect correction ---
# (1) scanorama
pip install scanorama
mamba install -c conda-forge python-annoy

# (2) harmony
pip install harmony-pytorch

# ######## [OPTIONAL] #########

# # ----- Jupyter -----
# mamba install ipykernel ipython

# # ----- R -----
# mamba install r-essentials r-base=4.2.2 r-devtools r-gam r-RColorBrewer r-BiocManager r-Seurat r-sctransform

# might want to install libpng because libpng --> ragg --> devtools/many other R tools above

# # might need a new compiler for r-devtools
# mamba install libgit2

# open R and run _R_env-setup.R

# # ----- install scib -----
# pip install scib

# # ----- rpy2 -----
# export CFLAGS=-std=c99 # for rpy2 [https://github.com/rpy2/rpy2/issues/963]
# pip install rpy2

# # might want to update gcc for building bbknn
# mamba install -c conda-forge gcc
# pip install bbknn
############################

# ####################
# #### CellBEnder ####
# ####################

# # ----- CellBender -----
# mamba create -n torch pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia
# # mamba install pytorch torchvision -c pytorch # ***** <-----don't do this as this installs the CPU version of torch
# mamba install -c anaconda pytables
# cd /home/zhangh5/work/prePOLAR/analysis/tools
# git clone git@github.com:broadinstitute/CellBender.git
# pip install -e CellBender


