################
#### SCANPY ####
################


# ----- scanpy and R -----
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
mamba create -n scanpy -c conda-forge scanpy python-igraph leidenalg 


r-essentials r-base=4.2.2

# mamab install r-devtools r-gam r-RColorBrewer r-BiocManager r-Seurat r-sctransform

# might want to install libpng because libpng --> ragg --> devtools/many other R tools above

# # might need a new compiler for r-devtools
# mamba install libgit2

# open R and run _R_env-setup.R


######## [OPTIONAL] #########
# ----- install scib -----
pip install scib

# optional packages
export CFLAGS=-std=c99 # for rpy2 [https://github.com/rpy2/rpy2/issues/963]
pip install rpy2 anndata2ri

# might want to update gcc for building bbknn
mamba install -c conda-forge gcc
pip install bbknn
############################

####################
#### CellBEnder ####
####################

# ----- CellBender -----
mamba create -n torch pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia
# mamba install pytorch torchvision -c pytorch # ***** <-----don't do this as this installs the CPU version of torch
mamba install -c anaconda pytables
cd /home/zhangh5/work/PREPOLAR/analysis/tools
git clone git@github.com:broadinstitute/CellBender.git
pip install -e CellBender


