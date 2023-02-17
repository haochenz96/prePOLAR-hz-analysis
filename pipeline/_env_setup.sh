# # ----- it's better to start from mamba-forge instead of Anaconda -----
# cd /home/zhangh5/work/
# wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
# bash Mambaforge-$(uname)-$(uname -m).sh


# # ----- install mamba if needed -----
# conda install mamba -n base -c conda-forge


################
#### SCANPY ####
################


# ----- scanpy and R -----
mamba create -n scanpy -c conda-forge scanpy python-igraph leidenalg r-essentials r-base=4.2.2

# we will directly install from Github
cd /home/zhangh5/work/prePOLAR/analysis/tools
git@github.com:scverse/scanpy.git
cd scanpy
pip install -e '.[dev,doc,test]'

# update gcc for scrublet: 
# mamba install gcc=12.2.0
pip install scrublet

# ######## [OPTIONAL] #########
# # ----- install scib -----
# pip install scib

# # optional packages
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


