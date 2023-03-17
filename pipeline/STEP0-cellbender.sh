#!/bin/bash 
#BSUB -J CellBender-PP22
#BSUB -n 4                # number of core(tasks) for parallel jobs          
#BSUB -R rusage[mem=16]    # expected resorce consumption for memory
#BSUB -q gpuqueue
#BSUB -gpu num=1:gmem=16
#BSUB -W 4:00            # run time limit (hours)
#BSUB -o /home/zhangh5/work/prePOLAR/analysis/CellBender-PP22.stdout
#BSUB -eo /home/zhangh5/work/prePOLAR/analysis/CellBender-PP22.stderr

if [ -f ~/.bashrc ] ; then
    . ~/.bashrc
fi

conda activate torch

DATA_DIR=/home/zhangh5/work/prePOLAR/raw_CR_outputs
OUTPUT_DIR=/home/zhangh5/work/prePOLAR/analysis/outputs/CellBender
mkdir -p ${OUTPUT_DIR}

# # ----- check GPU load if needed -----
# lsload -w -gpuload

# for H5_i in $(find $DATA_DIR -name "*raw_feature_bc_matrix.h5")
# do 

----- single-sample run -----
H5_i=/home/zhangh5/work/prePOLAR/raw_CR_outputs/EK-1898_12-045_NT_9_27_18_Baseline.cr.raw_feature_bc_matrix.h5
SAMPLE_ID=PP22

# SAMPLE_NAME=$(echo $H5_i | cut -d / -f 7) # e.g. EK-1663_PP01_IR0007
SAMPLE_NAME=$(basename $H5_i | cut -d . -f 1) # e.g. EK-1663_PP01_IR0007
F_NAME=${SAMPLE_ID}.cr.raw_feature_bc_matrix.h5
# rsync -avzu $H5_i zhangh5@lilac.mskcc.org:/home/zhangh5/work/prePOLAR/raw_CR_outputs/${F_NAME}


mkdir -p ${OUTPUT_DIR}/${SAMPLE_ID}
# run CellBender with default params
cellbender remove-background \
    --input ${H5_i} \
    --output ${OUTPUT_DIR}/${SAMPLE_ID}/${SAMPLE_ID}.cellbender_filtered.h5 \
    --cuda \
    --expected-cells 2000 \
    --total-droplets-included 7000 \
    --fpr 0.01 \
    --epochs 150
# done




