#!/bin/bash 
#BSUB -J combined-QC_norm_basic-plots
#BSUB -n 36                # number of core(tasks) for parallel jobs          
#BSUB -R rusage[mem=4]    # expected resorce consumption for memory
#BSUB -W 2:00            # run time limit (hours)
#BSUB -o /home/zhangh5/work/prePOLAR/analysis/combined-QC_norm_basic-plots.stdout
#BSUB -eo /home/zhangh5/work/prePOLAR/analysis/combined-QC_norm_basic-plots.stderr


if [ -f ~/.bashrc ] ; then
    . ~/.bashrc
fi

conda activate scanpy 

# combined
python /home/zhangh5/work/prePOLAR/analysis/pipeline/run-QC_and_norm.py

# # individual
# python /home/zhangh5/work/prePOLAR/analysis/pipeline/run-individual_QC_and_norm.py