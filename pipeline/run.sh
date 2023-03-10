#!/bin/bash 
#BSUB -J harmony_integrate
#BSUB -n 8                # number of core(tasks) for parallel jobs          
#BSUB -R rusage[mem=16]    # expected resorce consumption for memory
#BSUB -q gpuqueue
#BSUB -gpu num=1:gmem=16
#BSUB -W 2:00            # run time limit (hours)
#BSUB -o /home/zhangh5/work/prePOLAR/hz-analysis/outputs/harmony_integrate.stdout
#BSUB -eo /home/zhangh5/work/prePOLAR/hz-analysis/outputs/harmony_integrate.stderr


if [ -f ~/.bashrc ] ; then
    . ~/.bashrc
fi

conda activate scanpy 

# # combined
# python /home/zhangh5/work/prePOLAR/analysis/pipeline/run-QC_and_norm.py

# # individual
# python /home/zhangh5/work/prePOLAR/hz-analysis/pipeline/run-individual_QC_and_norm.py

# # scanorama
# python /home/zhangh5/work/prePOLAR/hz-analysis/pipeline/step1.2B-run-scanorama-combined_QC_and_norm.py

# harmony
python /home/zhangh5/work/prePOLAR/hz-analysis/pipeline/step1.2C-run-harmony-combined_QC_and_norm.py