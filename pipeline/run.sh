#!/bin/bash 
#BSUB -J scanorama_integrate+correct
#BSUB -sla IACOBUZC
#BSUB -n 36                # number of core(tasks) for parallel jobs          
#BSUB -R rusage[mem=16]    # expected resorce consumption for memory
#BSUB -W 2:00            # run time limit (hours)
#BSUB -o /home/zhangh5/work/prePOLAR/hz-analysis/outputs/scanorama_integrate+correct.stdout
#BSUB -eo /home/zhangh5/work/prePOLAR/hz-analysis/outputs/scanorama_integrate+correct.stderr


if [ -f ~/.bashrc ] ; then
    . ~/.bashrc
fi

conda activate scanpy 

# # combined
# python /home/zhangh5/work/prePOLAR/analysis/pipeline/run-QC_and_norm.py

# # individual
# python /home/zhangh5/work/prePOLAR/hz-analysis/pipeline/run-individual_QC_and_norm.py

# scanorama
python /home/zhangh5/work/prePOLAR/hz-analysis/pipeline/step1.2B-run-scanorama-combined_QC_and_norm.py