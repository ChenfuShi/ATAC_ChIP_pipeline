#!/bin/bash --login
#$ -pe smp.pe 8 
#$ -j y
#$ -o /mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/ATAC_ChIP_pipeline/logs

#$ -t 1-11
INDEX=$((SGE_TASK_ID-1))
# CD to directory
cd /mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/ATAC_ChIP_pipeline

# Inform the app how many cores we requested for our job. The app can use this many cores.
# The special $NSLOTS keyword is automatically set to the number used on the -pe line above.
export OMP_NUM_THREADS=$NSLOTS

# activate all neeeded modules and packages
# source activate personal_software
source activate /mnt/jw01-aruk-home01/projects/functional_genomics/bin/conda_ATAC_ChIP
# this contains fastp, , bedtools

module load tools/java/1.8.0



sleep $(($INDEX*20))
python ./main_CHIP.py

