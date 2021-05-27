#!/bin/bash --login
#$ -pe smp.pe 8 
#$ -j y
#$ -o /mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/ATAC_ChIP_pipeline/logs

#$ -t 1-1
INDEX=$((SGE_TASK_ID-1))
# CD to directory
cd /mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/ATAC_ChIP_pipeline


# activate all neeeded modules and packages
# source activate personal_software
source activate /mnt/iusers01/jw01/mdefscs4/communal_software/HiC-Pro/conda_hicpro3
# this contains fastp, ontad, bedtools

module load tools/java/1.8.0
module load apps/gcc/R/3.4.2



sleep $(($INDEX*20))
python ./main_ATAC.py -i initial_test_ATAC -s align_qc -s coverage