# ATAC_ChIP_pipeline
pipeline for processing of ATAC-seq and ChIP-seq data

Currently setup for ATAC-seq data only. 
Submit scripts with the run_atac.sh file

### How to use:

Generally there are a lot of parameters to set in the configuration.py file. These should work but you might need to reconfigure the location of the scratch folders for your account.

The way this software works is that it checks for folders in the reads_here directory and then checks which folders exist in the clean_alignments directory. It then choses one of the folders that are discrepant and analyses that. Once it choses which file to run it will create a folder in the clean_alignments directory to hold that file for the other parallel scripts. The delay of 20 seconds is there to make sure two scripts don't collide into eachother. If you want to reset the run you need to remove the folder from the clean_alignments directory or you can override the automatic process by inputing the folder as -i SAMPLE_PROTOCOL. You can also override steps to run by adding -s STEP_1 -s STEP_2 ....

#### Processing data

1. Create symbolic links to the input data folder (reads_here). Each sample needs to be within one folder, will all the reads for all the lanes in there. The folder name needs to end in one of the following: ATAC or CHIP. Set this based on the library you are processing. It is very important that the name of the files are the same for the paired end execpt the R1.fastq.gz and R2.fastq.gz. So your files should look like SAMPLENAME_lane_R1.fastq.gz and SAMPLENAME_lane_R2.fastq.gz.

2. Set in run_atac.sh the number of samples you are going to process. To do so set #$ -t 1-N with N being the number of samples you have.

3. Submit the script to the CSF.
   
#### QC metrics
QC metrics for sequencing (e.g. trimming) reads will be located in fastp_qc. 
QC metrics for everything else will be located in the qc folder. These include peak calling metrics, alignment metrics, and fraction of reads in TSS and peaks.