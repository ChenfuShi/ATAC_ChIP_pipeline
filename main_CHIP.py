########################################
# main script that calls all functions

# naming system:
# folder name needs to be sampleName_PROTOCOL
# within the folder name needs to be id+lane_R1.fastq.gz and R2

# protocols needs to be either CHIP or ATAC

# designed to work specifically for our library prepration methods
# only usable with paired end reads

########################################

from configuration import Config
import os
import glob
from random import random
from time import sleep
import argparse
import logging
from steps import trimming, align, coverage, genrich, macs2


if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Wrapper function for all steps of ChIP-seq analysis')

    parser.add_argument("-i",'--input', dest='infile', action='store', required=False,
                        help='input folder to force. Will overwrite all ouputs')
    parser.add_argument("-s",'--steps', dest='step', action='append', required=False,
                        help='chose steps instead of running everything')
    parser.add_argument("-c",'--no-control', dest='no_control', action='store_true', required=False, default=False,
                        help='disable the use of input background for chip')                        

    # parse arguments
    args = parser.parse_args()

    # set up configuration object for all steps. this sets up logging as well
    Configuration = Config()

    if args.infile == None:
        # First check if there are any ChIP input samples to run
        all_raw_inputs = [os.path.basename(x) for x in glob.glob(Configuration.RAW_input_dir + "/*CHIP_INPUT")]
        
        all_processed_inputs = [os.path.basename(x) for x in glob.glob(Configuration.cleaned_alignments_dir + "/*CHIP_INPUT")]
        # chose the first one of the ones that are still not processed and run 
        for i in all_raw_inputs:
            if i not in all_processed_inputs:
                os.makedirs(os.path.join(Configuration.cleaned_alignments_dir,i),exist_ok=True)
                Configuration.file_to_process = i
                Configuration.analysis_type = "CHIP_INPUT"
                break

        if Configuration.file_to_process == None:
            # if it reaches here it means that we processed all the input files and we can run the samples
            Configuration.analysis_type = "CHIP"
            all_raws_present = [os.path.basename(x) for x in glob.glob(Configuration.RAW_input_dir + "/*CHIP")]

            all_processed = [os.path.basename(x) for x in glob.glob(Configuration.cleaned_alignments_dir + "/*CHIP")]
            for i in all_raws_present:
                if i not in all_processed:
                    os.makedirs(os.path.join(Configuration.cleaned_alignments_dir,i),exist_ok=True)
                    Configuration.file_to_process = i
                    break

        if Configuration.file_to_process == None:
            logging.error("There were no new files to process")
            raise Exception
        
    else:
        Configuration.file_to_process = args.infile
        os.makedirs(os.path.join(Configuration.cleaned_alignments_dir,Configuration.file_to_process),exist_ok=True)
        if args.infile.endswith('CHIP_INPUT'):
            Configuration.analysis_type = "CHIP_INPUT"
        else:
            Configuration.analysis_type = "CHIP"

    
    logging.info(f"This script will run the file : {Configuration.file_to_process}")

    if Configuration.analysis_type == "CHIP" and args.no_control == False:
        # identify the control input for the sample if any
        controls = [os.path.basename(x) for x in glob.glob(os.path.join(Configuration.RAW_input_dir, Configuration.file_to_process) + "/*CHIP_INPUT")]
        if len(controls) == 1:
            Configuration.input_background = controls[0]
            logging.info(f"Identified background control file {Configuration.input_background}")

    if args.step == None:
        # call trimming
        trimming.run_fastp(Configuration)

        # run bowtie2 alignment
        align.align_bowtie(Configuration)

        # run QC and dedup
        align.dedup_QC_alignments(Configuration)

        # run coverage
        coverage.coverage(Configuration)

        # create file for genrich
        genrich.create_bam_for_genrich(Configuration)

        # create file for macs2
        macs2.create_bam_for_macs2_ATAC(Configuration)

        if Configuration.analysis_type == "CHIP":
            # run Genrich
            genrich.run_genrich_CHIP(Configuration)

            # run macs2
            macs2.run_macs2_CHIP(Configuration)

    else:
        if "trimming" in args.step:
            trimming.run_fastp(Configuration)
        if "align" in args.step:
            align.align_bowtie(Configuration)
        if "align_qc" in args.step:
            align.dedup_QC_alignments(Configuration)
        if "coverage" in args.step:
            coverage.coverage(Configuration)
        if "genrich" in args.step:
            genrich.create_bam_for_genrich(Configuration)
            genrich.run_genrich_CHIP(Configuration)
        if "macs2" in args.step:
            macs2.create_bam_for_macs2_ATAC(Configuration)
            macs2.run_macs2_CHIP(Configuration)