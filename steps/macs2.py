########################################
# script for running macs2



########################################

import os
import glob
import logging
import subprocess
from steps.helpers import clean_dir

def create_bam_for_macs2_ATAC(Configuration):
    """
    filters the data for macs2, this was done automatically in genrich but macs2 is dumb
    this file could also be useful for further analysis as it only contains proper reads
    """
    logging.info("creating cleaned bam file for macs2")
    cleaned_align_output_dir = os.path.join(Configuration.cleaned_alignments_dir, Configuration.file_to_process)
    dedup_alignment_file = cleaned_align_output_dir + f"/{Configuration.file_to_process}_align_dedup.bam"

    filtered_align_file = cleaned_align_output_dir + f"/{Configuration.file_to_process}_align_filtered_macs2.bam"

    # sort by name
    cmd = f"samtools view -h {dedup_alignment_file} | grep -v chrM | samtools view -h -q 30 - | samtools view -h -b -F 1804 -f 2 | samtools sort -O bam -o {filtered_align_file}"
    subprocess.run(cmd, shell=True)
    subprocess.run(["samtools", "index", filtered_align_file])

def run_macs2_ATAC(Configuration):
    """
    run macs2 application
    """
    logging.info("running macs2")
    cleaned_align_output_dir = os.path.join(Configuration.cleaned_alignments_dir, Configuration.file_to_process)
    filtered_align_file = cleaned_align_output_dir + f"/{Configuration.file_to_process}_align_filtered_macs2.bam"

    macs2_output_dir = os.path.join(Configuration.macs2_dir, Configuration.file_to_process)
    os.makedirs(macs2_output_dir, exist_ok = True)
    clean_dir(macs2_output_dir)
    subprocess.run(["macs2", "callpeak", "-f", "BAMPE", "-g", "hs", "--keep-dup", "all",
        "-n", Configuration.file_to_process, "-t", filtered_align_file, "--outdir", macs2_output_dir])

def run_macs2_CHIP(Configuration):
    """
    run macs2 application
    """
    logging.info("running macs2")
    cleaned_align_output_dir = os.path.join(Configuration.cleaned_alignments_dir, Configuration.file_to_process)
    filtered_align_file = cleaned_align_output_dir + f"/{Configuration.file_to_process}_align_filtered_macs2.bam"

    macs2_output_dir = os.path.join(Configuration.macs2_dir, Configuration.file_to_process)
    os.makedirs(macs2_output_dir, exist_ok = True)
    clean_dir(macs2_output_dir)

    if Configuration.input_background is not None:
        background_dir = os.path.join(Configuration.cleaned_alignments_dir, Configuration.input_background)
        background_bam = f"{background_dir}/{Configuration.input_background}_align_filtered_macs2.bam"    

        subprocess.run(["macs2", "callpeak", "-f", "BAMPE", "-g", "hs", "--keep-dup", "all",
            "-n", Configuration.file_to_process, "-t", filtered_align_file, "-c", background_bam,"--outdir", macs2_output_dir])
    else:
        subprocess.run(["macs2", "callpeak", "-f", "BAMPE", "-g", "hs", "--keep-dup", "all",
            "-n", Configuration.file_to_process, "-t", filtered_align_file, "--outdir", macs2_output_dir])