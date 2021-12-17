########################################
# script for running genrich



########################################

import os
import glob
import logging
import subprocess
from steps.helpers import clean_dir

def create_bam_for_genrich(Configuration):
    """
    this needs to take the original align file (not filtered/deduplicated) because genrich does it on its own.
    also needs to be sorted by readname for some stupid optical deduplication stuff
    """
    logging.info("creating sorted bam by name for genrich")
    align_output_dir = os.path.join(Configuration.aligned_dir, Configuration.file_to_process)
    aligned_bam = f"{align_output_dir}/{Configuration.file_to_process}_align.bam"

    # sort by name
    subprocess.run(["samtools", "sort", "-n", "-o", f"{align_output_dir}/{Configuration.file_to_process}_align_genrich.bam", aligned_bam])

def run_genric_ATAC(Configuration):
    """
    run genrich application with ATAC-seq settings
    """
    logging.info("running genrich")
    align_output_dir = os.path.join(Configuration.aligned_dir, Configuration.file_to_process)

    genrich_output_dir = os.path.join(Configuration.genrich_dir, Configuration.file_to_process)
    os.makedirs(genrich_output_dir, exist_ok = True)
    clean_dir(genrich_output_dir)
    genrich_output = genrich_output_dir + f"/{Configuration.file_to_process}_genrich.bed"
    aligned_bam_genrich = f"{align_output_dir}/{Configuration.file_to_process}_align_genrich.bam"

    subprocess.run(["Genrich", "-t", aligned_bam_genrich, "-o", genrich_output, "-j", "-y", "-r", "-e", "chrM"])


def run_genrich_CHIP(Configuration):
    """
    run genrich application with ChIP-seq settings
    supporting input control files
    """
    logging.info("running genrich")
    align_output_dir = os.path.join(Configuration.aligned_dir, Configuration.file_to_process)
    genrich_output_dir = os.path.join(Configuration.genrich_dir, Configuration.file_to_process)
    os.makedirs(genrich_output_dir, exist_ok = True)
    clean_dir(genrich_output_dir)
    genrich_output = genrich_output_dir + f"/{Configuration.file_to_process}_genrich.bed"
    aligned_bam_genrich = f"{align_output_dir}/{Configuration.file_to_process}_align_genrich.bam"

    if Configuration.input_background is not None:
        background_dir = os.path.join(Configuration.aligned_dir, Configuration.input_background)
        background_bam = f"{background_dir}/{Configuration.input_background}_align_genrich.bam"

        subprocess.run(["Genrich", "-t", aligned_bam_genrich, "-c", background_bam, "-o", genrich_output, "-r", "-e", "chrM"])
    else:
        subprocess.run(["Genrich", "-t", aligned_bam_genrich, "-o", genrich_output, "-r", "-e", "chrM"])
