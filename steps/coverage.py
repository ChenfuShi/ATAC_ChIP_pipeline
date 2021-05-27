########################################
# script for running bamcoverage



########################################

import os
import glob
import logging
import subprocess

def coverage(Configuration):
    """ 
    searches for files in the expected trimmed directory and runs bowtie2 in normal configuration for ATAC or ChIP
    there should be only one R1 file and one R2 file because they should have been merged
    """
    logging.info("starting bamcoverage")

    cleaned_align_output_dir = os.path.join(Configuration.cleaned_alignments_dir, Configuration.file_to_process)
    dedup_alignment_file = cleaned_align_output_dir + f"/{Configuration.file_to_process}_align_dedup.bam"    

    coverage_output_dir = os.path.join(Configuration.coverages_dir, Configuration.file_to_process)
    coverage_output_file = coverage_output_dir + f"/{Configuration.file_to_process}_coverage.bw"  

    subprocess.run(f"bamCoverage -p 4 -b {dedup_alignment_file} -of bigwig -o {coverage_output_file} --samFlagInclude 2 --samFlagExclude 1804 --minMappingQuality 30", shell = True)
