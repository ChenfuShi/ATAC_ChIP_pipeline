########################################
# script for running bamcoverage



########################################

import os
import glob
import logging
import subprocess
from steps.helpers import clean_dir

def coverage(Configuration):
    """ 

    """
    logging.info("starting bamcoverage")

    cleaned_align_output_dir = os.path.join(Configuration.cleaned_alignments_dir, Configuration.file_to_process)
    dedup_alignment_file = cleaned_align_output_dir + f"/{Configuration.file_to_process}_align_dedup.bam"    

    coverage_output_dir = os.path.join(Configuration.coverages_dir, Configuration.file_to_process)
    coverage_output_file = coverage_output_dir + f"/{Configuration.file_to_process}_coverage.bw"  
    os.makedirs(coverage_output_dir, exist_ok=True)
    clean_dir(coverage_output_dir)
    # exclude 1804
    # read unmapped (0x4)
    # mate unmapped (0x8)*
    # not primary alignment (0x100)
    # read fails platform/vendor quality checks (0x200)
    # read is PCR or optical duplicate (0x400)

    # include 2
    # read mapped in proper pair (0x2)*

    subprocess.run(f"bamCoverage -p 4 -b {dedup_alignment_file} -of bigwig -o {coverage_output_file} --samFlagInclude 2 --samFlagExclude 1804 --minMappingQuality 30", shell = True)
