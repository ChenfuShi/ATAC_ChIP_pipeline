########################################
# script for running qc for the individual sample



########################################

import os
import glob
import logging
import subprocess
from steps.helpers import clean_dir
import numpy as np
import pandas as pd
import pybedtools as pbed
import math
import statistics
import concurrent.futures
os.makedirs("/mnt/iusers01/jw01/mdefscs4/scratch/temp_pybedtools/", exist_ok = True)
pbed.helpers.set_tempdir("/mnt/iusers01/jw01/mdefscs4/scratch/temp_pybedtools/")
bed_genome_file = "/mnt/iusers01/jw01/mdefscs4/hg38.genome"

def run_qc(Configuration):
    tss_sites = "/mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/NEW_references/genes/gencode.v29.TSS_sites_protein_coding_sorted.bed"

    tss_bed = pbed.BedTool(tss_sites)
    tss_bed = tss_bed.slop(b=2000,g=bed_genome_file)

    def wccount(filename):
        out = subprocess.Popen(['wc', '-l', filename],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT
                            ).communicate()[0]
        return int(out.partition(b' ')[0])

    def samtools_count(filename):
        out = subprocess.Popen(['samtools','view','-c', filename],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT
                            ).communicate()[0]
        return int(out)

    def calculate_FRiP(tss_bed, genrich, macs, bam_file):
        bam_bed = pbed.BedTool(bam_file)
        reads_in_tss = bam_bed.intersect(tss_bed,bed=True,u=True).fn
        reads_in_macs = bam_bed.intersect(macs,bed=True,u=True).fn
        reads_in_genrich = bam_bed.intersect(genrich,bed=True,u=True).fn
        reads_in_tss = wccount(reads_in_tss)
        reads_in_macs = wccount(reads_in_macs)
        reads_in_genrich = wccount(reads_in_genrich)
        total_reads = samtools_count(bam_bed.fn)
        return reads_in_tss/total_reads, reads_in_macs/total_reads, reads_in_genrich/total_reads

    def call_frip(sample):
        return calculate_FRiP(tss_bed, 
                            pbed.BedTool(os.path.join(Configuration.genrich_dir,sample,sample + "_genrich.bed")), 
                            pbed.BedTool(os.path.join(Configuration.macs2_dir,sample,sample + "_peaks.narrowPeak")), 
                            os.path.join(Configuration.cleaned_alignments_dir,sample,sample + "_align_dedup.bam"))



    res = call_frip(Configuration.file_to_process)
    print(Configuration.file_to_process, f"reads in tss: {res[0]:.2f}, reads in macs peaks {res[1]:.2f}, reads in genrich peaks {res[2]:.2f}")
    quality_dir = os.path.join(Configuration.other_qc_dir,Configuration.file_to_process)
    with open(os.path.join(quality_dir, Configuration.file_to_process + "_FRIP_qc.txt"), "w") as output_qc:
        output_qc.write(f"{Configuration.file_to_process}, reads in tss: {res[0]:.2f}, reads in macs peaks {res[1]:.2f}, reads in genrich peaks {res[2]:.2f}")

