########################################
# script for running bowtie2



########################################

import os
import glob
import logging
import subprocess
from steps.helpers import clean_dir

def align_bowtie(Configuration):
    """ 
    searches for files in the expected trimmed directory and runs bowtie2 in normal configuration for ATAC or ChIP
    there should be only one R1 file and one R2 file because they should have been merged
    """
    logging.info("starting bowtie mapping")

    R1_file = glob.glob(os.path.join(Configuration.Trimmed_dir, Configuration.file_to_process) + "/*R1.fastq.gz")[0]
    R2_file = glob.glob(os.path.join(Configuration.Trimmed_dir, Configuration.file_to_process) + "/*R2.fastq.gz")[0]

    align_output_dir = os.path.join(Configuration.aligned_dir, Configuration.file_to_process)
    os.makedirs(align_output_dir, exist_ok=True)
    clean_dir(align_output_dir)

    # align
    # -k reports maximum 10 alignments and -X reports alignments within 2000bp
    aligner = subprocess.Popen(["bowtie2", "--very-sensitive", "-k", "10", "-X", "2000", "-x", Configuration.bowtie2_index,
                    "-1", R1_file, "-2", R2_file, "-p", "6"], stdout=subprocess.PIPE)
    sorter = subprocess.Popen(["samtools", "sort", "-o", f"{align_output_dir}/{Configuration.file_to_process}_align.bam", "-"],
                    stdin = aligner.stdout)

    sorter.wait()

    # index aligned reads
    subprocess.run(["samtools", "index", f"{align_output_dir}/{Configuration.file_to_process}_align.bam"])


def dedup_QC_alignments(Configuration):
    """
    retrieves the aligned file and runs MarkDuplicates and CollectAlignmentSummaryMetrics
    """
    align_output_dir = os.path.join(Configuration.aligned_dir, Configuration.file_to_process)
    alignment_file = align_output_dir + f"/{Configuration.file_to_process}_align.bam"
    
    cleaned_align_output_dir = os.path.join(Configuration.cleaned_alignments_dir, Configuration.file_to_process)
    os.makedirs(cleaned_align_output_dir, exist_ok=True)
    clean_dir(cleaned_align_output_dir)
    dedup_alignment_file = cleaned_align_output_dir + f"/{Configuration.file_to_process}_align_dedup.bam"

    qc_dir = os.path.join(Configuration.other_qc_dir, Configuration.file_to_process)
    os.makedirs(qc_dir, exist_ok=True)
    clean_dir(qc_dir)
    mark_duplicates_qc = qc_dir + f"/{Configuration.file_to_process}_markdup_qc.txt"
    logging.info("running mark duplicates")
    cmd = f"java -XX:ParallelGCThreads=4 -XX:ParallelCMSThreads=4 -Xmx8G -jar {Configuration.picard} MarkDuplicates QUIET=true REMOVE_DUPLICATES=true CREATE_INDEX=true I={alignment_file} O={dedup_alignment_file} M={mark_duplicates_qc}"
    subprocess.run(cmd, shell = True)

    logging.info("running CollectAlignmentSummaryMetrics")
    metrics_qc = qc_dir + f"/{Configuration.file_to_process}_alignment_metrics_qc.txt"
    cmd = f"java -XX:ParallelGCThreads=4 -XX:ParallelCMSThreads=4 -Xmx8G -jar {Configuration.picard} CollectAlignmentSummaryMetrics R={Configuration.genome_fasta} I={dedup_alignment_file} O={metrics_qc}"
    subprocess.run(cmd, shell = True)

    with open(qc_dir + f"/{Configuration.file_to_process}_idxstats.txt", "w") as idx_output:
        subprocess.run(f"samtools idxstats {dedup_alignment_file}",shell=True,stdout=idx_output)    
        
    with open(qc_dir + f"/{Configuration.file_to_process}_fragment_length_count.txt", "w") as idx_output:
        subprocess.run(f"samtools view {dedup_alignment_file} | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \\t]*//'",shell=True,stdout=idx_output)    