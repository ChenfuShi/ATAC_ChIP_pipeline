########################################
# utility that contains a class that gets the configuration



########################################

import logging
import datetime

class Config:
    """
    Class containing the parameters 
    """
    
    def __init__(self):

        self.RAW_input_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/reads_here"
        self.Trimmed_dir = "/mnt/iusers01/jw01/mdefscs4/scratch/master_ATAC_ChIP_analyzer/temp_trimming"
        self.aligned_dir = "/mnt/iusers01/jw01/mdefscs4/scratch/master_ATAC_ChIP_analyzer/temp_align"
        self.Reads_quality_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/fastp_qc"
        self.cleaned_alignments_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/clean_alignments"
        self.macs2_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/macs2"
        self.genrich_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/genrich"
        self.coverages_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/coverages"
        self.other_qc_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/qc"

        self.bowtie2_index = "/mnt/jw01-aruk-home01/projects/shared_resources/sequencing/data/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"
        self.genome_fasta = "/mnt/jw01-aruk-home01/projects/shared_resources/sequencing/data/Homo_sapiens/UCSC/hg38/Sequence/hg38/fasta/genome.fa"
        self.picard = "/mnt/jw01-aruk-home01/projects/functional_genomics/bin/picard/picard.jar"
        self.logs_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/ATAC_ChIP_pipeline/logs"

        self._init_logging()

        self.file_to_process = None
        self.analysis_type = None
        
    def _init_logging(self):
        cur_date = datetime.datetime.now()
        
        logging.basicConfig(
            level=logging.INFO,
            format="%(levelname)s - %(message)s",
            handlers=[
                logging.FileHandler("{0}/{1}.log".format(self.logs_dir, f"{cur_date.year}-{cur_date.month}-{cur_date.day}_{cur_date.hour}.{cur_date.minute}.{cur_date.second}"), mode="a"),
                logging.StreamHandler()]) 