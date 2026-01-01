############################################################
# Description:
# This script performs raw-data preprocessing and quality
# control for ChIP-seq FASTQ files.
#
# Major steps include:
# 1) Quality control of raw reads (FastQC + MultiQC)
# 2) Adapter trimming and quality filtering (Trim Galore)
# 3) Quality control of cleaned reads
#
# The script assumes all FASTQ files are located in the
# data directory of the HistoneAgingClock repository.
############################################################



# Clone the HistoneAgingClock repository and enter data directory
git clone https://github.com/zhixinniu/HistoneAgingClock.git  ## clone full data from HistoneAgingClock repo.
cd HistoneAgingClock/data



##### Pre-processing
##### QC before processing
# Perform quality control on raw FASTQ files

# Create directory for raw-read QC reports
mkdir QC_report

# Run FastQC on all compressed FASTQ files
fastqc --outdir QC_report --threads 30 *.fq.gz  ## generate QC report for all example data.

# Summarize FastQC results using MultiQC
multiqc --outdir QC_report QC_report  ## generate QC summary.



##### Remove adaptors and low-quality reads
# Trim adapters and low-quality bases from raw reads

# Create output directories for trimmed reads and QC
mkdir QC_report_post_trimming
mkdir Clean_reads

# Run Trim Galore with quality and length filtering
trim_galore -q 20 --phred33 --length 15 -e 0.1 --stringency 4 \
            --cores 8 \
            --fastqc \
            --fastqc_args "-t 30 -o QC_report_post_trimming" \
            -o Clean_reads *.fq.gz  ## trim and filter reads.

# Summarize QC results after trimming
multiqc --outdir QC_report_post_trimming QC_report_post_trimming  ## generate QC summary for clean reads.
