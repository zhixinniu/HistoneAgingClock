# Prerequisites:
# bedtools
# samtools
# trim_galore
# fastqc
# multiqc
# bowtie2
# macs2
# homer



git clone https://github.com/zhixinniu/HistoneAgingClock.git  ## clone full data from HistoneAgingClock repo.
cd HistoneAgingClock/data

# Pre-processing
# QC before processing
mkdir QC_report
fastqc --outdir QC_report --threads 30 *.fq.gz  ## generate QC report for all example data.
multiqc --outdir QC_report QC_report  ## generate QC summary.

# Remove adaptors and low-quality reads
mkdir QC_report_post_trimming
mkdir Clean_reads
trim_galore -q 20 --phred33 --length 15 -e 0.1 --stringency 4 --cores 8 --fastqc --fastqc_args "-t 30 -o QC_report_post_trimming" -o Clean_reads *.fq.gz ## trim and filter reads.
multiqc --outdir QC_report_post_trimming QC_report_post_trimming ## generate QC summary for clean reads.