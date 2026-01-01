############################################################
# Description:
# This script performs post-alignment quality control for
# ChIP-seq data using the ChIPQC package.
#
# It prepares metadata required by ChIPQC and computes
# standard QC metrics based on duplicate-removed BAM files
# and corresponding peak files.
############################################################



##### Post-alignment QC
# Install ChIPQC if not already available
# BiocManager::install("ChIPQC")

library(ChIPQC)
library(dplyr)



##### Prepare metadata
# Construct metadata table required by ChIPQC

# Follow ChIPQC manual:
# https://bioconductor.org/packages/release/bioc/html/ChIPQC.html
#
# Input:
# - metadata describing samples and BAM/peak file locations
#
# Output:
# - formatted metadata data frame for ChIPQC

chip_metadata <- read.table('../HistoneAgingClock_backup/data/metadata.txt',sep = '\t',header = T)

# Build ChIPQC-compatible metadata table
chipqc_metadata <- data.frame(SampleID=chip_metadata$IP,
                              Factor=chip_metadata$Age) %>%
  mutate(Tissue='Spleen',
         Condition='H3K27ac',
         bamReads=normalizePath(paste0('data/Alignment/',SampleID,'.unique_alignment_sorted_rd.bam')),
         Peaks=normalizePath(paste0('data/peakcalling/',SampleID,'_peaks.narrowPeak'))) %>%
  group_by(Factor) %>%
  mutate(Replicate=row_number()) %>%
  ungroup() %>%
  as.data.frame()



##### Perform quality control
# Run ChIPQC to calculate post-alignment QC metrics
chipqc.res <- ChIPQC(chip_metadata, annotation = 'hg19')
