##### Post-alignment QC
# BiocManager::install("ChIPQC")

library(ChIPQC)
library(dplyr)


##### Prepare metadata
# follow ChIPQC manual: https://bioconductor.org/packages/release/bioc/html/ChIPQC.html
# Input:
# metadata of .rd.bam files
#
# Output:
# QC metrics of samples

chip_metadata <- read.table('../HistoneAgingClock_backup/data/metadata.txt',sep = '\t',header = T)
chipqc_metadata <- data.frame(SampleID=chip_metadata$IP,Factor=chip_metadata$Age)%>%
  mutate(Tissue='Spleen',Condition='H3K27ac',
         bamReads=normalizePath(paste0('data/Alignment/',SampleID,'.unique_alignment_sorted_rd.bam')),
         Peaks=normalizePath(paste0('data/peakcalling/',SampleID,'_peaks.narrowPeak')))%>%
  group_by(Factor)%>%
  mutate(Replicate=row_number())%>%
  ungroup()%>%
  as.data.frame()

##### perform quality control
chipqc.res <- ChIPQC(chip_metadata,annotation = 'hg19')











