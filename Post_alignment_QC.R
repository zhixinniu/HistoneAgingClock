library(ChIPQC)


chip_metadata <- read.table('data/metadata.txt',sep = '\t',header = T)
# Generate ChIPQC metadata
chip_metadata <- data.frame(SampleID=chip_metadata$IP,Factor=chip_metadata$Age)%>%
  mutate(Tissue='Spleen',Condition='H3K27ac',
         bamReads=normalizePath(paste0('data/Alignment/',SampleID,'.unique_alignment_sorted_rd.bam')),
         Peaks=normalizePath(paste0('data/peakcalling/',SampleID,'_peaks.narrowPeak')))%>%
  group_by(Factor)%>%
  mutate(Replicate=row_number())%>%
  ungroup()%>%
  as.data.frame()

chipqc.res <- ChIPQC(chip_metadata,annotation = 'hg19')











