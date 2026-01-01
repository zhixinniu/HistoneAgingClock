############################################################
# Description:
# This script generates age-specific and final consensus peak sets
# for ChIP-seq data.
#
# Major steps include:
# 1) Identifying sub-consensus peak sets at each age/timepoint
#    using MSPC based on narrowPeak files
# 2) Merging age-specific consensus peaks into a final
#    unified consensus peak set for downstream analysis
############################################################


library(jsonlite)
library(dplyr)
library(stringr)



##### Identify sub-consensus peaksets at each timepoint
# Generate age-specific consensus peaks using MSPC

# Create JSON parser configuration for MSPC
# Note: MACS2 outputs -log10(p-value) in column 8 (0-based index = 7)
mspc_parser_config <- list(
  Chr = 0,
  Left = 1,
  Right = 2,
  Name = 3,
  Value = 7,
  Strand = -1,
  Summit = -1,
  Culture = "en-US",
  PValueFormat = 1,
  DefaultValue = 0.0001,
  DropPeakIfInvalidValue = TRUE
)

# Write MSPC parser configuration to JSON file
write(
  toJSON(mspc_parser_config, pretty = T, auto_unbox = T),
  file = "~/zniu_ws/project/epi_clock/encode/HistoneAgingClock/data/peakcalling/parser_config.json"
)


##### Run MSPC for each age group
# NarrowPeak files are grouped by age and processed independently

chip_metadata <- read.table('data/metadata.txt',sep = '\t',header = T)
peak_path <- 'data/peakcalling/'

# Split samples by age and run MSPC within each group
split(chip_metadata, chip_metadata$Age) %>%
  lapply(., function(x){

    peak <- paste0(peak_path, x$IP, '_peaks.narrowPeak', collapse = ' ')

    outpath <- paste0('consensus_peak_for_', unique(x$Age))

    system(paste0(
      '/PATH/TO/MSPC/mspc -i ', peak,
      ' -r bio -c 50% -w 1e-4 -s 1e-8',
      ' -p data/peakcalling/parser_config.json',
      ' -o data/peakcalling/', outpath
    ))
  })



##### Generate final consensus peakset
# Merge all age-specific consensus peak sets into one unified set

# Read all MSPC consensus peak files
all_peaks <- list.dirs("data/peakcalling", recursive = F) %>%
  paste0(., '/ConsensusPeaks.bed') %>%
  lapply(., function(x){

    peak <- read.table(x, sep = '\t', header = T)

    peak$name <- paste0(
      str_split(x, '/', simplify = T)[,3],
      '_peak_', 1:nrow(peak)
    )

    return(peak)
  }) %>%
  do.call(rbind, .) %>%
  arrange(chr, start)

write.table(
  all_peaks,
  'data/peakcalling/all_consensus_peak.bed',
  col.names = F,
  row.names = F,
  quote = F,
  sep = '\t'
)

# Merge overlapping peaks to generate final consensus peak set
system(
  '/PATH/TO/BEDTOOLS/bedtools merge -i data/peakcalling/all_consensus_peak.bed > data/peakcalling/final_consensus_peak_set.bed'
)
