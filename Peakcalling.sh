############################################################
# Description:
# This script performs peak calling and super-enhancer identification
# for ChIP-seq data.
#
# Major steps include:
# 1) Calling narrow or broad peaks using MACS2
# 2) Generating tag directories for downstream analysis
# 3) Identifying super-enhancers using HOMER
#
# The script assumes that alignment and duplicate removal
# have already been completed.
############################################################



##### Peak calling
# Create output directory for peak calling results
mkdir peakcalling

# MACS2 manual:
# https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html
#
# Input:
# - duplicate-removed BAM files (.rd.bam)
#
# Output:
# - narrowPeak or broadPeak files
# - peak summary (.xls)
# - peak summit files (.bed)

# For narrow peaks (e.g. H3K27ac, H3K4me3)
macs2 callpeak -t Alignment/${aligned_IP}.rd.bam \
               -c Alignment/${aligned_Input}.rd.bam \
               -g hs -f BAM -n ${prefix} --outdir peakcalling

# For broad peaks (e.g. H3K27me3, H3K36me3)
macs2 callpeak -t Alignment/${aligned_IP}.rd.bam \
               -c Alignment/${aligned_Input}.rd.bam \
               -g hs -f BAM --broad -n peakcalling/${prefix} --outdir peakcalling





##### Super enhancer identification
# Identify super-enhancers using HOMER

# Create tag directory from ChIP-seq alignments
# Input: duplicate-removed BAM file
# Output: tag directory
mkdir tag_directory
makeTagDirectory tag_directory/${prefix} ${aligned}.rd.bam


# Identify super-enhancers
# Input: tag directories for IP and Input
# Output: super-enhancer list
mkdir super_enhancer
findPeaks tag_directory/${tag_dict_IP} \
          -i tag_directory/${tag_dict_Input} \
          -style super \
          -o super_enhancer/${prefix}_SE.txt  # add -typical <file> to output typical enhancers
