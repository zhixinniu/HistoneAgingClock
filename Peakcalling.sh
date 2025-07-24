##### Peak calling
mkdir peakcalling

# MACS2 manual: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html
# Input:
# .rd.bam files
# Ouput:
# .narrow(broad)Peak
# .xls for identified peaks
# .bed for peak summits

# For narrow peaks (H3K27ac, H3K4me3,...)
macs2 callpeak -t Alignment/${aligned_IP}.rd.bam -c Alignment/${aligned_Input}.rd.bam -g hs -f BAM -n ${prefix} --outdir peakcalling

# For broad peaks (H3K27me3, H3K36me3,...)
macs2 callpeak -t Alignment/${aligned_IP}.rd.bam -c Alignment/${aligned_Input}.rd.bam -g hs -f BAM --broad -n peakcalling/${prefix} --outdir peakcalling





##### Super enhancer identification

# make create tag directory
# Input:
# .rd.bam
# Output:
# tag directory for super-enhancer identification

mkdir tag_directory
makeTagDirectory tag_directory/${prefix} ${aligned}.rd.bam


# identify super enhancer
# Input:
# tag directory
# Output:
# .txt file for identified super enhancer
mkdir super_enhancer
findPeaks tag_directory/${tag_dict_IP} -i tag_directory/${tag_dict_Input} -style super -o super_enhancer/${prefix}_SE.txt  # add -typical <filename> to output typical (non-super) enhancers


