# Peak calling
mkdir peakcalling
# For narrow peaks (H3K27ac, H3K4me3,...)
sed '1d' metadata.txt|while read line
do
  line=($line)
  ip=${line[0]}
  input=${line[1]}
  macs2 callpeak -t Alignment/${ip}.unique_alignment_sorted_rd.bam -c Alignment/${input}.unique_alignment_sorted_rd.bam -g hs -f BAM -n ${ip%%_*} --outdir peakcalling
done
# For broad peaks (H3K27me3, H3K36me3,...)
# sed '1d' metadata.txt|while read line
# do
#   line=($line)
#   ip=${line[0]}
#   input=${line[1]}
#   macs2 callpeak -t $ip -c $input -g hs -f BAM --broad -n peakcalling/${ip%%_*} --outdir peakcalling
# done


# Super enhancer identification
# make create tag directory
mkdir tag_directory
ls Alignment/*rd.bam|while read i
do
  makeTagDirectory tag_directory/$(basename ${i%%.*}) $i
done
# identify super enhancer
mkdir super_enhancer
sed '1d' metadata.txt|while read line
do
  line=($line)
  ip=${line[0]}
  input=${line[1]}
  findPeaks tag_directory/$ip -i tag_directory/$input -style super -o super_enhancer/${ip}_SE.txt
  # findPeaks tag_directory/$ip -i tag_directory/$input -style -o auto -typical typicalEnhancer.txt  # add -typical <filename> to output typical (non-super) enhancers
done
