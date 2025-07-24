##### Download Bowtie2 index for UCSC hg19
wget https://genome-idx.s3.amazonaws.com/bt/hg19.zip
unzip hg19.zip



##### Create the output directory of alignment files
mkdir Alignment  # data/Alignment



###### Define function
# alignment of single-end data
function align_chip(){
  INDEX=${1}
  FASTQ=${2}
  PREFIX=$(basename ${2})
  PREFIX=${PREFIX%%_*}
  # change "-p" to use more cores as desire.
  bowtie2 -p 1 -x ${INDEX} -U ${FASTQ} -S Alignment/${PREFIX}.sam  
  samtools view -@ 1 -h -F 268 -bS Alignment/${PREFIX}.sam > Alignment/${PREFIX}.unique_alignment.bam  ## change "-@" to use more cores as desire.
  samtools sort -o Alignment/${PREFIX}.unique_alignment_sorted.bam -@ 1 Alignment/${PREFIX}.unique_alignment.bam
  samtools rmdup Alignment/${PREFIX}.unique_alignment_sorted.bam Alignment/${PREFIX}.unique_alignment_sorted_rd.bam
  samtools index Alignment/${PREFIX}.unique_alignment_sorted_rd.bam
} 

# calculate NRF and PBC
# specify "SE" or "PE" based on your data type at first place.
function calc_PBC_NRF(){
  if [ -z $1 ];then
    echo "Please specify the library layout of the input file: SE or PE?"
  elif [ $1 == "SE" ];then
    sample=$(basename $2)
    sample=${sample%%_*}
    bedtools bamtobed -i $2|\
      awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' |  \
      grep -v 'chrM' | \
      sort |  \
      uniq -c |  \
      awk -v name="$sample" 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} \
      END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%s\t%f\t%f\t%f\n",name,m0/mt,m1/m0,m1_m2}'
  elif [ $1 == "PE" ];then
    sample=$(basename $2)
    sample=${sample%%_*}
    bedtools bamtobed -bedpe -i $2|\
      awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' |\
      grep -v 'chrM' | \
      sort | \
      uniq -c | \
      awk -v name="$sample" 'BEGIN{mt=0;m0=0;m1=0;m2=0}($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} \
  		END{printf "%s\t%.2f\t%.2f\t%.2f\n", name,m0/mt,m1/m0,m1/m2}'
  else
    echo "Please specify the library layout of the input file: SE or PE?"
  fi
}  ##  Calulate PBC and NRF



##### Alignment
# define the path to bowtie2 index, FASTQ file, and the prefix of output files.
# 
# Input: 
# cleaned fastq files
# 
# Output: 
# alignment files in sam and bam format (.sam and .bam)
# sorted alignment bam file (.sorted.bam)
# duplicates removed bam file (.rd.bam)
# index of bam file

align_chip hg19/hg19 ${clean}.fq.gz




##### PBC1, PBC2, and NRF calculation
# Input:
# .rd.bam
# 
# Output:
# calculated NRF and PBC file (NRF_PBC.txt)

mkdir QC_report_post_alignment
echo -e "sample\tNRF\tPBC1\tPBC2" > QC_report_post_alignment/NRF_PBC.txt
calc_PBC_NRF SE ${aligned}.bam >> QC_report_post_alignment/NRF_PBC.txt




###### Calculate RSC and NSC by phantompeakqualtools
# find run_spp.R in https://github.com/kundajelab/phantompeakqualtools.git
# Input:
# .sorted.bam file
# 
# Output:
# calculated NSC and RSC file (NSC_RSC.txt)

Rscript run_spp.R -c=${aligned}.bam -p=30 -out=QC_report_post_alignment/NSC_RSC.txt
sed -i '1i Filename\tnumReads\testFragLen\tcorr_estFragLen\tphantomPeak\tcorr_phantomPeak\targmin_corr\tmin_corr\tNSC\tRSC\tQualityTag' QC_report_post_alignment/NSC_RSC.txt


###### Please see the additional creteria of quality control, e.g. FRiP, in the R script "Post_alignment_QC.R"
