############################################################
# Description:
# This script performs standard ChIP-seq read alignment and
# post-alignment quality control.
#
# Major steps include:
# 1) Aligning single-end ChIP-seq reads to the hg19 genome using Bowtie2
# 2) Generating sorted, duplicate-removed BAM files
# 3) Assessing library complexity (NRF, PBC1, PBC2)
# 4) Calculating strand cross-correlation metrics (NSC, RSC)
#
# The script is designed for single-end ChIP-seq data by default,
# but includes functions supporting paired-end data.
############################################################


# Download and unpack pre-built Bowtie2 index
wget https://genome-idx.s3.amazonaws.com/bt/hg19.zip
unzip hg19.zip


# Store all alignment-related files
mkdir Alignment  # data/Alignment


###### Define function
# Align single-end ChIP-seq reads and generate processed BAM files
function align_chip(){
  INDEX=${1}          # Bowtie2 index prefix
  FASTQ=${2}          # Input FASTQ file
  PREFIX=$(basename ${2})
  PREFIX=${PREFIX%%.*}

  # Align reads to reference genome
  bowtie2 -p 1 -x ${INDEX} -U ${FASTQ} -S Alignment/${PREFIX}.sam  

  # Convert SAM to BAM and keep uniquely mapped reads
  samtools view -@ 1 -h -F 268 -bS Alignment/${PREFIX}.sam > Alignment/${PREFIX}.unique_alignment.bam  

  # Sort BAM file
  samtools sort -o Alignment/${PREFIX}.unique_alignment_sorted.bam -@ 1 Alignment/${PREFIX}.unique_alignment.bam

  # Remove PCR duplicates
  samtools rmdup Alignment/${PREFIX}.unique_alignment_sorted.bam Alignment/${PREFIX}.unique_alignment_sorted_rd.bam

  # Index final BAM file
  samtools index Alignment/${PREFIX}.unique_alignment_sorted_rd.bam
} 



# Calculate NRF and PBC metrics for library complexity
# Specify library type: SE (single-end) or PE (paired-end)
function calc_PBC_NRF(){

  # Check whether library type is provided
  if [ -z $1 ];then
    echo "Please specify the library layout of the input file: SE or PE?"

  # Single-end libraries
  elif [ $1 == "SE" ];then
    sample=$(basename $2)
    sample=${sample%%_*}

    # Convert BAM to BED and count duplicate reads
    bedtools bamtobed -i $2 | \
      awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | \
      grep -v 'chrM' | \
      sort | \
      uniq -c | \
      awk -v name="$sample" '
        BEGIN{mt=0;m0=0;m1=0;m2=0}
        ($1==1){m1=m1+1}
        ($1==2){m2=m2+1}
        {m0=m0+1}
        {mt=mt+$1}
        END{
          m1_m2=-1.0;
          if(m2>0) m1_m2=m1/m2;
          printf "%s\t%f\t%f\t%f\n",name,m0/mt,m1/m0,m1_m2
        }'

  # Paired-end libraries
  elif [ $1 == "PE" ];then
    sample=$(basename $2)
    sample=${sample%%_*}

    # Use BEDPE format for paired-end data
    bedtools bamtobed -bedpe -i $2 | \
      awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | \
      grep -v 'chrM' | \
      sort | \
      uniq -c | \
      awk -v name="$sample" '
        BEGIN{mt=0;m0=0;m1=0;m2=0}
        ($1==1){m1=m1+1}
        ($1==2){m2=m2+1}
        {m0=m0+1}
        {mt=mt+$1}
        END{
          printf "%s\t%.2f\t%.2f\t%.2f\n", name,m0/mt,m1/m0,m1/m2
        }'

  else
    echo "Please specify the library layout of the input file: SE or PE?"
  fi
}  ## Calculate PBC and NRF


# Run alignment for cleaned FASTQ file
align_chip hg19/hg19 ${clean}.fq.gz


##### PBC1, PBC2, and NRF calculation
# Compute library complexity metrics
mkdir QC_report_post_alignment
echo -e "sample\tNRF\tPBC1\tPBC2" > QC_report_post_alignment/NRF_PBC.txt
calc_PBC_NRF SE ${aligned}.bam >> QC_report_post_alignment/NRF_PBC.txt




###### Calculate RSC and NSC by phantompeakqualtools
# Compute strand cross-correlation metrics (NSC, RSC)
Rscript run_spp.R -c=${aligned}.bam -p=30 -out=QC_report_post_alignment/NSC_RSC.txt

# Add header to NSC/RSC output
sed -i '1i Filename\tnumReads\testFragLen\tcorr_estFragLen\tphantomPeak\tcorr_phantomPeak\targmin_corr\tmin_corr\tNSC\tRSC\tQualityTag' QC_report_post_alignment/NSC_RSC.txt



###### Additional QC
# Additional QC metrics (e.g. FRiP) are calculated in Post_alignment_QC.R
