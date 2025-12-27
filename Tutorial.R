# Step-by-step Tutorial

## Part I （Steps 1–3）

### Overview

This tutorial describes a modular workflow for constructing histone modification–based biological age clocks from ChIP-seq data. The pipeline is designed to be executed step by step, allowing users to inspect quality control metrics and intermediate results before proceeding to downstream analyses.

The first three steps cover:

1.  Raw read quality control and preprocessing

2.  Alignment and post-alignment quality assessment

3.  Post-alignment QC summarization in R

### Step 1. Raw read quality control and preprocessing

**Script:** `Preprocessing.sh`

#### Purpose

This step performs initial quality control of raw ChIP-seq FASTQ files, followed by adapter trimming and removal of low-quality reads.

#### Input

-   Raw ChIP-seq reads in gzipped or raw FASTQ format (`*.fq.gz, *.fq`)

-   Single-end or paired-end reads are supported (example shown for single-end)

#### Software requirements

-   FastQC

-   MultiQC

-   Trim Galore

-   bedtools, samtools, bowtie2 (required in later steps)

#### Procedure

First, clone the repository and navigate to the data directory:

```         
git clone https://github.com/zhixinniu/HistoneAgingClock.git cd HistoneAgingClock/data 
```

Run quality control on raw reads:

```         
mkdir QC_report 
fastqc --outdir QC_report --threads 30 *.fq.gz multiqc --outdir QC_report QC_report 
```

Trim adapters and low-quality bases:

```         
mkdir QC_report_post_trimming 
mkdir Clean_reads  

trim_galore -q 20 --phred33 --length 15 -e 0.1 --stringency 4 --cores 8 --fastqc --fastqc_args "-t 30 -o QC_report_post_trimming" -o Clean_reads *.fq.gz

multiqc --outdir QC_report_post_trimming QC_report_post_trimming 
```

#### Output

-   `QC_report/`: FastQC and MultiQC reports for raw reads

-   `Clean_reads/`: trimmed and filtered FASTQ files

-   `QC_report_post_trimming/`: QC reports after trimming

These cleaned FASTQ files are used as input for alignment in Step 2.

### Step 2. Alignment and post-alignment quality control

**Script:** `Alignment.sh`

#### Purpose

This step maps cleaned ChIP-seq reads to the reference genome, removes PCR duplicates, and computes post-alignment QC metrics, including NRF, PBC, NSC, and RSC.

#### Input

-   Cleaned FASTQ files from Step 1 (`Clean_reads/*.fq.gz`)

-   Bowtie2 index of the reference genome (example: hg19)

#### Reference genome preparation

Download and unzip the Bowtie2 index:

```         
wget https://genome-idx.s3.amazonaws.com/bt/hg19.zip unzip hg19.zip 
```

#### Alignment

Reads are aligned using Bowtie2, followed by filtering, sorting, duplicate removal, and indexing.\
Example command (single-end data):

```         
align_chip hg19/hg19 ${clean}.fq.gz 
```

> **Note:**
>
> -   `${clean}` refers to the prefix of the cleaned FASTQ file
>
> -   Users are expected to run this command iteratively for each sample
>
> -   CPU core usage (`-p`, `-@`) can be adjusted depending on available resources

#### Output

Alignment results are written to the `Alignment/` directory:

-   `.sam`: raw alignment

-   `.unique_alignment.bam`: filtered BAM

-   `.unique_alignment_sorted_rd.bam`: sorted, duplicate-removed BAM

-   `.bam.bai`: BAM index

#### Calculation of NRF and PBC

NRF and PBC metrics are computed from duplicate-removed BAM files.

```         
mkdir QC_report_post_alignment 
echo -e "sample\tNRF\tPBC1\tPBC2" > QC_report_post_alignment/NRF_PBC.txt 
calc_PBC_NRF SE ${aligned}.bam >> QC_report_post_alignment/NRF_PBC.txt 
```

#### Calculation of NSC and RSC

NSC and RSC are computed using `phantompeakqualtools`:

```         
Rscript run_spp.R -c=${aligned}.bam -p=30 -out=QC_report_post_alignment/NSC_RSC.txt 
```

#### Output

-   `QC_report_post_alignment/NRF_PBC.txt`

-   `QC_report_post_alignment/NSC_RSC.txt`

These metrics provide an initial assessment of ChIP-seq data quality and enrichment.

### Step 3. Post-alignment QC summarization in R

**Script:** `Post_alignment_QC.R`

#### Purpose

This step integrates alignment statistics, peak files, and metadata to compute additional ChIP-seq quality metrics using the `ChIPQC` package.

#### Input

-   Duplicate-removed BAM files (`*.unique_alignment_sorted_rd.bam`)

-   Peak files (`*.narrowPeak`)

-   A metadata table describing samples

Example metadata file (`metadata.txt`):

```         
IP    Age 
Sample1   30 
Sample2   45 
```

#### Procedure

Load required R packages:

```         
library(ChIPQC) 
library(dplyr) 
```

Prepare metadata for ChIPQC:

```         
chip_metadata <- read.table(
  'data/metadata.txt',
  sep = '\t',
  header = TRUE
)

chipqc_metadata <- data.frame(
  SampleID = chip_metadata$IP,
  Factor   = chip_metadata$Age
) %>%
  mutate(
    Tissue    = 'Spleen',  # replace with your tissue type
    Condition = 'H3K27ac', # replace with your histone mark
    bamReads  = normalizePath(
      paste0('data/Alignment/', SampleID,
             '.unique_alignment_sorted_rd.bam')
    ),
    Peaks     = normalizePath(
      paste0('data/peakcalling/', SampleID,
             '_peaks.narrowPeak')
    )
  ) %>%
  group_by(Factor) %>%
  mutate(Replicate = row_number()) %>%
  ungroup()
```

Run post-alignment quality control:

```         
chipqc.res <- ChIPQC(chipqc_metadata, annotation = 'hg19') 
```

#### Output

-   FRiP scores

-   Cross-correlation and enrichment summaries

-   A comprehensive QC object (`chipqc.res`) for downstream inspection

## Part II (Steps 4 - 6)

### Overview

This part describes the generation of peak-level features from aligned ChIP-seq data, construction of a consensus peak set, and development of histone modification–based biological age prediction models.

## Step 4. Peak calling and super-enhancer identification

**Script:** `Peakcalling.sh`

### Purpose

This step identifies ChIP-seq peaks for each sample using MACS2 and optionally detects super-enhancers using HOMER. Both narrow and broad histone marks are supported.

### Input

-   Duplicate-removed BAM files from Step 2\
    (`*.unique_alignment_sorted_rd.bam`)

-   Corresponding input/control BAM files

### Software requirements

-   MACS2

-   HOMER

### Peak calling with MACS2

Create an output directory:

```         
mkdir peakcalling 
```

#### Narrow peaks

(e.g. H3K27ac, H3K4me3)

```         
macs2 callpeak \
  -t Alignment/${aligned_IP}.rd.bam \
  -c Alignment/${aligned_Input}.rd.bam \
  -g hs \
  -f BAM \
  -n ${prefix} \
  --outdir peakcalling
```

#### Broad peaks

(e.g. H3K27me3, H3K36me3)

```         
macs2 callpeak \
  -t Alignment/${aligned_IP}.rd.bam \
  -c Alignment/${aligned_Input}.rd.bam \
  -g hs \
  -f BAM \
  --broad \
  -n ${prefix} \
  --outdir peakcalling
```

> **Note:**\
> Users are expected to run MACS2 iteratively for each IP–Input pair.\
> `${prefix}` should uniquely identify each sample.

### Output

-   `*_peaks.narrowPeak` or `*_peaks.broadPeak`

-   `*_summits.bed`

-   `*_peaks.xls`

These peak files are used as input for consensus peak construction in Step 5.

### Optional: super-enhancer identification with HOMER

Create a tag directory:

```         
mkdir tag_directory
makeTagDirectory tag_directory/${prefix} ${aligned}.rd.bam
```

Identify super-enhancers:

```         
mkdir super_enhancer
findPeaks tag_directory/${tag_dict_IP} \
  -i tag_directory/${tag_dict_Input} \
  -style super \
  -o super_enhancer/${prefix}_SE.txt
```

### Output

-   Super-enhancer annotation file (`*_SE.txt`)

## Step 5. Construction of consensus peak sets

**Script:** `Consensus_peakset.R`

### Purpose

This step generates age-specific sub-consensus peak sets using MSPC and then merges them into a final global consensus peak set across all ages.

### Input

-   MACS2 peak files (`*_peaks.narrowPeak`)

-   Sample metadata table (`metadata.txt`)

### Software requirements

-   R packages: `jsonlite`, `dplyr`, `stringr`

-   MSPC

-   BEDTools

### Generate MSPC parser configuration

A custom JSON parser is created to accommodate MACS2 p-values:

```         
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

write(
  toJSON(mspc_parser_config, pretty = TRUE, auto_unbox = TRUE),
  file = "data/peakcalling/parser_config.json"
)
```

### Identify age-specific sub-consensus peak sets

```         
chip_metadata <- read.table(
  'data/metadata.txt',
  sep = '\t',
  header = TRUE
)

peak_path <- 'data/peakcalling/'

split(chip_metadata, chip_metadata$Age) %>%
  lapply(function(x) {
    peak <- paste0(peak_path, x$IP, '_peaks.narrowPeak', collapse = ' ')
    outpath <- paste0('consensus_peak_for_', unique(x$Age))
    system(
      paste0(
        '/PATH/TO/MSPC/mspc -i ',
        peak,
        ' -r bio -c 50% -w 1e-4 -s 1e-8 ',
        '-p data/peakcalling/parser_config.json ',
        '-o data/peakcalling/', outpath
      )
    )
  })
```

### Output

-   Age-specific MSPC consensus peaks (`ConsensusPeaks.bed`)

### Generate final global consensus peak set

```         
all_peaks <- list.dirs("data/peakcalling", recursive = FALSE) %>%
  paste0('/ConsensusPeaks.bed') %>%
  lapply(function(x) {
    peak <- read.table(x, sep = '\t', header = TRUE)
    peak$name <- paste0(
      str_split(x, '/', simplify = TRUE)[,3],
      '_peak_', seq_len(nrow(peak))
    )
    peak
  }) %>%
  do.call(rbind, .) %>%
  arrange(chr, start)

write.table(
  all_peaks,
  'data/peakcalling/all_consensus_peak.bed',
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE,
  sep = '\t'
)

system(
  '/PATH/TO/BEDTOOLS/bedtools merge \
   -i data/peakcalling/all_consensus_peak.bed \
   > data/peakcalling/final_consensus_peak_set.bed'
)
```

### Output

-   `final_consensus_peak_set.bed`

This file defines the common feature space used for downstream modeling.

## Step 6. Feature extraction and biological age modeling

**Script:** `Modeling.R`

### Purpose

This step quantifies histone modification signal over the final consensus peak set, removes batch effects, identifies age-correlated peaks, and builds an elastic net regression model for biological age prediction.

### Calculate signal over consensus peaks

**Input**

-   Final consensus peak set

-   BAM files

-   Sample metadata

```         
final_consensus_peak_set <- read.table(
  'data/peakcalling/final_consensus_peak_set.bed',
  header = FALSE
)
final_consensus_peak_set$id <- paste0('peak_', seq_len(nrow(final_consensus_peak_set)))
colnames(final_consensus_peak_set)[1:3] <- c('chr','start','end')
```

Read counting and normalization using TCseq:

```         
tca <- TCA(design = design, genomicFeature = final_consensus_peak_set)
register(MulticoreParam(workers = 1))
tca <- TCseq::countReads(tca, dir = path_to_bam)
tca <- timecourseTable(tca, value = 'expression', norm.method = 'rpkm')

chip_all_signal <- TCseq::counts(
  tca,
  normalization = 'rpkm',
  log = TRUE
)
```

### Batch effect correction

```         
signal_combat_corrected <- ComBat(
  dat = chip_all_signal,
  batch = batch_metadata$Lab,
  mod = model.matrix(~ batch_metadata$Age),
  par.prior = TRUE,
  prior.plots = FALSE
)
```

PCA can be used to assess batch correction effectiveness.

### Identify age-correlated peaks

Spearman correlation is used to identify peaks associated with aging:

```         
chip_cor_peak <- cor_peak_ident(tca) 
```

### Output

-   Table of age-correlated peaks

-   Distribution of correlation coefficients

### Elastic net regression modeling

Define training and test sets:

```         
model.cva <- cva.glmnet(
  chip_all_signal.sub.train,
  as.numeric(rownames(chip_all_signal.sub.train)),
  alpha = seq(0.1, 0.9, 0.1),
  nfolds = 10
)
```

Fit the final model and predict age:

```         
model.glm <- glmnet(
  chip_all_signal.sub.train,
  as.numeric(rownames(chip_all_signal.sub.train)),
  family = 'gaussian',
  alpha = hp$alpha
)

predictions <- predict(
  model.glm,
  chip_all_signal.sub.test,
  s = hp$lambdaMin
)
```

### Output

-   Trained elastic net model

-   Predicted biological age

-   Performance metrics (RMSE, MAE, Pearson’s R)

## 🔑 Applying the clock to new datasets (summary)

To apply the trained clock to new ChIP-seq datasets:

1.  Process raw FASTQ files following Steps 1–3

2.  Quantify signal over `final_consensus_peak_set.bed`

3.  Apply the same normalization and batch correction

4.  Use the trained elastic net model to predict biological age

## 
