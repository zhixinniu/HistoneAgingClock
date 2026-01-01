# HistoneAgingClock

![](Main_fig.png)


---

This repository contains the code used in the manuscript:

> *“Histone Modification Clocks for Robust Cross-Species Biological Age Prediction and Elucidating Senescence Regulation.”*

The pipeline implements a complete ChIP-seq–based workflow, including data preprocessing, alignment, peak calling, quality control, consensus peak construction, and statistical modeling for biological age prediction.

---

## Overview of Scripts

Currently, the repository includes the following scripts, which should be executed sequentially:

* **Preprocessing.sh**
  Performs initial quality control of raw ChIP-seq FASTQ files using **FastQC** and **MultiQC**, followed by adapter trimming and quality filtering.

* **Alignment.sh**
  Aligns trimmed reads to the reference genome using **Bowtie2**, generates sorted and duplicate-removed BAM files, and performs part of the post-alignment quality control, including calculation of **NRF**, **PBC**, **RSC**, and **NSC** metrics.

* **Post_alignment_QC.R**
  Completes post-alignment quality control using the **ChIPQC** R package, including additional QC metrics such as **FRiP**.

* **Peakcalling.sh**
  Calls ChIP-seq peaks using **MACS2** and identifies super enhancers using **HOMER**.

* **Consensus_peakset.R**
  Identifies age-specific sub-consensus peak sets using **MSPC** and merges them into a final unified consensus peak set across all samples.

* **Modeling.R**
  Tunes hyperparameters and builds an elastic net regression model to predict biological age based on histone modification ChIP-seq signals. The script also includes batch correction, feature selection, and sample-size saturation analysis.
---

## Software Requirements

### R and R Packages

The following R packages are required to run the analysis scripts:

* **dplyr**
* **jsonlite**
* **stringr**
* **reshape2**
* **data.table**
* **ggplot2**
* **BiocParallel**
* **TCseq**
* **ChIPQC**
* **glmnet**
* **glmnetUtils**
* **Metrics**
* **sva**
* **factoextra**
* **FactoMineR**

Most Bioconductor packages can be installed via `BiocManager`, for example:

```r
BiocManager::install(c("ChIPQC", "TCseq", "BiocParallel"))
```
Other packages can be installed via CRAN, for example:

```
install.packages("dplyr")
```
---

### Command-Line Tools (Linux)

The following command-line tools must be installed and available in the system `PATH`:

* **bedtools**
* **samtools**
* **trim_galore**
* **fastqc**
* **multiqc**
* **bowtie2**
* **macs2**
* **homer**
* **multiqc**
  
These tools can be installed via `conda`, for example:
```
conda install -y bedtools
```
For installation of `MSPC`, please follow the user manual on https://genometric.github.io/MSPC/

---

## Notes

* All scripts assume standard ChIP-seq input formats (FASTQ, BAM, narrowPeak/broadPeak).
* Paths to external tools (e.g., MSPC, bedtools, HOMER) may need to be adjusted depending on the local installation.
* The pipeline was developed and tested on Linux-based systems.

---

