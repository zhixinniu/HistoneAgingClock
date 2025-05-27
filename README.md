# HistoneAgingClock

This repository contains code related to the manuscript *"Histone Modification Clocks for Robust Cross-Species Biological Age Prediction and Elucidating Senescence Regulation."* Currently, six scripts are included:

**Preprocessing.sh**: Performs quality control of ChIP-seq data using FastQC and MultiQC.\
**Alignment.sh**: Trims and filters low-quality reads, maps reads to the reference genome, and conducts part of the post-alignment quality control (calculating NRF, PBC, RSC, and NSC).\
**Post_alignment_QC.R**: Completes the remaining post-alignment quality control steps, including FRiP evaluation.\
**Peakcalling.sh**: Calls peaks using MACS2 and identifies super enhancers with HOMER.\
**Consensus_peakset.R**: Identifies consensus peak sets across samples.\
**Modeling.R**: Tunes hyperparameters and builds an elastic net regression model for biological age prediction.
