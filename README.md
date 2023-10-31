# RNAseq
This repository contain the scripts for the differential and enrichment analysis

# RNA-seq Data Analysis Script

## Overview

This script is designed for the analysis of RNA-seq data using the R programming language and several R packages, including DESeq2, data.table, ggplot2, and others. It performs various tasks such as quality control, differential gene expression analysis, data visualization, and result export.

## Prerequisites

Before running the script, make sure you have the following software and packages installed:

- R (version 4.2.2 or later)
- R packages: DESeq2, data.table, ggplot2, viridis, hrbrthemes, pheatmap, tidyverse, umap

## Usage

To run the script, you should provide two command line arguments:

1. `GSE_id`: The GEO identifier for the RNA-seq dataset you want to analyze.
2. `sample_info.tsv`: A tab-separated file containing sample information, including sample IDs and group labels.

For example, you can run the script using the following command:

```shell
Rscript script.R GSE123456 sample_info.tsv
