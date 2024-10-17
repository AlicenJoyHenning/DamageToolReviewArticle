# A review and evaluation of current damaged cell detection strategies for single cell RNA sequencing (scRNA-seq)

This repository contains the workflow, direct analysis code, and original data (where possible due to size constraints) for the results and figures generated in our article. For any concerns, questions or suggestions, please feel free to open a new  [![Issue](https://img.shields.io/badge/Issues-blue?style=flat&logo=github)](https://github.com/AlicenJoyHenning/DamageToolReviewArticle/issues)
<br>
<br>


### Repository Overview
⚪[  Analysis Overview](#-analysis-overview) | ⚪[  Analysis Scripts](#-analysis-scripts) | ⚪[  Data Availability](#-data-availability)   

<br>


## ⚪ Analysis Overview

![Alt text](https://github.com/AlicenJoyHenning/DamageToolReviewArticle/blob/main/Images/workflow-removebg-preview.png)

## ⚪ Analysis Scripts
Details for running each script are outlined here as well as within the script as comments.
<br>
<br>
| Analysis Stage | Code Link | Explanation | Input Needed | Output Expected |
|---|---|---|---|---|
| Literature Review | [R](https://github.com/AlicenJoyHenning/bioinformatics/blob/main/R/paper_cellQC.rmd) | For extracting and filtering scRNA-seq study search results from Google Scholar and plotting after manual review| NA | CSV of study titles and pie chart after manual review |
| Tool Search | [R]| For synthesizing automated searches in [scRNA-tools](https://www.scrna-tools.org/) with manual review | CSV of tool details [scRNA-tools.csv]() | CSV of potential tools, pie charts, and scatter plot for tools selected |
| Data Acquistion | [bash](SRA) | Downloading publicly available data through SRA | SRA accession number for dataset and sample of interest | FASTQ.gz files | 
| Build Alignment Indices | [bash](STAR) | Create reference genome indices for alignment | Species GTF and annotation files | STAR Index | 
| Data Processing | [bash](STARsolo) | Aligning and quanitfying FASTQ files to reference genome | FASTQ.gz files, Genome Index | matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz | 
| Data Processing | [R](Seurat) | Pre-processing the sample count matrices | matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz (ideally raw and filtered) | R object with sample scRNA-seq data (Seurat object) | 
| Tool Application | [R](Seurat) | Performing damaged detection strategies for each sample | Processed Seurat objects | CSVs of damaged labels (summary and for each dataset) | 
| Ground Truth Performance Metrics | [R](Seurat) | Calcalating the Precision, FNR, and PR-AUC for each ground truth dataset | CSV damaged labels (including truth) | Performance plots | 
| Ground Truth Comparision Metrics | [R](Seurat) | Calcalating damaged proportions, pairwise similarity scores, and unique damaged detection visualisation | CSV damaged labels (including truth) | Performance plots | 
| Non-Ground Truth Comparision Metrics | [R](Seurat) | Calcalating damaged proportions, pairwise similarity scores, and unique damaged detection visualisation | CSV damaged labels (including truth) | Performance plots | 
| Simulated Downstream Analysis | [R](Seurat) | Using a simulated data set and adding damaged cells, investigating the influence on downstream results including clustering, DGEA and enriched GO pathways | NA | Performance plots | 
| Remaining Questions | [R](Seurat) | Annotating the isolated damaged populations and integrating with control populations for visualisation, summary statistics, and DGEA for the different populations forming | R objects for ground truth samples | Performance plots | 

## ⚪ Data Availability   

<br>
<br>
#### ⚪ Data Availability



