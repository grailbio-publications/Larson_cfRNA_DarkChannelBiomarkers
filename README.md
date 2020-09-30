## Contents

- [Overview](#overview)
- [System requirements](#system-requirements)
- [Reproducibility](#reproducibility)
- [License](./LICENSE)

# Overview

Manuscript overview:
Cell-free RNA (cfRNA) is a promising analyte for cancer detection, and a comprehensive assessment of cfRNA in individuals with and without cancer has not been conducted. We performed the first transcriptome-wide characterization of cfRNA in cancer and non-cancer participants from the Circulating Cell-free Genome Atlas (CCGA), and identified tissue- and cancer-specific genes, defined as “dark channel biomarker” (DCB) genes, that were recurrently detected in individuals with cancer.

This code consists of an R package and a set of scripts used for the analysis and figure generation for the Larson, et al.

Patient clinical data are provided in Supplementary Table 5 of the manuscript. Sequencing data, summary gene expression counts by patient, and patient metadata will be deposited in the European Genome-phenome Archive (EGA) (accession number not yet available).

The CCGA study protocol is [available here](http://clinicaltrials.gov/ct2/show/NCT02889978).

# System requirements and installation

Code was developed for use on Ubuntu 16.04.

Follow instructions at [https://cran.r-project.org/bin/linux/ubuntu/README.html](https://cran.r-project.org/bin/linux/ubuntu/README.html) to install R (version 3.6.2 was used for this analysis).

We've supplied the R package `cellfreetranscriptome` to reproduce this analysis. Version dependencies are handled with packrat.

To add the required packages, open R in the project directory:

`packrat::restore()`

This pulls the versioned dependencies from CRAN and Bioconductor. This should take about 25 minutes on a standard computer. Then, run:

`packrat::init(infer.dependencies = FALSE, options = list(vcs.ignore.lib = TRUE, vcs.ignore.src = TRUE))`

Then exit R and run `R CMD INSTALL cellfreetranscriptome_0.1.0.tar.gz`. Installation should complete in several seconds.

# Reproducibility

This analysis relies on publicly available data sets, which are described below. The code dscribed below generates the figures as shown in Larson, et al.

## Data Access

TCGA transcriptome count data is downloaded from GDC Portal with [TCGAbiolinks, installation instructions here](https://github.com/BioinformaticsFMRP/TCGAbiolinks). HR and HER2 status was determined by IHC status from clinical metadata. The sample-wise counts were compiled into TSV files, where rows correspond to genes and columns correspond to TCGA samples. We collected primary tumor counts from the HT-Seq workflow for LUSC (lung squamous), LUAD (lung adenocarcinoma), and BRCA (HR+/HER2+, HR+/HER2-, HR-/HER2+, and TNBC).

For example, in R:

`query <- TCGAbiolinks::GDCquery(project = 'TCGA-LUAD', data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "HTSeq - Counts", legacy = FALSE)`

`TCGAbiolinks::GDCdownload(query)`

`data <- TCGAbiolinks::GDCprepare(query)`

`gene_count <- as.data.frame(SummarizedExperiment::assay(data))`

`write.table(gene_count, 'tcga_luad luad_tumor_counts.tsv', quote=F, sep= '\t')`

Gene-wise counts in healthy tissue compartments were downloaded from the [GTEx consortium](https://storage.googleapis.com/gtex_analysis_v4/rna_seq_data/GTEx_Analysis_V4_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz). Sample attributes with tissue localization is downloaded from [GTEx](https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt). Median gene-wise RPMs for each computed for each tissue and saved into an .RData file as a data frame (row names are Ensembl gene ID and column names are tissue compartments).

## Code

Code used to generate the analysis for the manuscript is found in the `exec` folder. Each script has a specific set of parameters pointing to file inputs which must be downloaded (from EGA, TCGA, and GTEx). To access instructions regarding the specific data parameters required to replicate each figure, try `Rscript figure2_deconvolution.R --help`. 

For example, to generate the heatmap using dummy data:

`Rscript figure3_whole_txome_heatmap.R --count_cfrna ../inst/extdata/dummy_spliced_counts.tsv`.
