# Bioconductor core
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)

# Annotation
library(biomaRt)

# Data handling
library(dplyr)
library(tibble)
library(stringr)

dir.create("data_raw", showWarnings = FALSE)
dir.create("data_processed", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("scripts", showWarnings = FALSE)