############################################################
# STEP 4B: Extract and save NORMAL raw STAR counts
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

library(SummarizedExperiment)

## ---- Project root ----
project_dir <- "C:/Users/abhir/Desktop/TNBC/TCGA_TNBC"
setwd(project_dir)

## ---- Load frozen TCGA object ----
brca_se <- readRDS("data_processed/brca_se.rds")

## ---- Identify normal samples ----
sample_types <- colData(brca_se)$sample_type

normal_se <- brca_se[, sample_types == "Solid Tissue Normal"]

cat("Normal samples identified:", ncol(normal_se), "\n")

stopifnot(ncol(normal_se) > 0)

## ---- Extract raw STAR counts ----
normal_counts <- assay(normal_se)

## ---- Save ----
saveRDS(
  normal_counts,
  "data_processed/counts_normal_raw_STAR.rds"
)

cat("✅ Normal raw STAR count matrix saved\n")
