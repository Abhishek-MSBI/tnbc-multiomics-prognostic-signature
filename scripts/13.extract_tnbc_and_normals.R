## =========================================
## Extract TNBC tumors and solid tissue normals
## =========================================

rm(list = ls())
options(stringsAsFactors = FALSE)

library(SummarizedExperiment)

## ---- Project root ----
project_dir <- "C:/Users/abhir/Desktop/TNBC/TCGA_TNBC"
setwd(project_dir)

## ---- Load frozen TCGA object ----
brca_se <- readRDS("data_processed/brca_se.rds")

## ---- Load finalized TNBC sample list (CANONICAL) ----
tnbc_ids <- readRDS(
  file.path("data_processed", "final_tnbc_tcga_samples.rds")
)

## ---- Sanity check ----
cat("TNBC samples in list:", length(tnbc_ids), "\n")

## ---- Identify sample types ----
sample_types <- colData(brca_se)$sample_type

tumor_se  <- brca_se[, sample_types == "Primary Tumor"]
normal_se <- brca_se[, sample_types == "Solid Tissue Normal"]

cat("Primary tumors:", ncol(tumor_se), "\n")
cat("Solid tissue normals:", ncol(normal_se), "\n")

## ---- Match TNBC tumors ----
tumor_barcodes_15 <- substr(colnames(tumor_se), 1, 15)

tnbc_se <- tumor_se[, tumor_barcodes_15 %in% tnbc_ids]

cat("TNBC tumors retained:", ncol(tnbc_se), "\n")

stopifnot(ncol(tnbc_se) > 0)
stopifnot(ncol(normal_se) > 0)

## ---- Extract raw STAR counts ----
tnbc_counts   <- assay(tnbc_se)
normal_counts <- assay(normal_se)

## ---- Save ----
dir.create("data_processed", showWarnings = FALSE)

saveRDS(tnbc_counts,   "data_processed/counts_tnbc_raw_STAR.rds")
saveRDS(normal_counts, "data_processed/counts_normal_raw_STAR.rds")

cat("✅ TNBC and normal STAR count matrices saved\n")
