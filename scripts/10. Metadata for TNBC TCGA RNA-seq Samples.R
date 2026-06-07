############################################################
# STEP 3 (FINAL): Create Metadata for TNBC TCGA RNA-seq
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

library(SummarizedExperiment)

## ---- Project root ----
project_dir <- "C:/Users/abhir/Desktop/TNBC/TCGA_TNBC"
setwd(project_dir)

## ---- Load frozen TCGA object ----
brca_se <- readRDS("data_processed/brca_se.rds")

## ---- Load finalized TNBC patient list ----
tnbc_ids <- readRDS("data_processed/final_tnbc_tcga_samples.rds")
# Format: TCGA-XX-YYYY-01

cat("TNBC samples in list:", length(tnbc_ids), "\n")

## ---- Extract colData ----
clinical_df <- as.data.frame(colData(brca_se))
clinical_df$RNAseq_Barcode <- colnames(brca_se)

## ---- Patient-level matching (CORRECT) ----
clinical_df$Patient_ID <- substr(clinical_df$RNAseq_Barcode, 1, 12)
tnbc_patients <- substr(tnbc_ids, 1, 12)

metadata_tnbc <- clinical_df[
  clinical_df$Patient_ID %in% tnbc_patients &
    clinical_df$sample_type == "Primary Tumor",
]

## ---- Hard safety checks ----
stopifnot(nrow(metadata_tnbc) == length(tnbc_ids))

## ---- Create clean metadata table ----
metadata_tnbc_clean <- data.frame(
  Sample_ID    = metadata_tnbc$RNAseq_Barcode,
  Patient_ID   = metadata_tnbc$Patient_ID,
  Sample_Type  = metadata_tnbc$sample_type,
  Tumor_Status = "Tumor",
  Cohort       = "TNBC",
  stringsAsFactors = FALSE
)

## ---- Inspect ----
head(metadata_tnbc_clean)
dim(metadata_tnbc_clean)
# Expected: 116 rows × 5 columns

## ---- Save ----
write.csv(
  metadata_tnbc_clean,
  "data_processed/metadata_tnbc_tcga.csv",
  row.names = FALSE
)

saveRDS(
  metadata_tnbc_clean,
  "data_processed/metadata_tnbc_tcga.rds"
)

sessionInfo()
