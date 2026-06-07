############################################################
# SUPPLEMENTARY SCRIPT S0
# Gene ID → Gene Symbol Conversion (ONE-TIME STEP)
# TCGA-BRCA / TNBC Project
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
})

## =========================================================
## Project setup
## =========================================================
project_dir <- "C:/Users/abhir/Desktop/TNBC/TCGA_TNBC"
setwd(project_dir)

dir.create("data_processed/gene_mapped",
           recursive = TRUE,
           showWarnings = FALSE)

## =========================================================
## Load raw counts (ENSEMBL IDs)
## =========================================================
counts <- readRDS("data_processed/counts_tnbc_normal_filtered.rds")

## Remove version numbers (ENSGxxxx.xx → ENSGxxxx)
rownames(counts) <- sub("\\..*", "", rownames(counts))

## =========================================================
## Connect to Ensembl (stable GRCh38)
## =========================================================
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = rownames(counts),
  mart       = mart
)

## Remove empty gene symbols
gene_map <- gene_map %>%
  filter(hgnc_symbol != "")

## =========================================================
## Merge counts with gene symbols
## =========================================================
counts_df <- as.data.frame(counts)
counts_df$ensembl_gene_id <- rownames(counts_df)

counts_merged <- counts_df %>%
  dplyr::inner_join(gene_map, by = "ensembl_gene_id") %>%
  dplyr::select(-ensembl_gene_id)

## =========================================================
## Handle duplicated gene symbols
## Keep gene with highest mean expression
## =========================================================
counts_symbol <- counts_merged %>%
  group_by(hgnc_symbol) %>%
  summarize(across(everything(), mean)) %>%
  ungroup()

counts_mat <- as.matrix(counts_symbol[, -1])
rownames(counts_mat) <- counts_symbol$hgnc_symbol

## =========================================================
## Save processed data (FINAL INPUT FOR ALL ANALYSES)
## =========================================================
saveRDS(
  counts_mat,
  "data_processed/gene_mapped/counts_TNBC_Normal_GeneSymbol.rds"
)

message("✅ Gene ID → Gene Symbol conversion completed successfully")

