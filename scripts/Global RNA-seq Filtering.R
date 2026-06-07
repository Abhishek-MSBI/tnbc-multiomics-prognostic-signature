############################################################
# STEP 5: Global RNA-seq Low-Expression Filtering
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

library(edgeR)

## ---- Project root ----
project_dir <- "C:/Users/abhir/Desktop/TNBC/TCGA_TNBC"
setwd(project_dir)

## ---- Load raw counts ----
tnbc_counts   <- readRDS("data_processed/counts_tnbc_raw_STAR.rds")
normal_counts <- readRDS("data_processed/counts_normal_raw_STAR.rds")

## ---- Ensure common gene universe ----
common_genes <- intersect(
  rownames(tnbc_counts),
  rownames(normal_counts)
)

tnbc_counts   <- tnbc_counts[common_genes, ]
normal_counts <- normal_counts[common_genes, ]

## ---- Merge tumor + normal ----
merged_counts <- cbind(tnbc_counts, normal_counts)

group <- factor(c(
  rep("TNBC", ncol(tnbc_counts)),
  rep("Normal", ncol(normal_counts))
))

## ---- RNA-seq global filtering ----
keep <- rowSums(merged_counts >= 10) >= ceiling(0.2 * ncol(merged_counts))

filtered_counts <- merged_counts[keep, ]

cat("Genes before filtering:", nrow(merged_counts), "\n")
cat("Genes after filtering:",  nrow(filtered_counts), "\n")

## ---- Save ----
saveRDS(
  filtered_counts,
  "data_processed/counts_tnbc_normal_filtered.rds"
)

saveRDS(
  group,
  "data_processed/group_labels.rds"
)
