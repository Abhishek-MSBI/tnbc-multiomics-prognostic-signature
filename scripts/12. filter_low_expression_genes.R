############################################################
# STEP 5: Filter Low-Expressed Genes (TNBC TCGA RNA-seq)
############################################################

library(SummarizedExperiment)

# Load raw TNBC count matrix
counts_tnbc <- readRDS("data_processed/counts_tnbc_raw_STAR.rds")

dim(counts_tnbc)
# Expected: ~60660 x 116

# Define filtering criteria
min_count <- 10
min_samples <- ceiling(0.20 * ncol(counts_tnbc))

min_samples
# Expected: 24

# Apply filtering
keep_genes <- rowSums(counts_tnbc >= min_count) >= min_samples

# Summary of filtering
table(keep_genes)

# Filtered count matrix
counts_tnbc_filtered <- counts_tnbc[keep_genes, ]

# Dimensions after filtering
dim(counts_tnbc_filtered)

# Percentage retained
pct_retained <- round(
  (nrow(counts_tnbc_filtered) / nrow(counts_tnbc)) * 100, 2
)

cat("Genes before filtering:", nrow(counts_tnbc), "\n")
cat("Genes after filtering:", nrow(counts_tnbc_filtered), "\n")
cat("Percentage retained:", pct_retained, "%\n")

# Save filtered counts
saveRDS(
  counts_tnbc_filtered,
  file = "data_processed/counts_tnbc_filtered.rds"
)

write.csv(
  as.data.frame(counts_tnbc_filtered),
  file = "data_processed/counts_tnbc_filtered.csv"
)

# Reproducibility
sessionInfo()
