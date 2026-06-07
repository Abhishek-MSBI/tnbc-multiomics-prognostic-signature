############################################################
# STEP 4: Extract Raw STAR Count Matrix for TNBC TCGA RNA-seq
############################################################

library(SummarizedExperiment)

# Load final TNBC TCGA RNA-seq sample IDs
tnbc_tcga_samples <- readRDS("data_processed/final_tnbc_tcga_samples.rds")

length(tnbc_tcga_samples)
# Expected: 116

# Inspect available assays (for record)
assayNames(brca_se)

# Extract raw counts from STAR pipeline (unstranded)
raw_counts <- assay(brca_se, "unstranded")

# Dimensions before subsetting
dim(raw_counts)

# Subset to TNBC samples
counts_tnbc <- raw_counts[, colnames(raw_counts) %in% tnbc_tcga_samples]

# Reorder columns to EXACTLY match TNBC sample order
counts_tnbc <- counts_tnbc[, tnbc_tcga_samples]

# HARD SAFETY CHECK
stopifnot(all(colnames(counts_tnbc) == tnbc_tcga_samples))

# Inspect counts
dim(counts_tnbc)
summary(as.numeric(counts_tnbc[,1]))

# Save raw count matrix
saveRDS(
  counts_tnbc,
  file = "data_processed/counts_tnbc_raw_STAR.rds"
)

write.csv(
  as.data.frame(counts_tnbc),
  file = "data_processed/counts_tnbc_raw_STAR.csv"
)

# Reproducibility
sessionInfo()
