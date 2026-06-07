########################################################
# STEP 2 (FINAL): Map cBioPortal TNBC to TCGA RNA-seq
########################################################

library(SummarizedExperiment)

# Load TNBC sample IDs (cBioPortal)
tnbc_cbio <- readRDS("data_processed/tnbc_sample_ids.rds")

length(tnbc_cbio)
# Expected: 117

# Extract TCGA RNA-seq sample barcodes
tcga_rnaseq_samples <- colnames(brca_se)

length(tcga_rnaseq_samples)
head(tcga_rnaseq_samples)

# IMPORTANT: Match using first 15 characters
tnbc_cbio_short  <- substr(tnbc_cbio, 1, 15)
tcga_rnaseq_short <- substr(tcga_rnaseq_samples, 1, 15)

# Find overlap
overlap_idx <- tcga_rnaseq_short %in% tnbc_cbio_short

tnbc_tcga_overlap <- tcga_rnaseq_samples[overlap_idx]

length(tnbc_tcga_overlap)
head(tnbc_tcga_overlap)

# Identify missing samples
missing_in_tcga <- setdiff(tnbc_cbio_short, tcga_rnaseq_short)

# Summary (for Methods)
cat("TNBC samples from cBioPortal:", length(tnbc_cbio), "\n")
cat("TNBC samples with TCGA RNA-seq:", length(tnbc_tcga_overlap), "\n")
cat("TNBC samples missing RNA-seq:", length(missing_in_tcga), "\n")

# Save final mapped TNBC RNA-seq samples
saveRDS(
  tnbc_tcga_overlap,
  file = "data_processed/final_tnbc_tcga_samples.rds"
)

write.csv(
  data.frame(Sample_ID = tnbc_tcga_overlap),
  file = "data_processed/final_tnbc_tcga_samples.csv",
  row.names = FALSE
)

# Save missing list (important for transparency)
write.csv(
  data.frame(Missing_Sample_ID = missing_in_tcga),
  file = "data_processed/tnbc_samples_missing_rnaseq.csv",
  row.names = FALSE
)

# Reproducibility
sessionInfo()
