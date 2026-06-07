## =====================================================
## TNBC Tumor-Only VST Normalization Script
## =====================================================

library(DESeq2)

## -----------------------------
## 1️⃣ Load Raw Counts
## -----------------------------

counts <- read.csv(
  "C:/Users/abhir/Desktop/TNBC/TCGA_TNBC/data_processed/counts_tnbc_filtered.csv",
  row.names = 1,
  check.names = FALSE
)

cat("Raw count matrix loaded.\n")
cat("Dimensions (genes x samples):", dim(counts), "\n")


## -----------------------------
## 2️⃣ Filter Very Low Count Genes (Recommended)
## -----------------------------

keep <- rowSums(counts >= 10) >= (0.1 * ncol(counts))
counts_filtered <- counts[keep, ]

cat("After filtering low counts:\n")
cat("Dimensions:", dim(counts_filtered), "\n")


## -----------------------------
## 3️⃣ Create Dummy Metadata
## -----------------------------

sample_info <- data.frame(
  row.names = colnames(counts_filtered),
  condition = rep("TNBC", ncol(counts_filtered))
)


## -----------------------------
## 4️⃣ Create DESeq2 Object
## -----------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_info,
  design = ~ 1
)


## -----------------------------
## 5️⃣ Perform VST Normalization
## -----------------------------

vsd <- vst(dds, blind = TRUE)

expr_vst <- assay(vsd)

cat("VST normalization completed.\n")
cat("Final VST dimensions:", dim(expr_vst), "\n")


## -----------------------------
## 6️⃣ Save Tumor-Only VST Matrix
## -----------------------------

write.csv(
  expr_vst,
  file = "tcga_tumor_vst_expression_matrix.csv"
)

cat("Tumor-only VST matrix saved successfully.\n")

