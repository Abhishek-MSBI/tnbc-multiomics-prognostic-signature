############################################################
# DEG Summary for DESeq2 Results (TNBC vs Normal)
############################################################

library(dplyr)
library(readr)

# -------------------------------
# 1. Load DESeq2 DEG results
# -------------------------------
deg_results <- read.csv(
  "results/DGE_DESeq2/DEG_DESeq2_TNBC_vs_Normal_all_genes.csv",
  stringsAsFactors = FALSE
)

# Inspect columns (optional)
colnames(deg_results)

# -------------------------------
# 2. Define thresholds
# -------------------------------
logFC_cutoff <- 1
padj_cutoff  <- 0.05

# -------------------------------
# 3. Sanity check for required columns
# -------------------------------
required_cols <- c("log2FoldChange", "padj")
missing_cols <- setdiff(required_cols, colnames(deg_results))

if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# -------------------------------
# 4. Classify genes
# -------------------------------
deg_results <- deg_results %>%
  mutate(
    Regulation = case_when(
      log2FoldChange >=  logFC_cutoff & padj < padj_cutoff ~ "Upregulated",
      log2FoldChange <= -logFC_cutoff & padj < padj_cutoff ~ "Downregulated",
      TRUE ~ "Not_Significant"
    )
  )

# -------------------------------
# 5. Count genes
# -------------------------------
deg_counts <- table(deg_results$Regulation)

print(deg_counts)

# -------------------------------
# 6. Extract significant genes
# -------------------------------
deg_significant <- deg_results %>%
  filter(Regulation != "Not_Significant")

# -------------------------------
# 7. Save outputs
# -------------------------------

# Full DEG table with regulation
write.csv(
  deg_results,
  "results/DGE_DESeq2/DEG_DESeq2_with_regulation.csv",
  row.names = FALSE
)

# Significant DEGs only
write.csv(
  deg_significant,
  "results/DGE_DESeq2/DEG_DESeq2_significant_genes.csv",
  row.names = FALSE
)

# DEG count summary
write.csv(
  as.data.frame(deg_counts),
  "results/DGE_DESeq2/DEG_DESeq2_gene_counts_summary.csv"
)

# -------------------------------
# 8. Console summary (for logs)
# -------------------------------
cat("\nDEG SUMMARY (DESeq2)\n")
cat("----------------------------\n")
cat("Upregulated genes   :", deg_counts["Upregulated"], "\n")
cat("Downregulated genes :", deg_counts["Downregulated"], "\n")
cat("Not significant     :", deg_counts["Not_Significant"], "\n")
cat("Total genes tested  :", nrow(deg_results), "\n")
cat("Total significant   :", nrow(deg_significant), "\n")
