############################################################
# TNBC External Validation Script - Combined Cohorts
# Author: Generated Script
# Date: 2026-02-16
# Description: External validation of LASSO Cox model on 
#              combined cohort1 (n=42) and cohort2 (n=15)
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set random seed for reproducibility
set.seed(123)

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(biomaRt)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(pROC)
  library(ComplexHeatmap)
  library(circlize)
})

############################################################
# CREATE OUTPUT DIRECTORIES
############################################################

dir.create("external_validation_results", showWarnings = FALSE)
dir.create("external_validation_results/plots", showWarnings = FALSE)
dir.create("external_validation_results/tables", showWarnings = FALSE)
dir.create("external_validation_results/rds_objects", showWarnings = FALSE)

############################################################
# 1️⃣ LOAD AND COMBINE RAW COUNT MATRICES
############################################################

cat("\n========================================\n")
cat("STEP 1: Loading Raw Count Matrices\n")
cat("========================================\n")

# Load cohort1
counts_cohort1 <- read.csv(
  "cohort1.csv",
  row.names = 1,
  check.names = FALSE
)
cat("Cohort1 loaded: ", nrow(counts_cohort1), "genes x", ncol(counts_cohort1), "samples\n")

# Load cohort2
counts_cohort2 <- read.csv(
  "cohort2.csv",
  row.names = 1,
  check.names = FALSE
)
cat("Cohort2 loaded: ", nrow(counts_cohort2), "genes x", ncol(counts_cohort2), "samples\n")

# Remove Ensembl version numbers if present
rownames(counts_cohort1) <- sub("\\..*", "", rownames(counts_cohort1))
rownames(counts_cohort2) <- sub("\\..*", "", rownames(counts_cohort2))

# Get common genes between cohorts
common_genes <- intersect(rownames(counts_cohort1), rownames(counts_cohort2))
cat("Common genes between cohorts:", length(common_genes), "\n")

# Subset to common genes
counts_cohort1 <- counts_cohort1[common_genes, ]
counts_cohort2 <- counts_cohort2[common_genes, ]

# Combine count matrices
counts_combined <- cbind(counts_cohort1, counts_cohort2)

# Create cohort labels
cohort_labels <- data.frame(
  Sample_ID = colnames(counts_combined),
  Cohort = c(
    rep("Cohort1", ncol(counts_cohort1)),
    rep("Cohort2", ncol(counts_cohort2))
  ),
  row.names = colnames(counts_combined)
)

cat("\nCombined matrix dimensions:", nrow(counts_combined), "genes x", ncol(counts_combined), "samples\n")
cat("Cohort1 samples:", sum(cohort_labels$Cohort == "Cohort1"), "\n")
cat("Cohort2 samples:", sum(cohort_labels$Cohort == "Cohort2"), "\n")

# Save combined raw counts
saveRDS(
  list(
    counts_combined = counts_combined,
    cohort_labels = cohort_labels,
    counts_cohort1 = counts_cohort1,
    counts_cohort2 = counts_cohort2
  ),
  "external_validation_results/rds_objects/01_raw_counts_combined.rds"
)

############################################################
# 2️⃣ FILTER LOW-COUNT GENES
############################################################

cat("\n========================================\n")
cat("STEP 2: Filtering Low-Count Genes\n")
cat("========================================\n")

# Apply same filtering rule as training
keep <- rowSums(counts_combined >= 10) >= (0.1 * ncol(counts_combined))
counts_filtered <- counts_combined[keep, ]

cat("Genes before filtering:", nrow(counts_combined), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")
cat("Genes removed:", nrow(counts_combined) - nrow(counts_filtered), "\n")

# Save filtered counts
saveRDS(
  counts_filtered,
  "external_validation_results/rds_objects/02_counts_filtered.rds"
)

############################################################
# 3️⃣ VST NORMALIZATION
############################################################

cat("\n========================================\n")
cat("STEP 3: VST Normalization\n")
cat("========================================\n")

# Create DESeq2 dataset
sample_info <- data.frame(
  row.names = colnames(counts_filtered),
  cohort = cohort_labels[colnames(counts_filtered), "Cohort"]
)

# Convert to numeric matrix first
counts_filtered <- as.matrix(counts_filtered)

# Ensure numeric mode
mode(counts_filtered) <- "numeric"

# Round (safe for raw counts)
counts_filtered <- round(counts_filtered)

# Convert to integer
mode(counts_filtered) <- "integer"


dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_info,
  design = ~ 1
)

# Apply VST
vsd <- vst(dds, blind = TRUE)
expr_vst <- assay(vsd)

cat("VST completed.\n")
cat("Expression range: [", round(min(expr_vst), 2), ",", round(max(expr_vst), 2), "]\n")
cat("Expression mean:", round(mean(expr_vst), 2), "\n")
cat("Expression SD:", round(sd(expr_vst), 2), "\n")

# Save VST data
saveRDS(
  list(
    vsd = vsd,
    expr_vst = expr_vst
  ),
  "external_validation_results/rds_objects/03_vst_normalized.rds"
)

############################################################
# 4️⃣ MAP ENSEMBL TO HGNC SYMBOLS
############################################################

cat("\n========================================\n")
cat("STEP 4: Gene Symbol Mapping\n")
cat("========================================\n")

# Connect to BioMart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get mapping
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(expr_vst),
  mart = mart
)

cat("Ensembl IDs queried:", nrow(expr_vst), "\n")
cat("Mappings retrieved:", nrow(mapping), "\n")

# Merge with expression data
expr_df <- data.frame(
  ensembl_gene_id = rownames(expr_vst),
  expr_vst,
  check.names = FALSE
)

expr_mapped <- merge(
  expr_df,
  mapping,
  by = "ensembl_gene_id"
)

# Remove empty symbols
expr_mapped <- expr_mapped[expr_mapped$hgnc_symbol != "", ]

cat("Genes with valid HGNC symbols:", nrow(expr_mapped), "\n")

############################################################
# 5️⃣ COLLAPSE DUPLICATE GENE SYMBOLS
############################################################

cat("\n========================================\n")
cat("STEP 5: Collapsing Duplicates\n")
cat("========================================\n")

expr_collapsed <- expr_mapped %>%
  group_by(hgnc_symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

cat("Unique genes after collapsing:", nrow(expr_collapsed), "\n")

# Save gene mapping
saveRDS(
  list(
    mapping = mapping,
    expr_mapped = expr_mapped,
    expr_collapsed = expr_collapsed
  ),
  "external_validation_results/rds_objects/04_gene_mapping.rds"
)

############################################################
# 6️⃣ LOAD TRAINED LASSO MODEL
############################################################

cat("\n========================================\n")
cat("STEP 6: Loading Trained Model\n")
cat("========================================\n")

model_objects <- readRDS(
  "results_TNBC/LASSO_Cox/08_model_objects/complete_model_objects.rds"
)

final_genes <- model_objects$final_genes
lasso_model <- model_objects$final_model
median_risk <- median(model_objects$risk_scores)

cat("Signature genes in model:", length(final_genes), "\n")
cat("Median risk score (training):", round(median_risk, 4), "\n")
cat("\nSignature genes:\n")
print(final_genes)

############################################################
# 7️⃣ EXTRACT SIGNATURE GENES
############################################################

cat("\n========================================\n")
cat("STEP 7: Extracting Signature Genes\n")
cat("========================================\n")

expr_sig <- expr_collapsed[
  expr_collapsed$hgnc_symbol %in% final_genes,
]

missing_genes <- setdiff(final_genes, expr_sig$hgnc_symbol)
if(length(missing_genes) > 0){
  cat("WARNING: Missing genes in validation data:\n")
  print(missing_genes)
  stop("ERROR: Not all signature genes found in validation dataset.")
}

cat("All", length(final_genes), "signature genes found!\n")

# Convert to matrix format
expr_sig <- as.data.frame(expr_sig)
rownames(expr_sig) <- expr_sig$hgnc_symbol
expr_sig$hgnc_symbol <- NULL

############################################################
# 8️⃣ PREPARE PREDICTION MATRIX
############################################################

cat("\n========================================\n")
cat("STEP 8: Preparing Prediction Matrix\n")
cat("========================================\n")

# Transpose to samples x genes
X_test <- t(as.matrix(expr_sig))

# Ensure column order matches training
X_test <- X_test[, final_genes]

cat("Prediction matrix dimensions:", nrow(X_test), "samples x", ncol(X_test), "genes\n")

############################################################
# 9️⃣ PREDICT RISK SCORES
############################################################

cat("\n========================================\n")
cat("STEP 9: Predicting Risk Scores\n")
cat("========================================\n")

# Calculate risk scores
risk_scores <- as.numeric(
  predict(lasso_model,
          newx = X_test,
          type = "link")
)

# Create results data frame
risk_results <- data.frame(
  Sample_ID = rownames(X_test),
  Cohort = cohort_labels[rownames(X_test), "Cohort"],
  Risk_Score = risk_scores,
  Risk_Group = ifelse(risk_scores > median_risk, "High", "Low"),
  row.names = rownames(X_test)
)

cat("\nRisk Score Summary:\n")
print(summary(risk_scores))

cat("\nRisk Group Distribution:\n")
print(table(risk_results$Risk_Group))

cat("\nRisk Groups by Cohort:\n")
print(table(risk_results$Cohort, risk_results$Risk_Group))

# Save risk scores
write.csv(
  risk_results,
  "external_validation_results/tables/combined_risk_scores.csv",
  row.names = FALSE
)

saveRDS(
  risk_results,
  "external_validation_results/rds_objects/05_risk_scores.rds"
)

############################################################
# 🔟 GENERATE VISUALIZATIONS
############################################################

cat("\n========================================\n")
cat("STEP 10: Generating Visualizations\n")
cat("========================================\n")

## Plot 1: Risk Score Distribution by Cohort
cat("Creating Plot 1: Risk Score Distribution...\n")

p1 <- ggplot(risk_results, aes(x = Risk_Score, fill = Cohort)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = median_risk, linetype = "dashed", 
             color = "red", size = 1) +
  annotate("text", x = median_risk, y = Inf, 
           label = paste0("Median (training) = ", round(median_risk, 3)),
           hjust = -0.1, vjust = 2, color = "red") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic(base_size = 14) +
  labs(
    title = "Risk Score Distribution - External Validation",
    subtitle = paste0("Combined Cohorts (n=", nrow(risk_results), ")"),
    x = "Risk Score",
    y = "Frequency",
    fill = "Cohort"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

tiff("external_validation_results/plots/01_risk_score_distribution.tiff",
     width = 10, height = 6, units = "in", res = 600, compression = "lzw")
print(p1)
dev.off()

## Plot 2: Risk Score Boxplot by Cohort and Risk Group
cat("Creating Plot 2: Risk Score Boxplots...\n")

p2 <- ggplot(risk_results, aes(x = Risk_Group, y = Risk_Score, fill = Risk_Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  facet_wrap(~ Cohort) +
  geom_hline(yintercept = median_risk, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("High" = "#E74C3C", "Low" = "#3498DB")) +
  theme_classic(base_size = 14) +
  labs(
    title = "Risk Score by Group and Cohort",
    x = "Risk Group",
    y = "Risk Score",
    fill = "Risk Group"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

tiff("external_validation_results/plots/02_risk_score_boxplot.tiff",
     width = 10, height = 6, units = "in", res = 600, compression = "lzw")
print(p2)
dev.off()

## Plot 3: Risk Score Waterfall Plot
cat("Creating Plot 3: Risk Score Waterfall...\n")

risk_ordered <- risk_results %>%
  arrange(Risk_Score) %>%
  mutate(Sample_Order = 1:n())

p3 <- ggplot(risk_ordered, aes(x = Sample_Order, y = Risk_Score, fill = Risk_Group)) +
  geom_bar(stat = "identity", width = 1) +
  geom_hline(yintercept = median_risk, linetype = "dashed", 
             color = "black", size = 1) +
  scale_fill_manual(values = c("High" = "#E74C3C", "Low" = "#3498DB")) +
  theme_classic(base_size = 14) +
  labs(
    title = "Risk Score Waterfall Plot - All Samples",
    subtitle = paste0("n=", nrow(risk_results), " (Cohort1: ", 
                      sum(risk_results$Cohort == "Cohort1"), 
                      ", Cohort2: ", sum(risk_results$Cohort == "Cohort2"), ")"),
    x = "Sample (ordered by risk score)",
    y = "Risk Score",
    fill = "Risk Group"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

tiff("external_validation_results/plots/03_risk_score_waterfall.tiff",
     width = 12, height = 6, units = "in", res = 600, compression = "lzw")
print(p3)
dev.off()

## Plot 4: Signature Gene Expression Heatmap
cat("Creating Plot 4: Signature Gene Heatmap...\n")

# Prepare heatmap data
heatmap_data <- expr_sig[final_genes, ]
heatmap_data <- as.matrix(heatmap_data)

# Order samples by risk score
sample_order <- risk_results %>%
  arrange(Risk_Score) %>%
  pull(Sample_ID)

heatmap_data <- heatmap_data[, sample_order]

# Column annotations
col_anno <- data.frame(
  Cohort = risk_results[sample_order, "Cohort"],
  Risk_Group = risk_results[sample_order, "Risk_Group"],
  row.names = sample_order
)

# Define colors
anno_colors <- list(
  Cohort = c("Cohort1" = "#66C2A5", "Cohort2" = "#FC8D62"),
  Risk_Group = c("High" = "#E74C3C", "Low" = "#3498DB")
)

# Create heatmap
tiff("external_validation_results/plots/04_signature_genes_heatmap.tiff",
     width = 14, height = 10, units = "in", res = 600, compression = "lzw")

pheatmap(
  heatmap_data,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = col_anno,
  annotation_colors = anno_colors,
  show_colnames = FALSE,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  border_color = NA,
  fontsize = 10,
  fontsize_row = 8,
  main = paste0("Signature Gene Expression Heatmap (n=", length(final_genes), " genes)"),
  angle_col = 90
)

dev.off()

## Plot 5: Individual Cohort Distributions
cat("Creating Plot 5: Individual Cohort Plots...\n")

p5a <- ggplot(risk_results %>% filter(Cohort == "Cohort1"), 
              aes(x = Risk_Score, fill = Risk_Group)) +
  geom_histogram(bins = 20, alpha = 0.7) +
  geom_vline(xintercept = median_risk, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("High" = "#E74C3C", "Low" = "#3498DB")) +
  theme_classic(base_size = 12) +
  labs(
    title = "Cohort 1 Risk Distribution",
    x = "Risk Score",
    y = "Frequency"
  )

p5b <- ggplot(risk_results %>% filter(Cohort == "Cohort2"), 
              aes(x = Risk_Score, fill = Risk_Group)) +
  geom_histogram(bins = 15, alpha = 0.7) +
  geom_vline(xintercept = median_risk, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("High" = "#E74C3C", "Low" = "#3498DB")) +
  theme_classic(base_size = 12) +
  labs(
    title = "Cohort 2 Risk Distribution",
    x = "Risk Score",
    y = "Frequency"
  )

tiff("external_validation_results/plots/05_individual_cohort_distributions.tiff",
     width = 12, height = 5, units = "in", res = 600, compression = "lzw")
grid.arrange(p5a, p5b, ncol = 2)
dev.off()

## Plot 6: Risk Score Comparison Between Cohorts
cat("Creating Plot 6: Cohort Comparison...\n")

p6 <- ggplot(risk_results, aes(x = Cohort, y = Risk_Score, fill = Cohort)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = 16) +
  geom_hline(yintercept = median_risk, linetype = "dashed", color = "red") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 14) +
  labs(
    title = "Risk Score Comparison Between Cohorts",
    subtitle = paste0("Median (training) = ", round(median_risk, 3)),
    x = "Cohort",
    y = "Risk Score"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  stat_compare_means(method = "wilcox.test", label.y = max(risk_results$Risk_Score) * 1.1)

tiff("external_validation_results/plots/06_cohort_comparison.tiff",
     width = 8, height = 7, units = "in", res = 600, compression = "lzw")
print(p6)
dev.off()

############################################################
# SUMMARY STATISTICS
############################################################

cat("\n========================================\n")
cat("GENERATING SUMMARY STATISTICS\n")
cat("========================================\n")

summary_stats <- risk_results %>%
  group_by(Cohort, Risk_Group) %>%
  summarise(
    N = n(),
    Mean_Risk = mean(Risk_Score),
    SD_Risk = sd(Risk_Score),
    Median_Risk = median(Risk_Score),
    Min_Risk = min(Risk_Score),
    Max_Risk = max(Risk_Score),
    .groups = "drop"
  )

write.csv(
  summary_stats,
  "external_validation_results/tables/summary_statistics.csv",
  row.names = FALSE
)

cat("\nSummary Statistics:\n")
print(summary_stats)

# Overall statistics
overall_stats <- data.frame(
  Total_Samples = nrow(risk_results),
  Cohort1_N = sum(risk_results$Cohort == "Cohort1"),
  Cohort2_N = sum(risk_results$Cohort == "Cohort2"),
  High_Risk_N = sum(risk_results$Risk_Group == "High"),
  Low_Risk_N = sum(risk_results$Risk_Group == "Low"),
  Mean_Risk_Score = mean(risk_results$Risk_Score),
  SD_Risk_Score = sd(risk_results$Risk_Score),
  Median_Risk_Score = median(risk_results$Risk_Score),
  Training_Median = median_risk
)

write.csv(
  overall_stats,
  "external_validation_results/tables/overall_statistics.csv",
  row.names = FALSE
)

############################################################
# SAVE SESSION INFO
############################################################

session_info <- sessionInfo()
saveRDS(
  session_info,
  "external_validation_results/rds_objects/06_session_info.rds"
)

cat("\n========================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("========================================\n")
cat("\nResults saved in: external_validation_results/\n")
cat("  - plots/: High-resolution TIFF images (600 DPI)\n")
cat("  - tables/: CSV summary tables\n")
cat("  - rds_objects/: R objects for reproducibility\n")
cat("\nTotal samples analyzed:", nrow(risk_results), "\n")
cat("  Cohort1:", sum(risk_results$Cohort == "Cohort1"), "\n")
cat("  Cohort2:", sum(risk_results$Cohort == "Cohort2"), "\n")
cat("\nHigh risk:", sum(risk_results$Risk_Group == "High"), 
    "(", round(100*mean(risk_results$Risk_Group == "High"), 1), "%)\n")
cat("Low risk:", sum(risk_results$Risk_Group == "Low"), 
    "(", round(100*mean(risk_results$Risk_Group == "Low"), 1), "%)\n")

############################################################
# END OF SCRIPT
############################################################