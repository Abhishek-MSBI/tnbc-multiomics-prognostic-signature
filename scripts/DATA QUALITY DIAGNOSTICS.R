#### =====================================================
#### DATA QUALITY DIAGNOSTICS
#### =====================================================

cat("\n========== DATA QUALITY CHECKS ==========\n\n")

#### CHECK 1: Load and inspect raw files
expr <- read.csv("tcga_expression_for_lasso_cox.csv")
surv <- read.csv("tcga_survival_for_lasso_cox.csv")

cat("--- FILE DIMENSIONS ---\n")
cat("Expression file:", dim(expr), "(rows × cols)\n")
cat("Survival file:", dim(surv), "(rows × cols)\n\n")

#### CHECK 2: Column names and structure
cat("--- EXPRESSION FILE STRUCTURE ---\n")
cat("Column names:", colnames(expr)[1:5], "...\n")
cat("First column name:", colnames(expr)[1], "\n")
cat("Is 'Sample_ID' present?", "Sample_ID" %in% colnames(expr), "\n\n")

cat("--- SURVIVAL FILE STRUCTURE ---\n")
cat("Column names:", colnames(surv), "\n")
cat("Required columns present?\n")
cat("  Sample_ID:", "Sample_ID" %in% colnames(surv), "\n")
cat("  time:", "time" %in% colnames(surv), "\n")
cat("  event:", "event" %in% colnames(surv), "\n\n")

#### CHECK 3: Sample ID matching
cat("--- SAMPLE ID MATCHING ---\n")
cat("Unique samples in expression:", length(unique(expr$Sample_ID)), "\n")
cat("Unique samples in survival:", length(unique(surv$Sample_ID)), "\n")
cat("Duplicate IDs in expression:", sum(duplicated(expr$Sample_ID)), "\n")
cat("Duplicate IDs in survival:", sum(duplicated(surv$Sample_ID)), "\n")

# Check overlap
common_samples <- intersect(expr$Sample_ID, surv$Sample_ID)
cat("Samples in both files:", length(common_samples), "\n")
cat("Samples only in expression:", length(setdiff(expr$Sample_ID, surv$Sample_ID)), "\n")
cat("Samples only in survival:", length(setdiff(surv$Sample_ID, expr$Sample_ID)), "\n\n")

#### CHECK 4: Missing values
cat("--- MISSING VALUES ---\n")
expr_genes <- expr[, -1]  # Exclude Sample_ID
cat("Missing values in expression data:", sum(is.na(expr_genes)), "\n")
cat("Columns with missing values:", sum(colSums(is.na(expr_genes)) > 0), "\n")
cat("Missing values in survival time:", sum(is.na(surv$time)), "\n")
cat("Missing values in survival event:", sum(is.na(surv$event)), "\n\n")

#### CHECK 5: Data types
cat("--- DATA TYPES ---\n")
cat("Expression data class:", class(expr_genes[,1]), "\n")
cat("Are all gene columns numeric?", 
    all(sapply(expr_genes, is.numeric)), "\n")
cat("Survival time class:", class(surv$time), "\n")
cat("Survival event class:", class(surv$event), "\n\n")

#### CHECK 6: Survival data validity
cat("--- SURVIVAL DATA VALIDATION ---\n")
cat("Event values:", paste(unique(surv$event), collapse = ", "), "\n")
cat("Event counts:\n")
print(table(surv$event))
cat("Time range:", min(surv$time, na.rm = TRUE), "to", 
    max(surv$time, na.rm = TRUE), "\n")
cat("Negative time values:", sum(surv$time < 0, na.rm = TRUE), "\n")
cat("Zero time values:", sum(surv$time == 0, na.rm = TRUE), "\n\n")

#### CHECK 7: Expression data distribution
cat("--- EXPRESSION DATA DISTRIBUTION ---\n")

# After merge (your actual analysis data)
expr_clean <- expr[!duplicated(expr$Sample_ID), ]
data <- merge(expr_clean, surv, by = "Sample_ID")
X <- as.matrix(data[, !(colnames(data) %in% c("Sample_ID","time","event"))])

cat("Expression matrix dimensions:", dim(X), "\n")
cat("Expression value range:", min(X, na.rm = TRUE), "to", 
    max(X, na.rm = TRUE), "\n")
cat("Expression value mean:", mean(X, na.rm = TRUE), "\n")
cat("Expression value SD:", sd(X, na.rm = TRUE), "\n\n")

# Check for zero variance genes
zero_var_genes <- apply(X, 2, var, na.rm = TRUE) == 0
cat("Zero variance genes:", sum(zero_var_genes), "\n")
if(sum(zero_var_genes) > 0) {
  cat("  Genes with zero variance:", 
      paste(colnames(X)[zero_var_genes], collapse = ", "), "\n")
}

# Check for extreme outliers
extreme_values <- abs(scale(X)) > 10
cat("Extreme outliers (>10 SD):", sum(extreme_values, na.rm = TRUE), "\n\n")

#### CHECK 8: Gene variance distribution
gene_vars <- apply(X, 2, var, na.rm = TRUE)
cat("--- GENE VARIANCE ---\n")
cat("Min variance:", min(gene_vars, na.rm = TRUE), "\n")
cat("Median variance:", median(gene_vars, na.rm = TRUE), "\n")
cat("Max variance:", max(gene_vars, na.rm = TRUE), "\n")
cat("Low variance genes (<0.01):", sum(gene_vars < 0.01, na.rm = TRUE), "\n\n")

#### CHECK 9: Correlation structure
cat("--- CORRELATION STRUCTURE ---\n")
if(ncol(X) <= 100) {
  cor_matrix <- cor(X, use = "pairwise.complete.obs")
  high_cor <- sum(abs(cor_matrix[upper.tri(cor_matrix)]) > 0.95)
  cat("Gene pairs with |r| > 0.95:", high_cor, "\n")
} else {
  # Sample 50 genes for quick check
  sample_genes <- sample(1:ncol(X), min(50, ncol(X)))
  cor_matrix <- cor(X[, sample_genes], use = "pairwise.complete.obs")
  high_cor <- sum(abs(cor_matrix[upper.tri(cor_matrix)]) > 0.95)
  cat("Gene pairs with |r| > 0.95 (sample of 50 genes):", high_cor, "\n")
}
cat("\n")

#### CHECK 10: Sample outliers
cat("--- SAMPLE OUTLIERS ---\n")
sample_means <- rowMeans(X, na.rm = TRUE)
sample_sds <- apply(X, 1, sd, na.rm = TRUE)
cat("Sample mean range:", min(sample_means), "to", max(sample_means), "\n")
cat("Sample SD range:", min(sample_sds), "to", max(sample_sds), "\n")

outlier_samples <- abs(scale(sample_means)) > 3
cat("Outlier samples (>3 SD from mean):", sum(outlier_samples), "\n")
if(sum(outlier_samples) > 0) {
  cat("  Outlier sample IDs:", 
      paste(data$Sample_ID[outlier_samples], collapse = ", "), "\n")
}
cat("\n")

#### CHECK 11: Verify the 18 events
cat("--- FINAL MERGED DATA VALIDATION ---\n")
cat("Final sample size:", nrow(data), "\n")
cat("Final gene count:", ncol(X), "\n")
cat("Event distribution:\n")
print(table(data$event))
cat("Survival time summary:\n")
print(summary(data$time))
cat("\n")

# Check if events are concentrated in specific time periods
if(sum(data$event == 1) > 0) {
  cat("Death time distribution:\n")
  print(summary(data$time[data$event == 1]))
  cat("\nCensored time distribution:\n")
  print(summary(data$time[data$event == 0]))
}

cat("\n========== DIAGNOSTICS COMPLETE ==========\n")

#### VISUALIZATION: Create diagnostic plots
dir.create("results_TNBC/LASSO_Cox/diagnostics", 
           recursive = TRUE, showWarnings = FALSE)

# Plot 1: Expression distribution
pdf("results_TNBC/LASSO_Cox/diagnostics/expression_distribution.pdf", 
    width = 10, height = 6)
par(mfrow = c(2, 2))
hist(X, breaks = 50, main = "All Expression Values", 
     xlab = "Expression", col = "lightblue")
boxplot(X, main = "Expression by Gene", 
        xlab = "Genes (showing first 50)", las = 2)
hist(gene_vars, breaks = 30, main = "Gene Variance Distribution",
     xlab = "Variance", col = "lightgreen")
hist(sample_means, breaks = 30, main = "Sample Mean Expression",
     xlab = "Mean Expression", col = "lightyellow")
dev.off()

# Plot 2: Survival data
pdf("results_TNBC/LASSO_Cox/diagnostics/survival_distribution.pdf", 
    width = 8, height = 6)
par(mfrow = c(2, 2))
hist(data$time, breaks = 30, main = "Survival Time Distribution",
     xlab = "Time", col = "lightblue")
hist(data$time[data$event == 1], breaks = 20, 
     main = "Death Times Only", xlab = "Time", col = "red")
hist(data$time[data$event == 0], breaks = 20, 
     main = "Censored Times Only", xlab = "Time", col = "green")
plot(data$time, data$event, main = "Time vs Event",
     xlab = "Time", ylab = "Event (0=censored, 1=death)", pch = 19)
dev.off()

cat("\nDiagnostic plots saved to results_TNBC/LASSO_Cox/diagnostics/\n")

