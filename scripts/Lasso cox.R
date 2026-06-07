#### =====================================================
#### PUBLICATION-READY LASSO-COX REGRESSION ANALYSIS
#### Triple Negative Breast Cancer (TNBC) Prognostic Signature
#### 
#### Comprehensive workflow including:
#### - Data quality diagnostics
#### - Univariate Cox pre-screening
#### - Consensus LASSO selection (100 iterations)
#### - Pre-filtered LASSO approach
#### - Bootstrap validation
#### - Risk stratification & survival analysis
#### 
#### All outputs: 600 dpi TIFF format
#### =====================================================

#### Load required packages
required_packages <- c("tidyverse", "glmnet", "survival", "survminer", 
                       "pheatmap", "RColorBrewer", "gridExtra")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

#### Set global parameters
set.seed(123)
N_CONSENSUS_ITERATIONS <- 100
N_BOOTSTRAP_VALIDATION <- 200
CONSENSUS_THRESHOLD <- 0.50  # 50% selection frequency
PREFILTER_TOP_N <- 20  # Top N genes for pre-filtering approach

#### Create output directories
dir.create("results_TNBC/LASSO_Cox/01_diagnostics", recursive = TRUE, showWarnings = FALSE)
dir.create("results_TNBC/LASSO_Cox/02_univariate", recursive = TRUE, showWarnings = FALSE)
dir.create("results_TNBC/LASSO_Cox/03_consensus", recursive = TRUE, showWarnings = FALSE)
dir.create("results_TNBC/LASSO_Cox/04_prefiltered", recursive = TRUE, showWarnings = FALSE)
dir.create("results_TNBC/LASSO_Cox/05_final_model", recursive = TRUE, showWarnings = FALSE)
dir.create("results_TNBC/LASSO_Cox/06_validation", recursive = TRUE, showWarnings = FALSE)

#### Start logging
sink("results_TNBC/LASSO_Cox/analysis_log.txt", split = TRUE)
cat("\n==========================================================\n")
cat("LASSO-COX ANALYSIS FOR TNBC PROGNOSTIC SIGNATURE\n")
cat("Analysis started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==========================================================\n\n")

#### =====================================================
#### STEP 1: DATA LOADING AND QUALITY CHECKS
#### =====================================================

cat("\n========== STEP 1: DATA LOADING ==========\n")

expr <- read.csv("tcga_expression_for_lasso_cox.csv", stringsAsFactors = FALSE)
surv <- read.csv("tcga_survival_for_lasso_cox.csv", stringsAsFactors = FALSE)

cat("Expression file dimensions:", dim(expr), "\n")
cat("Survival file dimensions:", dim(surv), "\n")

#### Remove duplicates
expr <- expr[!duplicated(expr$Sample_ID), ]
cat("After removing duplicates:", nrow(expr), "samples\n")

#### Merge datasets
data <- merge(expr, surv, by = "Sample_ID")
cat("Merged dataset:", nrow(data), "samples\n")

#### Create expression matrix
X <- as.matrix(data[, !(colnames(data) %in% c("Sample_ID", "time", "event"))])
y <- Surv(time = data$time, event = data$event)

cat("\n--- Data Summary ---\n")
cat("Samples:", nrow(X), "\n")
cat("Genes:", ncol(X), "\n")
cat("Events:", sum(data$event), "\n")
cat("Censored:", sum(data$event == 0), "\n")
cat("Event rate:", round(100 * sum(data$event) / nrow(data), 1), "%\n")
cat("Events per variable (EPV):", round(sum(data$event) / ncol(X), 3), "\n")

#### Data quality diagnostics
cat("\n--- Data Quality Checks ---\n")
cat("Missing values in X:", sum(is.na(X)), "\n")
cat("Missing values in survival:", sum(is.na(data$time)) + sum(is.na(data$event)), "\n")
cat("Expression range:", round(min(X), 2), "to", round(max(X), 2), "\n")
cat("Survival time range:", round(min(data$time), 1), "to", round(max(data$time), 1), "\n")

zero_var_genes <- apply(X, 2, var) == 0
cat("Zero variance genes:", sum(zero_var_genes), "\n")

if(sum(zero_var_genes) > 0) {
  cat("WARNING: Removing", sum(zero_var_genes), "zero-variance genes\n")
  X <- X[, !zero_var_genes]
}

#### Save diagnostic information
diagnostic_info <- data.frame(
  Metric = c("Total_Samples", "Total_Genes", "Events", "Censored", 
             "Event_Rate_Percent", "EPV", "Min_Expression", "Max_Expression",
             "Min_Time", "Max_Time", "Zero_Var_Genes"),
  Value = c(nrow(X), ncol(X), sum(data$event), sum(data$event == 0),
            round(100 * sum(data$event) / nrow(data), 2),
            round(sum(data$event) / ncol(X), 3),
            round(min(X), 2), round(max(X), 2),
            round(min(data$time), 1), round(max(data$time), 1),
            sum(zero_var_genes))
)

write.csv(diagnostic_info,
          "results_TNBC/LASSO_Cox/01_diagnostics/data_summary.csv",
          row.names = FALSE)

#### Diagnostic plots
tiff("results_TNBC/LASSO_Cox/01_diagnostics/data_quality_plots.tiff",
     width = 12, height = 10, units = "in", res = 600, compression = "lzw")

par(mfrow = c(3, 3), mar = c(4, 4, 3, 2))

# Expression distribution
hist(X, breaks = 50, main = "Expression Value Distribution",
     xlab = "Expression", col = "lightblue", border = "white")

# Gene variance
gene_vars <- apply(X, 2, var)
hist(gene_vars, breaks = 30, main = "Gene Variance Distribution",
     xlab = "Variance", col = "lightgreen", border = "white")

# Sample means
sample_means <- rowMeans(X)
hist(sample_means, breaks = 30, main = "Sample Mean Expression",
     xlab = "Mean Expression", col = "lightyellow", border = "white")

# Survival time
hist(data$time, breaks = 30, main = "Survival Time Distribution",
     xlab = "Time (days)", col = "lightcoral", border = "white")

# Death times only
hist(data$time[data$event == 1], breaks = 20, 
     main = paste("Death Times (n =", sum(data$event), ")"),
     xlab = "Time (days)", col = "red", border = "white")

# Censored times only
hist(data$time[data$event == 0], breaks = 20,
     main = paste("Censored Times (n =", sum(data$event == 0), ")"),
     xlab = "Time (days)", col = "green", border = "white")

# Kaplan-Meier curve (preliminary)
fit_all <- survfit(y ~ 1)
plot(fit_all, main = "Overall Survival", xlab = "Time (days)",
     ylab = "Survival Probability", col = "blue", lwd = 2)

# Time vs Event
plot(data$time, jitter(data$event, amount = 0.05),
     main = "Event Distribution Over Time",
     xlab = "Time (days)", ylab = "Event (0=censored, 1=death)",
     pch = 19, col = ifelse(data$event == 1, "red", "green"), cex = 0.8)

# Gene-wise expression boxplot (sample)
set.seed(123)
sample_genes <- sample(1:ncol(X), min(20, ncol(X)))
boxplot(X[, sample_genes], las = 2, main = "Expression by Gene (sample)",
        xlab = "", ylab = "Expression", col = "lightblue", border = "gray50",
        cex.axis = 0.7)

dev.off()

cat("Diagnostic plots saved.\n")

#### =====================================================
#### STEP 2: UNIVARIATE COX ANALYSIS
#### =====================================================

cat("\n========== STEP 2: UNIVARIATE COX ANALYSIS ==========\n")

univariate_results <- data.frame(
  gene = colnames(X),
  coef = NA,
  hr = NA,
  se_coef = NA,
  z = NA,
  pval = NA,
  lower_95 = NA,
  upper_95 = NA
)

cat("Running univariate Cox for", ncol(X), "genes...\n")

for(i in 1:ncol(X)) {
  tryCatch({
    fit <- coxph(y ~ X[, i])
    coef_summary <- summary(fit)$coefficients
    conf_int <- summary(fit)$conf.int
    
    univariate_results$coef[i] <- coef_summary[1]
    univariate_results$hr[i] <- conf_int[1]
    univariate_results$se_coef[i] <- coef_summary[3]
    univariate_results$z[i] <- coef_summary[4]
    univariate_results$pval[i] <- coef_summary[5]
    univariate_results$lower_95[i] <- conf_int[3]
    univariate_results$upper_95[i] <- conf_int[4]
  }, error = function(e) {
    # Skip genes that fail to converge
  })
}

# Remove genes with failed fits
univariate_results <- univariate_results[!is.na(univariate_results$pval), ]

# Sort by p-value
univariate_results <- univariate_results[order(univariate_results$pval), ]

# Add significance indicators
univariate_results$sig <- ifelse(univariate_results$pval < 0.001, "***",
                                 ifelse(univariate_results$pval < 0.01, "**",
                                        ifelse(univariate_results$pval < 0.05, "*",
                                               ifelse(univariate_results$pval < 0.10, ".", ""))))

cat("Genes with p < 0.05:", sum(univariate_results$pval < 0.05), "\n")
cat("Genes with p < 0.10:", sum(univariate_results$pval < 0.10), "\n")

write.csv(univariate_results,
          "results_TNBC/LASSO_Cox/02_univariate/univariate_cox_results.csv",
          row.names = FALSE)

#### Univariate volcano plot
tiff("results_TNBC/LASSO_Cox/02_univariate/univariate_volcano.tiff",
     width = 8, height = 7, units = "in", res = 600, compression = "lzw")

univariate_results$neg_log10_p <- -log10(univariate_results$pval)

ggplot(univariate_results, aes(x = coef, y = neg_log10_p)) +
  geom_point(aes(color = pval < 0.05), alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("gray50", "red"), 
                     labels = c("NS", "p < 0.05")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Univariate Cox Regression: Volcano Plot",
    subtitle = paste(sum(univariate_results$pval < 0.05), 
                     "genes with p < 0.05"),
    x = "Coefficient",
    y = "-log10(p-value)",
    color = "Significance"
  ) +
  theme(legend.position = "top")

dev.off()

#### Top genes forest plot
top_20_genes <- head(univariate_results, 20)

tiff("results_TNBC/LASSO_Cox/02_univariate/top20_forest_plot.tiff",
     width = 10, height = 8, units = "in", res = 600, compression = "lzw")

top_20_genes$gene <- factor(top_20_genes$gene, 
                            levels = rev(top_20_genes$gene))

ggplot(top_20_genes, aes(x = hr, y = gene)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = lower_95, xmax = upper_95), height = 0.3) +
  geom_point(size = 3, color = "blue") +
  scale_x_log10() +
  theme_bw(base_size = 12) +
  labs(
    title = "Top 20 Genes by Univariate Cox P-value",
    x = "Hazard Ratio (95% CI)",
    y = ""
  )

dev.off()

cat("Univariate analysis plots saved.\n")

#### =====================================================
#### STEP 3: CONSENSUS LASSO (100 ITERATIONS)
#### =====================================================

cat("\n========== STEP 3: CONSENSUS LASSO ANALYSIS ==========\n")
cat("Running", N_CONSENSUS_ITERATIONS, "iterations...\n")

gene_selection_matrix <- matrix(0, 
                                nrow = N_CONSENSUS_ITERATIONS, 
                                ncol = ncol(X))
colnames(gene_selection_matrix) <- colnames(X)

lambda_min_values <- numeric(N_CONSENSUS_ITERATIONS)
lambda_1se_values <- numeric(N_CONSENSUS_ITERATIONS)
n_genes_selected <- numeric(N_CONSENSUS_ITERATIONS)

pb <- txtProgressBar(min = 0, max = N_CONSENSUS_ITERATIONS, style = 3)

for(i in 1:N_CONSENSUS_ITERATIONS) {
  set.seed(123 + i)
  
  cvfit <- cv.glmnet(X, y, family = "cox", alpha = 1, 
                     nfolds = 10, standardize = TRUE)
  
  lambda_min_values[i] <- cvfit$lambda.min
  lambda_1se_values[i] <- cvfit$lambda.1se
  
  coef_min <- coef(cvfit, s = "lambda.min")
  selected_genes <- rownames(coef_min)[coef_min[, 1] != 0]
  n_genes_selected[i] <- length(selected_genes)
  
  if(length(selected_genes) > 0) {
    gene_selection_matrix[i, selected_genes] <- 1
  }
  
  setTxtProgressBar(pb, i)
}

close(pb)

#### Calculate selection frequencies
selection_frequency <- colSums(gene_selection_matrix) / N_CONSENSUS_ITERATIONS

consensus_summary <- data.frame(
  Gene = names(selection_frequency[selection_frequency > 0]),
  Selection_Frequency = selection_frequency[selection_frequency > 0]
) %>%
  arrange(desc(Selection_Frequency))

# Add univariate stats to consensus results
consensus_summary <- merge(consensus_summary, 
                           univariate_results[, c("gene", "coef", "hr", "pval")],
                           by.x = "Gene", by.y = "gene", all.x = TRUE)

consensus_summary <- consensus_summary %>%
  arrange(desc(Selection_Frequency))

write.csv(consensus_summary,
          "results_TNBC/LASSO_Cox/03_consensus/consensus_gene_selection.csv",
          row.names = FALSE)

cat("\nConsensus results:\n")
cat("Total unique genes selected:", nrow(consensus_summary), "\n")
cat("Genes selected ≥80% of time:", sum(consensus_summary$Selection_Frequency >= 0.80), "\n")
cat("Genes selected ≥60% of time:", sum(consensus_summary$Selection_Frequency >= 0.60), "\n")
cat("Genes selected ≥50% of time:", sum(consensus_summary$Selection_Frequency >= 0.50), "\n")
cat("Mean genes per iteration:", round(mean(n_genes_selected), 1), "\n")
cat("Range:", min(n_genes_selected), "-", max(n_genes_selected), "\n")

#### Consensus iteration summary
iteration_summary <- data.frame(
  Iteration = 1:N_CONSENSUS_ITERATIONS,
  N_Genes = n_genes_selected,
  Lambda_Min = lambda_min_values,
  Lambda_1SE = lambda_1se_values
)

write.csv(iteration_summary,
          "results_TNBC/LASSO_Cox/03_consensus/iteration_summary.csv",
          row.names = FALSE)

#### Plot: Consensus selection frequency
tiff("results_TNBC/LASSO_Cox/03_consensus/selection_frequency_barplot.tiff",
     width = 12, height = 8, units = "in", res = 600, compression = "lzw")

top_consensus <- head(consensus_summary, min(30, nrow(consensus_summary)))
top_consensus$Gene <- factor(top_consensus$Gene, 
                             levels = rev(top_consensus$Gene))

ggplot(top_consensus, aes(x = Gene, y = Selection_Frequency * 100)) +
  geom_bar(stat = "identity", 
           aes(fill = Selection_Frequency >= CONSENSUS_THRESHOLD)) +
  geom_hline(yintercept = CONSENSUS_THRESHOLD * 100, 
             linetype = "dashed", color = "red", size = 1) +
  coord_flip() +
  scale_fill_manual(values = c("gray70", "steelblue"),
                    labels = c(paste0("< ", CONSENSUS_THRESHOLD * 100, "%"),
                               paste0("≥ ", CONSENSUS_THRESHOLD * 100, "%"))) +
  theme_bw(base_size = 12) +
  labs(
    title = paste("Gene Selection Frequency Across", N_CONSENSUS_ITERATIONS, "Iterations"),
    subtitle = paste(sum(consensus_summary$Selection_Frequency >= CONSENSUS_THRESHOLD),
                     "genes selected ≥", CONSENSUS_THRESHOLD * 100, "% of time"),
    x = "Gene",
    y = "Selection Frequency (%)",
    fill = "Frequency"
  ) +
  theme(legend.position = "top")

dev.off()

#### Plot: Number of genes per iteration
tiff("results_TNBC/LASSO_Cox/03_consensus/genes_per_iteration.tiff",
     width = 10, height = 6, units = "in", res = 600, compression = "lzw")

par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

hist(n_genes_selected, breaks = 20, 
     main = "Distribution of Genes Selected",
     xlab = "Number of Genes", ylab = "Frequency",
     col = "lightblue", border = "white")
abline(v = mean(n_genes_selected), col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste("Mean =", round(mean(n_genes_selected), 1)),
       col = "red", lty = 2, lwd = 2, bty = "n")

plot(1:N_CONSENSUS_ITERATIONS, n_genes_selected, 
     type = "l", col = "blue", lwd = 1.5,
     main = "Genes Selected Over Iterations",
     xlab = "Iteration", ylab = "Number of Genes")
abline(h = mean(n_genes_selected), col = "red", lwd = 2, lty = 2)

dev.off()

#### Plot: Lambda stability
tiff("results_TNBC/LASSO_Cox/03_consensus/lambda_stability.tiff",
     width = 10, height = 6, units = "in", res = 600, compression = "lzw")

par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

hist(lambda_min_values, breaks = 30,
     main = "Distribution of Lambda.min",
     xlab = "Lambda.min", ylab = "Frequency",
     col = "lightgreen", border = "white")

plot(lambda_min_values, lambda_1se_values,
     pch = 19, col = alpha("blue", 0.5),
     main = "Lambda.min vs Lambda.1se",
     xlab = "Lambda.min", ylab = "Lambda.1se")
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)

dev.off()

#### =====================================================
#### STEP 4: PRE-FILTERED LASSO APPROACH
#### =====================================================

cat("\n========== STEP 4: PRE-FILTERED LASSO APPROACH ==========\n")

# Select top N genes by univariate p-value
top_genes <- head(univariate_results$gene, PREFILTER_TOP_N)
X_filtered <- X[, top_genes]

cat("Pre-filtered to", PREFILTER_TOP_N, "genes\n")
cat("P-value range:", 
    round(min(univariate_results$pval[1:PREFILTER_TOP_N]), 4), "to",
    round(max(univariate_results$pval[1:PREFILTER_TOP_N]), 4), "\n")

#### Run LASSO on filtered genes
set.seed(123)
cvfit_filtered <- cv.glmnet(X_filtered, y, family = "cox", 
                            alpha = 1, nfolds = 10, standardize = TRUE)

# Save CV curve for filtered LASSO
tiff("results_TNBC/LASSO_Cox/04_prefiltered/cv_curve_filtered.tiff",
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
plot(cvfit_filtered, main = paste("Cross-Validation Curve (Top", PREFILTER_TOP_N, "Genes)"))
dev.off()

# Extract selected genes
coef_filtered <- coef(cvfit_filtered, s = "lambda.min")
selected_filtered <- rownames(coef_filtered)[coef_filtered[, 1] != 0]

cat("Genes selected from filtered set (lambda.min):", length(selected_filtered), "\n")

#### Stability analysis on filtered genes
cat("Testing stability of pre-filtered LASSO...\n")

stability_filtered_matrix <- matrix(0, nrow = 50, ncol = length(top_genes))
colnames(stability_filtered_matrix) <- top_genes

for(i in 1:50) {
  set.seed(123 + i)
  cvfit_stab <- cv.glmnet(X_filtered, y, family = "cox", alpha = 1, nfolds = 10)
  coef_stab <- coef(cvfit_stab, s = "lambda.min")
  selected_stab <- rownames(coef_stab)[coef_stab[, 1] != 0]
  if(length(selected_stab) > 0) {
    stability_filtered_matrix[i, selected_stab] <- 1
  }
}

stability_filtered_freq <- colSums(stability_filtered_matrix) / 50

filtered_stability_df <- data.frame(
  Gene = names(stability_filtered_freq[stability_filtered_freq > 0]),
  Selection_Frequency = stability_filtered_freq[stability_filtered_freq > 0]
) %>%
  arrange(desc(Selection_Frequency))

write.csv(filtered_stability_df,
          "results_TNBC/LASSO_Cox/04_prefiltered/prefiltered_stability.csv",
          row.names = FALSE)

cat("Genes selected ≥60% in pre-filtered approach:", 
    sum(filtered_stability_df$Selection_Frequency >= 0.6), "\n")

#### =====================================================
#### STEP 5: FINAL MODEL SELECTION
#### =====================================================

cat("\n========== STEP 5: FINAL MODEL CONSTRUCTION ==========\n")

# Select final genes based on consensus threshold
final_genes <- consensus_summary$Gene[consensus_summary$Selection_Frequency >= CONSENSUS_THRESHOLD]

cat("Final gene signature (consensus ≥", CONSENSUS_THRESHOLD * 100, "%):", 
    length(final_genes), "genes\n")

if(length(final_genes) == 0) {
  cat("WARNING: No genes meet", CONSENSUS_THRESHOLD * 100, "% threshold.\n")
  cat("Using top 5 genes by selection frequency instead.\n")
  final_genes <- head(consensus_summary$Gene, 5)
}

cat("Final genes:", paste(final_genes, collapse = ", "), "\n")

# Fit final LASSO model with all data
X_final <- X[, final_genes, drop = FALSE]

set.seed(123)
cvfit_final <- cv.glmnet(X_final, y, family = "cox", 
                         alpha = 1, nfolds = 10, standardize = TRUE)

# Get final coefficients
final_model <- glmnet(X_final, y, family = "cox", 
                      alpha = 1, lambda = cvfit_final$lambda.min)

final_coefs <- coef(final_model, s = cvfit_final$lambda.min)
final_coefs_df <- data.frame(
  Gene = rownames(final_coefs)[final_coefs[, 1] != 0],
  Coefficient = final_coefs[final_coefs[, 1] != 0, 1]
)

final_coefs_df$HR <- exp(final_coefs_df$Coefficient)

# Add selection frequency and univariate stats
final_coefs_df <- merge(final_coefs_df, 
                        consensus_summary[, c("Gene", "Selection_Frequency")],
                        by = "Gene", all.x = TRUE)

final_coefs_df <- merge(final_coefs_df,
                        univariate_results[, c("gene", "pval", "lower_95", "upper_95")],
                        by.x = "Gene", by.y = "gene", all.x = TRUE)

final_coefs_df <- final_coefs_df %>%
  arrange(desc(abs(Coefficient)))

write.csv(final_coefs_df,
          "results_TNBC/LASSO_Cox/05_final_model/final_gene_signature.csv",
          row.names = FALSE)

cat("\nFinal gene signature:\n")
print(final_coefs_df)

#### Final coefficient plot
tiff("results_TNBC/LASSO_Cox/05_final_model/final_coefficients.tiff",
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")

final_coefs_df$Gene <- factor(final_coefs_df$Gene,
                              levels = final_coefs_df$Gene[order(final_coefs_df$Coefficient)])

ggplot(final_coefs_df, aes(x = Gene, y = Coefficient, fill = Coefficient > 0)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("steelblue", "firebrick"),
                    labels = c("Protective", "Risk")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Final LASSO-Cox Gene Signature",
    subtitle = paste(nrow(final_coefs_df), "genes"),
    x = "",
    y = "LASSO Coefficient",
    fill = "Direction"
  ) +
  theme(legend.position = "top")

dev.off()

#### =====================================================
#### STEP 6: RISK SCORE CALCULATION
#### =====================================================

cat("\n========== STEP 6: RISK SCORE & STRATIFICATION ==========\n")

# Calculate risk scores
risk_score <- predict(final_model, newx = X_final, 
                      s = cvfit_final$lambda.min, type = "link")
data$risk_score <- as.numeric(risk_score)

# Risk group assignment (median split)
median_risk <- median(data$risk_score)
data$risk_group <- ifelse(data$risk_score > median_risk, "High", "Low")
data$risk_group <- factor(data$risk_group, levels = c("Low", "High"))

cat("Risk score range:", round(min(data$risk_score), 3), "to",
    round(max(data$risk_score), 3), "\n")
cat("Median risk score:", round(median_risk, 3), "\n")
cat("High risk group:", sum(data$risk_group == "High"), "samples\n")
cat("Low risk group:", sum(data$risk_group == "Low"), "samples\n")

# Events by risk group
risk_table <- table(data$risk_group, data$event)
cat("\nEvents by risk group:\n")
print(risk_table)
cat("High risk event rate:", 
    round(100 * risk_table["High", "1"] / sum(risk_table["High", ]), 1), "%\n")
cat("Low risk event rate:", 
    round(100 * risk_table["Low", "1"] / sum(risk_table["Low", ]), 1), "%\n")

# Save risk scores
risk_scores_df <- data.frame(
  Sample_ID = data$Sample_ID,
  Risk_Score = data$risk_score,
  Risk_Group = data$risk_group,
  Time = data$time,
  Event = data$event
)

write.csv(risk_scores_df,
          "results_TNBC/LASSO_Cox/05_final_model/risk_scores.csv",
          row.names = FALSE)

#### Kaplan-Meier survival analysis
fit_risk <- survfit(Surv(time, event) ~ risk_group, data = data)

# Log-rank test
logrank_test <- survdiff(Surv(time, event) ~ risk_group, data = data)
logrank_pval <- 1 - pchisq(logrank_test$chisq, 1)

cat("\nLog-rank test p-value:", format.pval(logrank_pval, digits = 3), "\n")

# Median survival times
median_survival <- summary(fit_risk)$table[, "median"]
cat("Median survival - High risk:", median_survival["risk_group=High"], "days\n")
cat("Median survival - Low risk:", median_survival["risk_group=Low"], "days\n")

#### Kaplan-Meier plot
tiff("results_TNBC/LASSO_Cox/05_final_model/kaplan_meier_curve.tiff",
     width = 9, height = 8, units = "in", res = 600, compression = "lzw")

km_plot <- ggsurvplot(
  fit_risk,
  data = data,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  palette = c("#4DBBD5", "#E64B35"),
  title = "Kaplan-Meier Survival Curves by Risk Group",
  xlab = "Time (days)",
  ylab = "Survival Probability",
  legend.title = "Risk Group",
  legend.labs = c("Low Risk", "High Risk"),
  font.main = c(16, "bold"),
  font.x = c(14, "plain"),
  font.y = c(14, "plain"),
  font.tickslab = c(12, "plain"),
  ggtheme = theme_bw(base_size = 12)
)

print(km_plot)
dev.off()

#### Risk score distribution plot
tiff("results_TNBC/LASSO_Cox/05_final_model/risk_score_distribution.tiff",
     width = 10, height = 6, units = "in", res = 600, compression = "lzw")

par(mfrow = c(1, 2))

# Histogram
hist(data$risk_score, breaks = 30, 
     main = "Risk Score Distribution",
     xlab = "Risk Score", ylab = "Frequency",
     col = "lightblue", border = "white")
abline(v = median_risk, col = "red", lwd = 2, lty = 2)
legend("topright", legend = "Median", col = "red", lty = 2, lwd = 2, bty = "n")

# Boxplot by risk group and event
boxplot(risk_score ~ risk_group + event, data = data,
        names = c("Low/Censored", "High/Censored", "Low/Event", "High/Event"),
        col = c("#4DBBD5", "#E64B35", "#4DBBD5", "#E64B35"),
        main = "Risk Score by Group and Outcome",
        xlab = "", ylab = "Risk Score", las = 2)
abline(h = median_risk, col = "red", lwd = 2, lty = 2)

dev.off()

#### =====================================================
#### STEP 7: BOOTSTRAP VALIDATION
#### =====================================================

cat("\n========== STEP 7: BOOTSTRAP VALIDATION ==========\n")
cat("Running", N_BOOTSTRAP_VALIDATION, "bootstrap iterations...\n")

# Apparent C-index
cox_apparent <- coxph(Surv(time, event) ~ risk_score, data = data)
c_apparent <- summary(cox_apparent)$concordance[1]

cat("Apparent C-index:", round(c_apparent, 3), "\n")

# Bootstrap validation
c_bootstrap <- numeric(N_BOOTSTRAP_VALIDATION)
c_original <- numeric(N_BOOTSTRAP_VALIDATION)

pb <- txtProgressBar(min = 0, max = N_BOOTSTRAP_VALIDATION, style = 3)

for(i in 1:N_BOOTSTRAP_VALIDATION) {
  set.seed(456 + i)
  
  # Bootstrap sample
  boot_idx <- sample(1:nrow(X_final), replace = TRUE)
  X_boot <- X_final[boot_idx, , drop = FALSE]
  y_boot <- y[boot_idx]
  
  # Fit model on bootstrap
  tryCatch({
    cvfit_boot <- cv.glmnet(X_boot, y_boot, family = "cox", 
                            alpha = 1, nfolds = 10)
    
    # C-index on bootstrap sample
    risk_boot <- predict(cvfit_boot, newx = X_boot, 
                         s = "lambda.min", type = "link")
    cox_boot <- coxph(y_boot ~ risk_boot)
    c_bootstrap[i] <- summary(cox_boot)$concordance[1]
    
    # C-index on original sample (using bootstrap model)
    risk_orig <- predict(cvfit_boot, newx = X_final, 
                         s = "lambda.min", type = "link")
    cox_orig <- coxph(y ~ risk_orig)
    c_original[i] <- summary(cox_orig)$concordance[1]
  }, error = function(e) {
    c_bootstrap[i] <<- NA
    c_original[i] <<- NA
  })
  
  setTxtProgressBar(pb, i)
}

close(pb)

# Remove failed bootstraps
valid_boots <- !is.na(c_bootstrap) & !is.na(c_original)
c_bootstrap <- c_bootstrap[valid_boots]
c_original <- c_original[valid_boots]

# Calculate optimism
optimism <- mean(c_bootstrap - c_original)
c_corrected <- c_apparent - optimism

# 95% CI for corrected C-index
c_bootstrap_corrected <- c_bootstrap - (c_bootstrap - c_original)
ci_lower <- quantile(c_bootstrap_corrected, 0.025)
ci_upper <- quantile(c_bootstrap_corrected, 0.975)

cat("\n--- Bootstrap Validation Results ---\n")
cat("Successful bootstraps:", sum(valid_boots), "/", N_BOOTSTRAP_VALIDATION, "\n")
cat("Apparent C-index:", round(c_apparent, 3), "\n")
cat("Mean optimism:", round(optimism, 3), "\n")
cat("Optimism-corrected C-index:", round(c_corrected, 3), "\n")
cat("95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3), "\n")

# Save validation results
validation_results <- data.frame(
  Metric = c("Apparent_C_Index", "Mean_Optimism", 
             "Corrected_C_Index", "CI_Lower", "CI_Upper"),
  Value = c(c_apparent, optimism, c_corrected, ci_lower, ci_upper)
)

write.csv(validation_results,
          "results_TNBC/LASSO_Cox/06_validation/bootstrap_validation.csv",
          row.names = FALSE)

# Save all bootstrap values
bootstrap_values <- data.frame(
  Bootstrap_ID = 1:length(c_bootstrap),
  C_Bootstrap = c_bootstrap,
  C_Original = c_original,
  Optimism = c_bootstrap - c_original
)

write.csv(bootstrap_values,
          "results_TNBC/LASSO_Cox/06_validation/bootstrap_values.csv",
          row.names = FALSE)

#### Bootstrap validation plots
tiff("results_TNBC/LASSO_Cox/06_validation/bootstrap_validation.tiff",
     width = 12, height = 6, units = "in", res = 600, compression = "lzw")

par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

# Histogram of optimism
hist(c_bootstrap - c_original, breaks = 30,
     main = "Distribution of Optimism",
     xlab = "Optimism (Bootstrap - Original)",
     ylab = "Frequency",
     col = "lightcoral", border = "white")
abline(v = optimism, col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste("Mean =", round(optimism, 3)),
       col = "red", lty = 2, lwd = 2, bty = "n")

# Bootstrap vs original C-index
plot(c_original, c_bootstrap,
     pch = 19, col = alpha("blue", 0.3),
     main = "Bootstrap vs Original C-index",
     xlab = "C-index on Original Data",
     ylab = "C-index on Bootstrap Data",
     xlim = c(0.4, 1), ylim = c(0.4, 1))
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)

# C-index distribution
hist(c_bootstrap_corrected, breaks = 30,
     main = "Optimism-Corrected C-index Distribution",
     xlab = "Corrected C-index",
     ylab = "Frequency",
     col = "lightgreen", border = "white")
abline(v = c_corrected, col = "red", lwd = 2, lty = 2)
abline(v = c(ci_lower, ci_upper), col = "blue", lwd = 2, lty = 3)
legend("topleft", 
       legend = c(paste("Mean =", round(c_corrected, 3)),
                  paste("95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3))),
       col = c("red", "blue"), lty = c(2, 3), lwd = 2, bty = "n", cex = 0.9)

dev.off()

#### =====================================================
#### STEP 8: COMPREHENSIVE SUMMARY REPORT
#### =====================================================

cat("\n========== STEP 8: GENERATING SUMMARY REPORT ==========\n")

# Create comprehensive summary
summary_report <- list(
  Dataset = list(
    Total_Samples = nrow(data),
    Total_Genes = ncol(X),
    Events = sum(data$event),
    Censored = sum(data$event == 0),
    Event_Rate = round(100 * sum(data$event) / nrow(data), 1),
    EPV = round(sum(data$event) / ncol(X), 3)
  ),
  
  Univariate = list(
    Genes_p_less_0.05 = sum(univariate_results$pval < 0.05),
    Genes_p_less_0.10 = sum(univariate_results$pval < 0.10),
    Top_Gene = univariate_results$gene[1],
    Top_Gene_Pval = univariate_results$pval[1]
  ),
  
  Consensus = list(
    Iterations = N_CONSENSUS_ITERATIONS,
    Threshold = CONSENSUS_THRESHOLD,
    Total_Unique_Genes = nrow(consensus_summary),
    Genes_Above_Threshold = sum(consensus_summary$Selection_Frequency >= CONSENSUS_THRESHOLD),
    Mean_Genes_Per_Iteration = round(mean(n_genes_selected), 1),
    Range_Genes = paste(min(n_genes_selected), "-", max(n_genes_selected))
  ),
  
  Final_Model = list(
    N_Genes = nrow(final_coefs_df),
    Genes = paste(final_coefs_df$Gene, collapse = ", "),
    Lambda_Min = round(cvfit_final$lambda.min, 4)
  ),
  
  Risk_Stratification = list(
    High_Risk_N = sum(data$risk_group == "High"),
    Low_Risk_N = sum(data$risk_group == "Low"),
    High_Risk_Events = sum(data$risk_group == "High" & data$event == 1),
    Low_Risk_Events = sum(data$risk_group == "Low" & data$event == 1),
    LogRank_Pval = format.pval(logrank_pval, digits = 3)
  ),
  
  Performance = list(
    Apparent_C_Index = round(c_apparent, 3),
    Optimism = round(optimism, 3),
    Corrected_C_Index = round(c_corrected, 3),
    CI_95 = paste0(round(ci_lower, 3), " - ", round(ci_upper, 3))
  )
)

# Save as JSON-like text file
capture.output(
  summary_report,
  file = "results_TNBC/LASSO_Cox/ANALYSIS_SUMMARY.txt"
)

# Save as CSV for easy reading
summary_flat <- data.frame(
  Category = c("Dataset", "Dataset", "Dataset", "Dataset", "Dataset", "Dataset",
               "Univariate", "Univariate", "Univariate", "Univariate",
               "Consensus", "Consensus", "Consensus", "Consensus", "Consensus", "Consensus",
               "Final_Model", "Final_Model", "Final_Model",
               "Risk_Stratification", "Risk_Stratification", "Risk_Stratification", 
               "Risk_Stratification", "Risk_Stratification",
               "Performance", "Performance", "Performance", "Performance"),
  Metric = c("Total_Samples", "Total_Genes", "Events", "Censored", "Event_Rate_%", "EPV",
             "Genes_p<0.05", "Genes_p<0.10", "Top_Gene", "Top_Gene_Pval",
             "Iterations", "Threshold", "Total_Unique_Genes", "Genes_Above_Threshold",
             "Mean_Genes_Per_Iter", "Range_Genes",
             "N_Genes", "Genes", "Lambda_Min",
             "High_Risk_N", "Low_Risk_N", "High_Risk_Events", "Low_Risk_Events", "LogRank_Pval",
             "Apparent_C_Index", "Optimism", "Corrected_C_Index", "95%_CI"),
  Value = c(nrow(data), ncol(X), sum(data$event), sum(data$event == 0),
            round(100 * sum(data$event) / nrow(data), 1), round(sum(data$event) / ncol(X), 3),
            sum(univariate_results$pval < 0.05), sum(univariate_results$pval < 0.10),
            as.character(univariate_results$gene[1]), univariate_results$pval[1],
            N_CONSENSUS_ITERATIONS, CONSENSUS_THRESHOLD, nrow(consensus_summary),
            sum(consensus_summary$Selection_Frequency >= CONSENSUS_THRESHOLD),
            round(mean(n_genes_selected), 1), 
            paste(min(n_genes_selected), "-", max(n_genes_selected)),
            nrow(final_coefs_df), paste(final_coefs_df$Gene, collapse = ", "),
            round(cvfit_final$lambda.min, 4),
            sum(data$risk_group == "High"), sum(data$risk_group == "Low"),
            sum(data$risk_group == "High" & data$event == 1),
            sum(data$risk_group == "Low" & data$event == 1),
            format.pval(logrank_pval, digits = 3),
            round(c_apparent, 3), round(optimism, 3), round(c_corrected, 3),
            paste0(round(ci_lower, 3), " - ", round(ci_upper, 3)))
)

write.csv(summary_flat,
          "results_TNBC/LASSO_Cox/ANALYSIS_SUMMARY.csv",
          row.names = FALSE)

#### =====================================================
#### FINAL OUTPUTS
#### =====================================================

cat("\n==========================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("==========================================================\n\n")

cat("Key Findings:\n")
cat("-------------\n")
cat("• Dataset:", nrow(data), "samples,", ncol(X), "genes,", 
    sum(data$event), "events\n")
cat("• EPV ratio:", round(sum(data$event) / ncol(X), 3), "(underpowered)\n")
cat("• Consensus genes (≥", CONSENSUS_THRESHOLD * 100, "%):", 
    sum(consensus_summary$Selection_Frequency >= CONSENSUS_THRESHOLD), "\n")
cat("• Final signature:", nrow(final_coefs_df), "genes\n")
cat("• Log-rank p-value:", format.pval(logrank_pval, digits = 3), "\n")
cat("• Optimism-corrected C-index:", round(c_corrected, 3), 
    "(95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3), ")\n\n")

cat("All results saved to: results_TNBC/LASSO_Cox/\n")
cat("\nKey output files:\n")
cat("  01_diagnostics/     - Data quality checks\n")
cat("  02_univariate/      - Univariate Cox results\n")
cat("  03_consensus/       - Consensus gene selection\n")
cat("  04_prefiltered/     - Pre-filtered LASSO approach\n")
cat("  05_final_model/     - Final signature & KM curves\n")
cat("  06_validation/      - Bootstrap validation\n")
cat("  ANALYSIS_SUMMARY.csv - Complete summary report\n\n")

cat("Analysis completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==========================================================\n")

sink()
#### =====================================================
#### SAVE R OBJECTS (RDS FORMAT)
#### =====================================================

cat("\n========== SAVING RDS OBJECTS ==========\n")

# Create directory for RDS files
dir.create("results_TNBC/LASSO_Cox/RDS_objects", showWarnings = FALSE)

# 1️⃣ Final cross-validated LASSO object
saveRDS(cvfit_final,
        file = "results_TNBC/LASSO_Cox/RDS_objects/cvfit_final.rds")

# 2️⃣ Final glmnet model (lambda.min)
saveRDS(final_model,
        file = "results_TNBC/LASSO_Cox/RDS_objects/final_glmnet_model.rds")

# 3️⃣ Final gene signature table
saveRDS(final_coefs_df,
        file = "results_TNBC/LASSO_Cox/RDS_objects/final_gene_signature_df.rds")

# 4️⃣ Consensus selection matrix
saveRDS(gene_selection_matrix,
        file = "results_TNBC/LASSO_Cox/RDS_objects/consensus_selection_matrix.rds")

# 5️⃣ Consensus summary
saveRDS(consensus_summary,
        file = "results_TNBC/LASSO_Cox/RDS_objects/consensus_summary.rds")

# 6️⃣ Bootstrap validation values
saveRDS(list(
  c_apparent = c_apparent,
  c_bootstrap = c_bootstrap,
  c_original = c_original,
  c_corrected = c_corrected,
  ci_lower = ci_lower,
  ci_upper = ci_upper
),
file = "results_TNBC/LASSO_Cox/RDS_objects/bootstrap_validation_results.rds")

# 7️⃣ Entire analysis workspace (Highly Recommended)
save.image("results_TNBC/LASSO_Cox/RDS_objects/full_workspace.RData")
model_objects <- readRDS("results_TNBC/LASSO_Cox/08_model_objects/complete_model_objects.rds")

final_genes <- model_objects$final_genes
lasso_model <- model_objects$final_model
cox_model <- model_objects$cox_model
median_risk <- median(model_objects$risk_scores)


cat("All RDS files saved successfully.\n")

cat("\n✅ ANALYSIS COMPLETE! Check results_TNBC/LASSO_Cox/ for all outputs.\n")
