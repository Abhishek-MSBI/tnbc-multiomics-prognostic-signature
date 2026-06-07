#### =====================================================
#### TIME-DEPENDENT ROC ANALYSIS FOR TNBC LASSO-COX MODEL
#### 
#### Performs time-dependent ROC analysis at 1, 3, and 5 years
#### Calculates AUC with confidence intervals
#### Generates publication-quality plots (600 dpi TIFF)
#### Saves model objects in RDS format for reproducibility
#### 
#### =====================================================

#### Load required packages
required_packages <- c("tidyverse", "glmnet", "survival", "survminer", 
                       "timeROC", "survivalROC", "pROC", "gridExtra")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

#### Create output directories
dir.create("results_TNBC/LASSO_Cox/07_timeROC", recursive = TRUE, showWarnings = FALSE)
dir.create("results_TNBC/LASSO_Cox/08_model_objects", recursive = TRUE, showWarnings = FALSE)

#### Start logging
sink("results_TNBC/LASSO_Cox/timeROC_analysis_log.txt", split = TRUE)
cat("\n==========================================================\n")
cat("TIME-DEPENDENT ROC ANALYSIS FOR TNBC PROGNOSTIC SIGNATURE\n")
cat("Analysis started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==========================================================\n\n")

#### =====================================================
#### STEP 1: LOAD DATA AND PREVIOUS RESULTS
#### =====================================================

cat("\n========== STEP 1: LOADING DATA ==========\n")

# Check if previous results exist
if(file.exists("results_TNBC/LASSO_Cox/05_final_model/risk_scores.csv")) {
  
  cat("Loading previous analysis results...\n")
  
  # Load data
  expr <- read.csv("tcga_expression_for_lasso_cox.csv", stringsAsFactors = FALSE)
  surv <- read.csv("tcga_survival_for_lasso_cox.csv", stringsAsFactors = FALSE)
  
  expr <- expr[!duplicated(expr$Sample_ID), ]
  data <- merge(expr, surv, by = "Sample_ID")
  
  # Load risk scores
  risk_scores <- read.csv("results_TNBC/LASSO_Cox/05_final_model/risk_scores.csv")
  
  # Load final gene signature
  final_genes_df <- read.csv("results_TNBC/LASSO_Cox/05_final_model/final_gene_signature.csv")
  final_genes <- final_genes_df$Gene
  
  cat("Loaded:", nrow(risk_scores), "samples\n")
  cat("Final signature:", length(final_genes), "genes\n")
  
  # Merge risk scores with data
  data <- merge(data, risk_scores[, c("Sample_ID", "Risk_Score", "Risk_Group")], 
                by = "Sample_ID")
  
  # Create expression matrix for final genes
  X <- as.matrix(data[, final_genes])
  y <- Surv(time = data$time, event = data$event)
  
} else {
  
  cat("ERROR: Previous analysis results not found!\n")
  cat("Please run the main LASSO-Cox analysis script first.\n")
  stop("Previous results required for time-ROC analysis")
  
}

cat("Samples:", nrow(data), "\n")
cat("Events:", sum(data$event), "\n")
cat("Follow-up range:", round(min(data$time), 1), "to", 
    round(max(data$time), 1), "days\n")

# Convert days to years for easier interpretation
data$time_years <- data$time / 365.25

cat("Follow-up range (years):", round(min(data$time_years), 2), "to",
    round(max(data$time_years), 2), "\n\n")

#### =====================================================
#### STEP 2: DEFINE TIME POINTS FOR ROC ANALYSIS
#### =====================================================

cat("\n========== STEP 2: DEFINING TIME POINTS ==========\n")

# Time points in days (1, 3, 5 years)
time_points_years <- c(1, 3, 5)
time_points_days <- time_points_years * 365.25

cat("Time points for ROC analysis:\n")
for(i in 1:length(time_points_years)) {
  n_events_by_time <- sum(data$event == 1 & data$time <= time_points_days[i])
  n_at_risk <- sum(data$time >= time_points_days[i])
  cat(sprintf("  %d-year: %.1f days (%d events by this time, %d still at risk)\n",
              time_points_years[i], time_points_days[i], 
              n_events_by_time, n_at_risk))
}

# Check if we have enough follow-up for each time point
max_followup_years <- max(data$time_years)
cat("\nMaximum follow-up:", round(max_followup_years, 2), "years\n")

# Adjust time points if necessary
valid_time_points <- time_points_years[time_points_years <= max_followup_years]
valid_time_days <- valid_time_points * 365.25

if(length(valid_time_points) < length(time_points_years)) {
  cat("WARNING: Insufficient follow-up for some time points.\n")
  cat("Using time points:", paste(valid_time_points, "years"), "\n")
  time_points_years <- valid_time_points
  time_points_days <- valid_time_days
}

cat("\n")

#### =====================================================
#### STEP 3: TIME-DEPENDENT ROC USING timeROC PACKAGE
#### =====================================================

cat("\n========== STEP 3: TIME-DEPENDENT ROC ANALYSIS ==========\n")

# Prepare data for timeROC
# Note: timeROC requires complete cases
complete_cases <- complete.cases(data$time, data$event, data$Risk_Score)
data_complete <- data[complete_cases, ]

cat("Complete cases for analysis:", nrow(data_complete), "\n")

# Run timeROC analysis
cat("Running timeROC analysis...\n")

set.seed(123)

tryCatch({
  
  roc_obj <- timeROC(
    T = data_complete$time,
    delta = data_complete$event,
    marker = data_complete$Risk_Score,
    cause = 1,  # Event of interest
    weighting = "marginal",  # Marginal weighting (recommended)
    times = time_points_days,
    iid = TRUE  # Calculate confidence intervals
  )
  
  cat("✓ timeROC analysis completed successfully\n\n")
  
  # Extract AUC values and CIs
  auc_results <- data.frame(
    Time_Point_Years = time_points_years,
    Time_Point_Days = time_points_days,
    AUC = roc_obj$AUC,
    SE = sqrt(roc_obj$inference$vect_sd_1^2),
    CI_Lower = roc_obj$AUC - 1.96 * sqrt(roc_obj$inference$vect_sd_1^2),
    CI_Upper = roc_obj$AUC + 1.96 * sqrt(roc_obj$inference$vect_sd_1^2)
  )
  
  cat("Time-Dependent AUC Results:\n")
  cat("---------------------------\n")
  for(i in 1:nrow(auc_results)) {
    cat(sprintf("%d-year AUC: %.3f (95%% CI: %.3f - %.3f)\n",
                auc_results$Time_Point_Years[i],
                auc_results$AUC[i],
                auc_results$CI_Lower[i],
                auc_results$CI_Upper[i]))
  }
  cat("\n")
  
  # Save AUC results
  write.csv(auc_results,
            "results_TNBC/LASSO_Cox/07_timeROC/timeROC_AUC_summary.csv",
            row.names = FALSE)
  
  # Save complete ROC object
  saveRDS(roc_obj, "results_TNBC/LASSO_Cox/08_model_objects/timeROC_object.rds")
  
}, error = function(e) {
  cat("ERROR in timeROC analysis:", e$message, "\n")
  cat("Attempting alternative method...\n\n")
  roc_obj <<- NULL
})

#### =====================================================
#### STEP 4: ALTERNATIVE ROC USING survivalROC PACKAGE
#### =====================================================

cat("\n========== STEP 4: ALTERNATIVE ROC ANALYSIS (survivalROC) ==========\n")

# survivalROC uses a different approach - calculate for each time point separately
survival_roc_results <- list()

for(i in 1:length(time_points_years)) {
  
  cat(sprintf("Computing %d-year ROC...\n", time_points_years[i]))
  
  tryCatch({
    
    sroc <- survivalROC(
      Stime = data_complete$time,
      status = data_complete$event,
      marker = data_complete$Risk_Score,
      predict.time = time_points_days[i],
      method = "KM"  # Kaplan-Meier method
    )
    
    survival_roc_results[[i]] <- sroc
    
    cat(sprintf("  AUC: %.3f\n", sroc$AUC))
    
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    survival_roc_results[[i]] <<- NULL
  })
}

cat("\n")

# Create summary table
surv_auc_summary <- data.frame(
  Time_Point_Years = time_points_years,
  Time_Point_Days = time_points_days,
  AUC_survivalROC = sapply(survival_roc_results, function(x) {
    if(!is.null(x)) x$AUC else NA
  })
)

write.csv(surv_auc_summary,
          "results_TNBC/LASSO_Cox/07_timeROC/survivalROC_AUC_summary.csv",
          row.names = FALSE)

saveRDS(survival_roc_results, 
        "results_TNBC/LASSO_Cox/08_model_objects/survivalROC_objects.rds")

#### =====================================================
#### STEP 5: TIME-DEPENDENT ROC PLOTS
#### =====================================================

cat("\n========== STEP 5: GENERATING ROC PLOTS ==========\n")

#### Plot 1: Combined time-dependent ROC curves
tiff("results_TNBC/LASSO_Cox/07_timeROC/time_dependent_ROC_curves.tiff",
     width = 10, height = 8, units = "in", res = 600, compression = "lzw")

if(!is.null(roc_obj)) {
  
  # Create color palette
  colors <- c("#E64B35", "#4DBBD5", "#00A087")
  
  # Plot using base R for better control
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate (1 - Specificity)",
       ylab = "True Positive Rate (Sensitivity)",
       main = "Time-Dependent ROC Curves",
       cex.lab = 1.3, cex.main = 1.5, cex.axis = 1.2)
  
  # Add diagonal reference line
  abline(a = 0, b = 1, lty = 2, col = "gray50", lwd = 2)
  
  # Plot ROC curves for each time point
  for(i in 1:length(time_points_years)) {
    lines(roc_obj$FP[, i], roc_obj$TP[, i], 
          col = colors[i], lwd = 3, type = "l")
  }
  
  # Add legend with AUC values
  legend_text <- sapply(1:length(time_points_years), function(i) {
    sprintf("%d-year AUC = %.3f (%.3f-%.3f)",
            time_points_years[i],
            auc_results$AUC[i],
            auc_results$CI_Lower[i],
            auc_results$CI_Upper[i])
  })
  
  legend("bottomright", 
         legend = legend_text,
         col = colors[1:length(time_points_years)],
         lwd = 3,
         bty = "n",
         cex = 1.1)
  
  # Add grid
  grid(col = "gray90", lty = 1)
  
} else {
  plot.new()
  text(0.5, 0.5, "timeROC analysis not available", cex = 1.5)
}

dev.off()
cat("✓ Combined ROC curves saved\n")

#### Plot 2: Individual ROC curves (survivalROC method)
tiff("results_TNBC/LASSO_Cox/07_timeROC/individual_ROC_curves.tiff",
     width = 12, height = 4, units = "in", res = 600, compression = "lzw")

par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

for(i in 1:length(time_points_years)) {
  
  if(!is.null(survival_roc_results[[i]])) {
    
    sroc <- survival_roc_results[[i]]
    
    plot(sroc$FP, sroc$TP, type = "l", lwd = 3, col = "#E64B35",
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = "1 - Specificity",
         ylab = "Sensitivity",
         main = sprintf("%d-Year ROC Curve", time_points_years[i]),
         cex.lab = 1.2, cex.main = 1.3, cex.axis = 1.1)
    
    abline(a = 0, b = 1, lty = 2, col = "gray50", lwd = 2)
    
    # Add AUC text
    text(0.6, 0.2, 
         sprintf("AUC = %.3f", sroc$AUC),
         cex = 1.3, font = 2)
    
    grid(col = "gray90", lty = 1)
    
  } else {
    plot.new()
    text(0.5, 0.5, "Analysis not available", cex = 1.2)
  }
}

dev.off()
cat("✓ Individual ROC curves saved\n")

#### Plot 3: AUC comparison across time points
tiff("results_TNBC/LASSO_Cox/07_timeROC/AUC_over_time.tiff",
     width = 10, height = 7, units = "in", res = 600, compression = "lzw")

if(!is.null(roc_obj)) {
  
  par(mar = c(5, 5, 4, 2))
  
  # Create plot
  plot(time_points_years, auc_results$AUC,
       type = "b", pch = 19, cex = 2, lwd = 2, col = "#E64B35",
       ylim = c(0.5, 1),
       xlab = "Time (years)",
       ylab = "AUC",
       main = "Time-Dependent AUC with 95% Confidence Intervals",
       cex.lab = 1.3, cex.main = 1.5, cex.axis = 1.2)
  
  # Add error bars (95% CI)
  arrows(time_points_years, auc_results$CI_Lower,
         time_points_years, auc_results$CI_Upper,
         angle = 90, code = 3, length = 0.1, lwd = 2, col = "#E64B35")
  
  # Add reference line at AUC = 0.5
  abline(h = 0.5, lty = 2, col = "gray50", lwd = 2)
  
  # Add reference line at AUC = 0.7 (acceptable)
  abline(h = 0.7, lty = 3, col = "blue", lwd = 1.5)
  
  # Add text labels
  text(time_points_years, auc_results$AUC + 0.05,
       sprintf("%.3f", auc_results$AUC),
       cex = 1.1, font = 2)
  
  # Add legend
  legend("bottomleft",
         legend = c("AUC with 95% CI", "Random classifier", "Acceptable (0.7)"),
         col = c("#E64B35", "gray50", "blue"),
         lty = c(1, 2, 3),
         lwd = c(2, 2, 1.5),
         pch = c(19, NA, NA),
         bty = "n",
         cex = 1.1)
  
  grid(col = "gray90", lty = 1)
  
} else {
  plot.new()
  text(0.5, 0.5, "timeROC analysis not available", cex = 1.5)
}

dev.off()
cat("✓ AUC over time plot saved\n")

#### Plot 4: Sensitivity and Specificity over time
tiff("results_TNBC/LASSO_Cox/07_timeROC/sensitivity_specificity_time.tiff",
     width = 10, height = 7, units = "in", res = 600, compression = "lzw")

if(!is.null(roc_obj)) {
  
  # Find optimal cutoff for each time point (Youden's index)
  optimal_cutoffs <- data.frame(
    Time_Years = time_points_years,
    Cutoff = NA,
    Sensitivity = NA,
    Specificity = NA,
    Youden = NA
  )
  
  for(i in 1:length(time_points_years)) {
    
    # Calculate Youden's index for each threshold
    youden <- roc_obj$TP[, i] + (1 - roc_obj$FP[, i]) - 1
    max_idx <- which.max(youden)
    
    # Get the marker value at this threshold
    # Note: timeROC doesn't directly provide marker values, 
    # so we'll use the index position
    optimal_cutoffs$Sensitivity[i] <- roc_obj$TP[max_idx, i]
    optimal_cutoffs$Specificity[i] <- 1 - roc_obj$FP[max_idx, i]
    optimal_cutoffs$Youden[i] <- youden[max_idx]
  }
  
  par(mar = c(5, 5, 4, 2))
  
  plot(time_points_years, optimal_cutoffs$Sensitivity,
       type = "b", pch = 19, cex = 2, lwd = 2, col = "#E64B35",
       ylim = c(0, 1),
       xlab = "Time (years)",
       ylab = "Proportion",
       main = "Optimal Sensitivity and Specificity over Time",
       cex.lab = 1.3, cex.main = 1.5, cex.axis = 1.2)
  
  lines(time_points_years, optimal_cutoffs$Specificity,
        type = "b", pch = 17, cex = 2, lwd = 2, col = "#4DBBD5")
  
  legend("bottomleft",
         legend = c("Sensitivity", "Specificity"),
         col = c("#E64B35", "#4DBBD5"),
         pch = c(19, 17),
         lty = 1,
         lwd = 2,
         bty = "n",
         cex = 1.2)
  
  grid(col = "gray90", lty = 1)
  
  # Save optimal cutoffs
  write.csv(optimal_cutoffs,
            "results_TNBC/LASSO_Cox/07_timeROC/optimal_cutoffs.csv",
            row.names = FALSE)
  
} else {
  plot.new()
  text(0.5, 0.5, "timeROC analysis not available", cex = 1.5)
}

dev.off()
cat("✓ Sensitivity/Specificity plot saved\n")

#### =====================================================
#### STEP 6: CALIBRATION ANALYSIS
#### =====================================================

cat("\n========== STEP 6: CALIBRATION ANALYSIS ==========\n")

# Calculate predicted probabilities at each time point
cat("Calculating predicted survival probabilities...\n")

# Fit Cox model with risk score
cox_model <- coxph(Surv(time, event) ~ Risk_Score, data = data_complete)

# Create calibration data
calibration_results <- list()

for(i in 1:length(time_points_years)) {
  
  cat(sprintf("Calibration at %d years...\n", time_points_years[i]))
  
  # Calculate baseline survival at this time
  base_surv <- survfit(cox_model)
  
  # Find survival at target time
  time_idx <- which.min(abs(base_surv$time - time_points_days[i]))
  S0_t <- base_surv$surv[time_idx]
  
  # Calculate predicted survival for each patient
  linear_pred <- predict(cox_model, type = "lp")
  pred_surv <- S0_t ^ exp(linear_pred)
  
  # Calculate observed survival using Kaplan-Meier
  # Divide patients into risk groups
  risk_groups <- cut(data_complete$Risk_Score, 
                     breaks = quantile(data_complete$Risk_Score, 
                                       probs = seq(0, 1, 0.2)),
                     include.lowest = TRUE,
                     labels = FALSE)
  
  obs_surv <- numeric(5)
  pred_surv_mean <- numeric(5)
  
  for(g in 1:5) {
    group_idx <- which(risk_groups == g)
    
    if(length(group_idx) > 0) {
      # Observed survival
      km_fit <- survfit(Surv(time, event) ~ 1, 
                        data = data_complete[group_idx, ])
      
      time_idx_km <- which.min(abs(km_fit$time - time_points_days[i]))
      if(length(time_idx_km) > 0) {
        obs_surv[g] <- km_fit$surv[time_idx_km]
      } else {
        obs_surv[g] <- NA
      }
      
      # Mean predicted survival
      pred_surv_mean[g] <- mean(pred_surv[group_idx])
    } else {
      obs_surv[g] <- NA
      pred_surv_mean[g] <- NA
    }
  }
  
  calibration_results[[i]] <- data.frame(
    Time_Years = time_points_years[i],
    Risk_Group = 1:5,
    Predicted_Survival = pred_surv_mean,
    Observed_Survival = obs_surv
  )
}

# Combine all calibration results
calibration_df <- do.call(rbind, calibration_results)
calibration_df <- calibration_df[complete.cases(calibration_df), ]

write.csv(calibration_df,
          "results_TNBC/LASSO_Cox/07_timeROC/calibration_data.csv",
          row.names = FALSE)

# Calibration plot
tiff("results_TNBC/LASSO_Cox/07_timeROC/calibration_plot.tiff",
     width = 12, height = 4, units = "in", res = 600, compression = "lzw")

par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

for(i in 1:length(time_points_years)) {
  
  calib_data <- calibration_df[calibration_df$Time_Years == time_points_years[i], ]
  
  if(nrow(calib_data) > 0) {
    
    plot(calib_data$Predicted_Survival, calib_data$Observed_Survival,
         pch = 19, cex = 2, col = "#E64B35",
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = "Predicted Survival Probability",
         ylab = "Observed Survival Probability",
         main = sprintf("%d-Year Calibration", time_points_years[i]),
         cex.lab = 1.2, cex.main = 1.3, cex.axis = 1.1)
    
    # Add perfect calibration line
    abline(a = 0, b = 1, lty = 2, col = "gray50", lwd = 2)
    
    # Add smoothed calibration curve
    if(nrow(calib_data) >= 3) {
      smooth_fit <- loess(Observed_Survival ~ Predicted_Survival, 
                          data = calib_data, span = 0.75)
      pred_seq <- seq(min(calib_data$Predicted_Survival), 
                      max(calib_data$Predicted_Survival), 
                      length.out = 100)
      smooth_pred <- predict(smooth_fit, newdata = pred_seq)
      lines(pred_seq, smooth_pred, col = "#4DBBD5", lwd = 2)
    }
    
    grid(col = "gray90", lty = 1)
    
    legend("bottomright",
           legend = c("Perfect calibration", "Observed"),
           col = c("gray50", "#4DBBD5"),
           lty = c(2, 1),
           lwd = 2,
           bty = "n",
           cex = 1)
  }
}

dev.off()
cat("✓ Calibration plots saved\n")

#### =====================================================
#### STEP 7: SAVE ALL MODEL OBJECTS
#### =====================================================

cat("\n========== STEP 7: SAVING MODEL OBJECTS ==========\n")

# Load final model from previous analysis or reconstruct
if(file.exists("results_TNBC/LASSO_Cox/05_final_model/final_gene_signature.csv")) {
  
  # Reconstruct final LASSO model
  X_final <- as.matrix(data[, final_genes])
  y <- Surv(time = data$time, event = data$event)
  
  set.seed(123)
  cvfit_final <- cv.glmnet(X_final, y, family = "cox", 
                           alpha = 1, nfolds = 10, standardize = TRUE)
  
  final_model <- glmnet(X_final, y, family = "cox", 
                        alpha = 1, lambda = cvfit_final$lambda.min)
  
  # Save objects
  model_objects <- list(
    final_genes = final_genes,
    final_coefs = final_genes_df,
    X_final = X_final,
    y = y,
    cvfit = cvfit_final,
    final_model = final_model,
    cox_model = cox_model,
    data = data,
    risk_scores = data$Risk_Score,
    risk_groups = data$Risk_Group
  )
  
  saveRDS(model_objects, 
          "results_TNBC/LASSO_Cox/08_model_objects/complete_model_objects.rds")
  
  cat("✓ Complete model objects saved\n")
  
  # Also save individual components
  saveRDS(final_model, 
          "results_TNBC/LASSO_Cox/08_model_objects/final_lasso_model.rds")
  saveRDS(cvfit_final, 
          "results_TNBC/LASSO_Cox/08_model_objects/cv_lasso_model.rds")
  saveRDS(cox_model, 
          "results_TNBC/LASSO_Cox/08_model_objects/cox_model.rds")
  
  cat("✓ Individual model components saved\n")
}

#### =====================================================
#### STEP 8: COMPREHENSIVE SUMMARY
#### =====================================================

cat("\n========== STEP 8: COMPREHENSIVE SUMMARY ==========\n")

# Create comprehensive summary table
time_roc_summary <- data.frame(
  Time_Point_Years = time_points_years,
  Time_Point_Days = time_points_days,
  AUC = if(!is.null(roc_obj)) auc_results$AUC else NA,
  AUC_SE = if(!is.null(roc_obj)) auc_results$SE else NA,
  AUC_Lower_95 = if(!is.null(roc_obj)) auc_results$CI_Lower else NA,
  AUC_Upper_95 = if(!is.null(roc_obj)) auc_results$CI_Upper else NA,
  AUC_survivalROC = surv_auc_summary$AUC_survivalROC,
  Events_By_Time = sapply(time_points_days, function(t) {
    sum(data_complete$event == 1 & data_complete$time <= t)
  }),
  At_Risk = sapply(time_points_days, function(t) {
    sum(data_complete$time >= t)
  })
)

write.csv(time_roc_summary,
          "results_TNBC/LASSO_Cox/07_timeROC/TIME_ROC_SUMMARY.csv",
          row.names = FALSE)

cat("\nTime-Dependent ROC Summary:\n")
cat("===========================\n")
print(time_roc_summary)
cat("\n")

# Create interpretation guide
interpretation <- data.frame(
  AUC_Range = c("0.90-1.00", "0.80-0.90", "0.70-0.80", "0.60-0.70", "0.50-0.60"),
  Interpretation = c("Excellent", "Good", "Acceptable", "Poor", "Fail")
)

cat("AUC Interpretation Guide:\n")
print(interpretation)
cat("\n")

# Assess model performance
get_auc_interpretation <- function(x){
  if(x >= 0.9) return("Excellent")
  else if(x >= 0.8) return("Good")
  else if(x >= 0.7) return("Acceptable")
  else if(x >= 0.6) return("Poor")
  else return("Fail")
}

for(i in 1:nrow(time_roc_summary)) {
  auc_val <- time_roc_summary$AUC[i]
  if(!is.na(auc_val)) {
    interp <- get_auc_interpretation(auc_val)
    cat(sprintf("%d-year AUC: %.3f [%s]\n", 
                time_roc_summary$Time_Point_Years[i],
                auc_val,
                interp))
  }
}


#### =====================================================
#### STEP 9: ADDITIONAL DIAGNOSTIC PLOTS
#### =====================================================

cat("\n========== STEP 9: ADDITIONAL DIAGNOSTIC PLOTS ==========\n")

#### Plot: Risk score distribution by event status and time
tiff("results_TNBC/LASSO_Cox/07_timeROC/risk_distribution_by_outcome.tiff",
     width = 12, height = 8, units = "in", res = 600, compression = "lzw")

par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))

# Overall distribution
hist(data_complete$Risk_Score, breaks = 30,
     main = "Overall Risk Score Distribution",
     xlab = "Risk Score",
     ylab = "Frequency",
     col = "lightblue",
     border = "white",
     cex.lab = 1.2, cex.main = 1.3)

# By event status
boxplot(Risk_Score ~ event, data = data_complete,
        names = c("Censored", "Event"),
        main = "Risk Score by Event Status",
        ylab = "Risk Score",
        col = c("lightgreen", "lightcoral"),
        cex.lab = 1.2, cex.main = 1.3)

# Density plot
plot(density(data_complete$Risk_Score[data_complete$event == 0]),
     main = "Risk Score Density by Event Status",
     xlab = "Risk Score",
     ylab = "Density",
     col = "green", lwd = 2,
     cex.lab = 1.2, cex.main = 1.3)
lines(density(data_complete$Risk_Score[data_complete$event == 1]),
      col = "red", lwd = 2)
legend("topright",
       legend = c("Censored", "Event"),
       col = c("green", "red"),
       lwd = 2,
       bty = "n")

# Risk score vs time
plot(data_complete$time_years, data_complete$Risk_Score,
     pch = ifelse(data_complete$event == 1, 19, 1),
     col = ifelse(data_complete$event == 1, "red", "blue"),
     main = "Risk Score vs Follow-up Time",
     xlab = "Time (years)",
     ylab = "Risk Score",
     cex.lab = 1.2, cex.main = 1.3)
legend("topright",
       legend = c("Event", "Censored"),
       pch = c(19, 1),
       col = c("red", "blue"),
       bty = "n")

dev.off()
cat("✓ Risk distribution plots saved\n")

#### =====================================================
#### FINAL SUMMARY
#### =====================================================

cat("\n==========================================================\n")
cat("TIME-DEPENDENT ROC ANALYSIS COMPLETE!\n")
cat("==========================================================\n\n")

cat("Key Results:\n")
cat("------------\n")

if(!is.null(roc_obj)) {
  for(i in 1:nrow(auc_results)) {
    cat(sprintf("• %d-year AUC: %.3f (95%% CI: %.3f - %.3f)\n",
                auc_results$Time_Point_Years[i],
                auc_results$AUC[i],
                auc_results$CI_Lower[i],
                auc_results$CI_Upper[i]))
  }
} else {
  cat("• timeROC analysis: Not available\n")
  cat("• survivalROC results saved separately\n")
}

cat("\nAll results saved to: results_TNBC/LASSO_Cox/07_timeROC/\n")
cat("Model objects saved to: results_TNBC/LASSO_Cox/08_model_objects/\n\n")

cat("Output files:\n")
cat("  07_timeROC/\n")
cat("    - time_dependent_ROC_curves.tiff\n")
cat("    - individual_ROC_curves.tiff\n")
cat("    - AUC_over_time.tiff\n")
cat("    - sensitivity_specificity_time.tiff\n")
cat("    - calibration_plot.tiff\n")
cat("    - risk_distribution_by_outcome.tiff\n")
cat("    - TIME_ROC_SUMMARY.csv\n")
cat("    - timeROC_AUC_summary.csv\n")
cat("    - calibration_data.csv\n\n")

cat("  08_model_objects/\n")
cat("    - complete_model_objects.rds\n")
cat("    - final_lasso_model.rds\n")
cat("    - cv_lasso_model.rds\n")
cat("    - cox_model.rds\n")
cat("    - timeROC_object.rds\n")
cat("    - survivalROC_objects.rds\n\n")

cat("Analysis completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==========================================================\n")

sink()

cat("\n✅ TIME-DEPENDENT ROC ANALYSIS COMPLETE!\n")
cat("Check results_TNBC/LASSO_Cox/07_timeROC/ for all outputs.\n\n")

#### =====================================================
#### BONUS: HOW TO LOAD AND USE SAVED MODELS
#### =====================================================

# Create a usage guide file
usage_guide <- "
==========================================================
HOW TO LOAD AND USE SAVED MODEL OBJECTS
==========================================================

# Load complete model objects
model_objects <- readRDS('results_TNBC/LASSO_Cox/08_model_objects/complete_model_objects.rds')

# Access components:
final_genes <- model_objects$final_genes
final_model <- model_objects$final_model
cox_model <- model_objects$cox_model
data <- model_objects$data

# Load individual models
lasso_model <- readRDS('results_TNBC/LASSO_Cox/08_model_objects/final_lasso_model.rds')
cv_model <- readRDS('results_TNBC/LASSO_Cox/08_model_objects/cv_lasso_model.rds')
cox_model <- readRDS('results_TNBC/LASSO_Cox/08_model_objects/cox_model.rds')

# Load ROC objects
roc_obj <- readRDS('results_TNBC/LASSO_Cox/08_model_objects/timeROC_object.rds')
surv_roc <- readRDS('results_TNBC/LASSO_Cox/08_model_objects/survivalROC_objects.rds')

# PREDICT ON NEW DATA:
# 1. Prepare new data with same genes
new_data <- ... # Your new expression data
X_new <- as.matrix(new_data[, final_genes])

# 2. Calculate risk scores
risk_scores_new <- predict(lasso_model, newx = X_new, type = 'link')

# 3. Classify into risk groups
median_risk <- median(model_objects$risk_scores)
risk_groups_new <- ifelse(risk_scores_new > median_risk, 'High', 'Low')

# 4. Predict survival probabilities
surv_prob_new <- predict(cox_model, 
                          newdata = data.frame(Risk_Score = risk_scores_new),
                          type = 'survival',
                          times = c(365.25, 365.25*3, 365.25*5))

==========================================================
"

writeLines(usage_guide, "results_TNBC/LASSO_Cox/08_model_objects/HOW_TO_USE_MODELS.txt")

cat("📖 Usage guide saved: results_TNBC/LASSO_Cox/08_model_objects/HOW_TO_USE_MODELS.txt\n")

