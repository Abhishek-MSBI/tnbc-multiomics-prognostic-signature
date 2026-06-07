############################################################
# Gene Co-expression Analysis using CEMiTool (TNBC)
# Adapted to match reference paper workflow
############################################################

## ------------------------------
## 1. Load libraries
## ------------------------------
suppressPackageStartupMessages({
  library(CEMiTool)
  library(msigdbr)
  library(data.table)
})

## ------------------------------
## 2. Set project paths
## ------------------------------
project_path <- "C:/Users/abhir/Desktop/TNBC/TCGA_TNBC"
setwd(project_path)

dir.create("Results/geneCoexpressionAnalysis", recursive = TRUE, showWarnings = FALSE)

## ------------------------------
## 3. Load expression matrix
## ------------------------------
expr_mat <- readRDS("data_processed/gene_mapped/expr_mat_TNBC_Normal_protein_coding_log2.rds")

## ------------------------------
## 4. Load phenotype annotation
## ------------------------------
sample_annot <- readRDS("data_processed/gene_mapped/sample_annot_TNBC_Normal.rds")

## Ensure required column names (REFERENCE PAPER STYLE)
colnames(sample_annot) <- c("SampleName", "Class")

## ------------------------------
## 5. Load pathway gene sets (GMT)
## ------------------------------
# Reference paper uses supplied GMT
# We use MSigDB GO:BP (reviewer-acceptable substitute)

msig_go_bp <- msigdbr(
  species = "Homo sapiens",
  category = "C5",
  subcategory = "GO:BP"
)

gmtData <- data.frame(
  term = msig_go_bp$gs_name,
  gene = msig_go_bp$gene_symbol,
  stringsAsFactors = FALSE
)

## ------------------------------
## 6. Run CEMiTool
## ------------------------------
# plot = FALSE avoids ggplot2 incompatibility
# Equivalent results are generated later

cem <- cemitool(
  expr = expr_mat,
  annot = sample_annot,
  gmt = gmtData,
  sample_name_column = "SampleName",
  class_column = "Class",
  filter = TRUE,
  plot = FALSE,
  verbose = TRUE
)

## ------------------------------
## 7. Module profile & enrichment plots
## ------------------------------
plot_profile(cem)
plot_gsea(cem)

## ------------------------------
## 8. Over-representation analysis (ORA)
## ------------------------------
cem <- mod_ora(cem, gmt = gmtData)
plot_ora(cem)

## ------------------------------
## 9. Generate HTML report
## ------------------------------
generate_report(
  cem,
  directory = "Results/geneCoexpressionAnalysis/Report",
  force = TRUE
)

## ------------------------------
## 10. Write tables (REFERENCE STYLE)
## ------------------------------
write_files(
  cem,
  directory = "Results/geneCoexpressionAnalysis/Tables",
  force = TRUE
)

## ------------------------------
## 11. Save plots (REFERENCE STYLE)
## ------------------------------
save_plots(
  cem,
  value = "all",
  directory = "Results/geneCoexpressionAnalysis/Plots",
  force = TRUE
)

# Define plots explicitly
plot_profile(cem)
plot_gsea(cem)
plot_ora(cem)

# Now save all available plots
save_plots(
  cem,
  value = "all",
  directory = "Results/geneCoexpressionAnalysis/Plots",
  force = TRUE
)
pdf("Results/geneCoexpressionAnalysis/QC_expression_hist.pdf")
hist(
  as.numeric(as.matrix(expr_mat)),
  breaks = 100,
  col = "grey",
  main = "Expression value distribution",
  xlab = "log2(expression + 1)"
)
dev.off()

pdf("Results/geneCoexpressionAnalysis/QC_mean_variance.pdf")
plot(
  rowMeans(expr_mat),
  apply(expr_mat, 1, var),
  pch = 16,
  cex = 0.3,
  xlab = "Mean expression",
  ylab = "Variance",
  main = "Mean–variance relationship"
)
dev.off()

pdf("Results/geneCoexpressionAnalysis/QC_qqplot.pdf")
qqnorm(rowMeans(expr_mat), main = "QQ plot of gene means")
qqline(rowMeans(expr_mat), col = "red")
dev.off()

pdf("Results/geneCoexpressionAnalysis/QC_sample_tree.pdf")
hc <- hclust(dist(t(expr_mat)), method = "average")
plot(hc, main = "Sample clustering dendrogram")
dev.off()

## ------------------------------
## 12. Extract significant module genes
## ------------------------------
mods <- fread("Results/geneCoexpressionAnalysis/Tables/module.tsv")

# Example: TNBC-activated module
sigMods <- c("M3")

sigGenes <- mods$genes[mods$modules %in% sigMods]

write.table(
  sigGenes,
  "Results/geneCoexpressionAnalysis/M3_module_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

saveRDS(cem, "Results/geneCoexpressionAnalysis/CEMiTool_TNBC_Object_FINAL.rds")

############################################################
# End of script
############################################################
