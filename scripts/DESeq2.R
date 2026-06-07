############################################################
# SUPPLEMENTARY SCRIPT S1
# DESeq2 Differential Expression Analysis
# Tumor vs Normal (TNBC)
# TCGA-BRCA Project
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

## =========================================================
## Project setup
## =========================================================
project_dir <- "C:/Users/abhir/Desktop/TNBC/TCGA_TNBC"
setwd(project_dir)

dir.create("results/DGE_DESeq2/plots",
           recursive = TRUE,
           showWarnings = FALSE)

## =========================================================
## Load processed data (GENE SYMBOL LEVEL)
## =========================================================
counts <- readRDS(
  "data_processed/gene_mapped/counts_TNBC_Normal_GeneSymbol.rds"
)

group <- readRDS(
  "data_processed/group_labels.rds"
)

stopifnot(ncol(counts) == length(group))

coldata <- data.frame(
  row.names = colnames(counts),
  condition = factor(group, levels = c("Normal", "TNBC"))
)

## =========================================================
## Construct DESeq2 dataset
## =========================================================
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = coldata,
  design    = ~ condition
)

## =========================================================
## Pre-filtering (low-count genes)
## =========================================================
dds <- dds[rowSums(counts(dds)) >= 10, ]

## =========================================================
## Run DESeq2
## =========================================================
dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "TNBC", "Normal"),
  alpha = 0.05
)

## Shrink log2FC for stability
res <- lfcShrink(
  dds,
  coef = "condition_TNBC_vs_Normal",
  res  = res,
  type = "apeglm"
)

res_df <- as.data.frame(res)
res_df$gene_symbol <- rownames(res_df)

## =========================================================
## DEG classification
## =========================================================
res_df$regulation <- "NS"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange >=  1] <- "Up"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "Down"

## =========================================================
## Save DEG tables
## =========================================================
write.csv(
  res_df,
  "results/DGE_DESeq2/DEG_DESeq2_TNBC_vs_Normal_all_genes.csv",
  row.names = FALSE
)

saveRDS(
  res_df,
  "results/DGE_DESeq2/DEG_DESeq2_TNBC_vs_Normal_all_genes.rds"
)

## =========================================================
## Volcano plot (HIGH-QUALITY, NO GENE LOSS)
## =========================================================

## Cap extreme values for visualization ONLY
res_df$log2FC_plot <- pmax(pmin(res_df$log2FoldChange, 8), -8)

res_df$negLog10Padj <- -log10(res_df$padj)
res_df$negLog10Padj[is.infinite(res_df$negLog10Padj)] <-
  max(res_df$negLog10Padj[is.finite(res_df$negLog10Padj)]) + 1

label_up <- res_df %>%
  filter(regulation == "Up") %>%
  arrange(padj) %>%
  slice_head(n = 15)

label_down <- res_df %>%
  filter(regulation == "Down") %>%
  arrange(padj) %>%
  slice_head(n = 15)


volcano <- ggplot(
  res_df,
  aes(log2FC_plot, negLog10Padj, color = regulation)
) +
  geom_point(
    data = subset(res_df, regulation == "NS"),
    size = 0.4, alpha = 0.35
  ) +
  geom_point(
    data = subset(res_df, regulation != "NS"),
    size = 1.1, alpha = 0.9
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed", linewidth = 0.4
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed", linewidth = 0.4
  ) +
  geom_text_repel(
    data = label_up,
    aes(label = gene_symbol),
    size = 3,
    color = "#d73027",
    box.padding = 0.45,
    point.padding = 0.35,
    max.overlaps = 30
  ) +
  geom_text_repel(
    data = label_down,
    aes(label = gene_symbol),
    size = 3,
    color = "#4575b4",
    box.padding = 0.45,
    point.padding = 0.35,
    max.overlaps = 30
  ) +
  scale_color_manual(values = c(
    "Up"   = "#d73027",
    "Down" = "#4575b4",
    "NS"   = "grey70"
  )) +
  theme_classic(base_size = 13) +
  labs(
    title = "Differential Expression: TNBC vs Normal",
    subtitle = "DESeq2 | |log2FC| ≥ 1 | adj. p < 0.05",
    x = "Log2 Fold Change (capped)",
    y = expression(-log[10](Adjusted~p~value))
  )

ggsave(
  "results/DGE_DESeq2/plots/Volcano_plot_DESeq2.tif",
  volcano,
  device = "tiff",
  width = 9,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

## =========================================================
## MA plot
## =========================================================
## =========================================================
## MA plot (REVIEWER-GRADE, DISTINCT FROM VOLCANO)
## =========================================================

ma_df <- as.data.frame(res)
ma_df$gene <- rownames(ma_df)

## Define DEGs conservatively
ma_df$DEG <- "NS"
ma_df$DEG[ma_df$padj < 0.05 & abs(ma_df$log2FoldChange) >= 1] <- "DEG"

## Remove zero baseMean (log scale safety)
ma_df <- ma_df[ma_df$baseMean > 0, ]

tiff(
  "results/DGE_DESeq2/plots/MA_plot_DESeq2.tif",
  width = 4200,
  height = 3000,
  res = 600,
  compression = "lzw"
)

plot(
  log10(ma_df$baseMean),
  ma_df$log2FoldChange,
  pch = 16,
  cex = 0.25,
  col = rgb(0.7, 0.7, 0.7, 0.4),
  xlab = "Log10 Mean Expression",
  ylab = "Log2 Fold Change",
  main = "MA Plot: TNBC vs Normal"
)

## Overlay DEGs only (subtle highlight)
with(
  subset(ma_df, DEG == "DEG"),
  points(
    log10(baseMean),
    log2FoldChange,
    pch = 16,
    cex = 0.35,
    col = rgb(0.8, 0.2, 0.2, 0.6)
  )
)

abline(h = c(-1, 1), lty = 2, col = "black", lwd = 1)
abline(h = 0, lwd = 1)

dev.off()

## =========================================================
## Dispersion plot
## =========================================================
tiff(
  "results/DGE_DESeq2/plots/Dispersion_plot_DESeq2.tif",
  width = 4200,
  height = 3200,
  res = 600,
  compression = "lzw"
)

plotDispEsts(
  dds,
  main = "Dispersion Estimates: TNBC vs Normal"
)

dev.off()

## =========================================================
## Heatmap: Top 50 DEGs
## =========================================================
vsd <- vst(dds, blind = FALSE)

top50 <- res_df %>%
  filter(regulation != "NS") %>%
  arrange(padj) %>%
  slice_head(n = 50)

heat_mat <- assay(vsd)[top50$gene_symbol, ]

## Order samples: Normal → TNBC
sample_order <- order(coldata$condition)
heat_mat <- heat_mat[, sample_order]

annotation_col <- data.frame(
  condition = coldata$condition[sample_order]
)
rownames(annotation_col) <- colnames(heat_mat)

tiff(
  "results/DGE_DESeq2/plots/Top50_DEG_heatmap_DESeq2.tif",
  width = 4800,
  height = 3000,
  res = 600,
  compression = "lzw"
)

pheatmap(
  heat_mat,
  scale = "row",
  annotation_col = annotation_col,
  cluster_cols = FALSE,   # 🔥 CRITICAL
  clustering_distance_rows = "correlation",
  show_colnames = FALSE,
  fontsize_row = 7,
  main = "Top 50 Differentially Expressed Genes (TNBC vs Normal)"
)

dev.off()


## =========================================================
## Reproducibility
## =========================================================
sessionInfo()
