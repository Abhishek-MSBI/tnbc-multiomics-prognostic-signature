############################################################
# SUPPLEMENTARY SCRIPT S2
# Functional Enrichment Analysis (GO & KEGG)
# TNBC vs Normal | TCGA
# Reference-matched methodology (Scientific Reports 2024)
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(tidyverse)
  library(readr)
})

############################################################
# 1. Set WORKING DIRECTORY (IMPORTANT)
############################################################

setwd("C:/Users/abhir/Desktop/TNBC/TCGA_TNBC")

project_dir <- getwd()

############################################################
# 2. Define directories
############################################################

dge_dir    <- file.path(project_dir, "results", "DGE_DESeq2")
enrich_dir <- file.path(project_dir, "results", "enrichment")
fig_dir    <- file.path(project_dir, "results", "figures")

dir.create(file.path(enrich_dir, "GO"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(enrich_dir, "KEGG"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "GO"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "KEGG"), recursive = TRUE, showWarnings = FALSE)

############################################################
# 3. Load DEG data
############################################################

deg_file <- file.path(
  dge_dir,
  "DEG_DESeq2_protein_coding_significant.csv"
)

deg <- read_csv(deg_file)

# REQUIRED column check
stopifnot("gene_symbol" %in% colnames(deg))

gene_symbols <- unique(deg$gene_symbol)

############################################################
# 4. Convert Gene Symbols → Entrez IDs
############################################################

gene_entrez <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

entrez_ids <- unique(gene_entrez$ENTREZID)

## GO enrichment
ontologies <- c("BP")

for (ont in ontologies) {
  
 
  go_raw <- enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    ont           = ont,
    keyType       = "ENTREZID",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  ## Convert to data frame
  go_df <- as.data.frame(go_raw)
  
  ## Remove empty results safely
  if (nrow(go_df) == 0) next
  
  ## Rank by statistical significance
  go_df <- go_df |>
    arrange(desc(Count))
  
  ## Remove duplicate / highly similar descriptions
  go_df <- go_df |>
    distinct(Description, .keep_all = TRUE)
  
  ## Keep top terms only (manuscript standard)
  go_df <- go_df |>
    slice_head(n = 10)
  
  ## Wrap long GO names
  go_df$Description <- str_wrap(go_df$Description, width = 45)
  
  ## Save clean table
  write.csv(
    go_df,
    file = file.path(
      enrich_dir,
      "GO",
      paste0("GO_", ont, "_enrichment_results_clean.csv")
    ),
    row.names = FALSE
  )
  
  ## Reconstruct enrichResult for plotting
  go_plot <- go_raw
  go_plot@result <- go_plot@result[
    go_plot@result$ID %in% go_df$ID, ]
  
  ## Order GO terms by gene count (same behavior as KEGG)
  go_plot@result <- go_plot@result |>
    arrange(Count)
  
  ## Fix y-axis order explicitly
  go_plot@result$Description <- factor(
    go_plot@result$Description,
    levels = go_plot@result$Description
  )
  
  ## High-resolution TIFF
  tiff(
    filename = file.path(
      fig_dir,
      "GO",
      paste0("GO_", ont, "_dotplot_publication.tiff")
    ),
    width = 8,
    height = 6,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  
  p <- dotplot(
    go_plot,
    showCategory = 10,
    orderBy = "Count",
    color = "p.adjust"
  ) +
    scale_color_continuous(
      low  = "#B2182B",
      high = "#2166AC",
      name = "Adjusted P-value",
      guide = guide_colorbar(reverse = TRUE)
    ) +
    labs(
      title = paste("GO", ont, "Enrichment Analysis"),
      subtitle = "TNBC vs Normal",
      x = "Gene Ratio",
      y = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title       = element_text(face = "bold", size = 16),
      plot.subtitle    = element_text(size = 12),
      axis.text.y      = element_text(size = 12),
      axis.text.x      = element_text(size = 11),
      legend.title     = element_text(size = 12),
      legend.text      = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  dev.off()
  
############################################################
# 6. KEGG Pathway Enrichment
############################################################

kegg_enrich <- enrichKEGG(
  gene          = entrez_ids,
  organism      = "hsa",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

kegg_enrich <- setReadable(
  kegg_enrich,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

write.csv(
  as.data.frame(kegg_enrich),
  file = file.path(enrich_dir, "KEGG",
                   "KEGG_enrichment_results.csv"),
  row.names = FALSE
)

tiff(
  filename = file.path(fig_dir, "KEGG", "KEGG_dotplot.tiff"),
  width = 7,
  height = 5,
  units = "in",
  res = 600,
  compression = "lzw"
)

print(
  dotplot(kegg_enrich,
          showCategory = 15,
          font.size = 12,
          title = "KEGG Pathway Enrichment – TNBC vs Normal")
)

dev.off()

############################################################
# 7. Save session info
############################################################

writeLines(
  capture.output(sessionInfo()),
  con = file.path(project_dir,
                  "results",
                  "GO_KEGG_sessionInfo.txt")
)

############################################################
# END OF SCRIPT
############################################################
