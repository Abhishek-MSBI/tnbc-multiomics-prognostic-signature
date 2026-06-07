############################################################
## Script: 05_DEG_CEMiTool_Intersection_FINAL.R
## Project: TCGA_TNBC
##
## Purpose:
##   Identify candidate genes by intersecting:
##   - TRUE DEGs (Up + Down)
##   - TNBC-associated genes identified by CEMiTool
##
## Output:
##   Saved intersection tables for downstream analysis
##
## Reference:
##   Scientific Reports (2024)
############################################################

## ==============================
## 1. Load required libraries
## ==============================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

## ==============================
## 2. Set project root
## ==============================
setwd("C:/Users/abhir/Desktop/TNBC/TCGA_TNBC")

## ==============================
## 3. Define input paths
## ==============================

DEG_FILE <- file.path(
  "results",
  "DGE_DESeq2",
  "DEG_DESeq2_protein_coding_significant.csv"
)

CEMITOOL_OBJECT <- file.path(
  "results",
  "geneCoexpressionAnalysis",
  "CEMiTool_TNBC_Object_FINAL.rds"
)

OUTPUT_DIR <- file.path(
  "results",
  "Intersection_DEG_CEMiTool"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ==============================
## ==============================
## 4. Load DESeq2 results
## ==============================
deg_df <- read_csv(DEG_FILE, show_col_types = FALSE)

## ---- Robust gene column detection ----
possible_gene_cols <- c(
  "gene_symbol",
  "gene",
  "Gene",
  "symbol",
  "GeneSymbol",
  "ensembl_gene_id"
)

gene_col <- intersect(possible_gene_cols, colnames(deg_df))

if (length(gene_col) > 0) {
  # Use detected gene column
  deg_df <- deg_df %>%
    rename(gene = !!gene_col[1])
  message("Using gene column: ", gene_col[1])
} else if (!is.null(rownames(deg_df)) && !all(rownames(deg_df) == "")) {
  # Fall back to rownames
  deg_df <- deg_df %>%
    mutate(gene = rownames(deg_df))
  message("Using rownames as gene identifiers.")
} else {
  stop(
    "ERROR: No gene identifier column found.\n",
    "Expected one of: ",
    paste(possible_gene_cols, collapse = ", "),
    " or non-empty rownames."
  )
}

## ---- Sanity check ----
required_cols <- c("gene", "log2FoldChange", "padj", "regulation")
missing_cols <- setdiff(required_cols, colnames(deg_df))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "ERROR: DEG file missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

## ==============================
## 5. Select TRUE DEGs (paper definition)
## ==============================
deg_filtered <- deg_df %>%
  filter(regulation %in% c("Up", "Down"))

deg_genes <- unique(deg_filtered$gene)

message("True DEGs (Up + Down): ", length(deg_genes))  # 5012

## ==============================
## 6. Load CEMiTool object
## ==============================
cem <- readRDS(CEMITOOL_OBJECT)

if (is.null(cem@selected_genes)) {
  stop("ERROR: selected_genes not found in CEMiTool object.")
}

## IMPORTANT:
## In your CEMiTool object, `selected_genes` already
## represents the TNBC-associated module genes
tnbc_module_genes <- unique(cem@selected_genes)

message(
  "Genes in TNBC-associated CEMiTool module: ",
  length(tnbc_module_genes)
)  # 1447

## ==============================
## 7. DEG ∩ TNBC-module intersection
## ==============================
intersect_genes <- intersect(deg_genes, tnbc_module_genes)

message(
  "Intersected genes (DEG ∩ TNBC module): ",
  length(intersect_genes)
)  # 1030

## ==============================
## 8. Save intersection results
## ==============================

# Full annotation table
intersect_df <- deg_filtered %>%
  filter(gene %in% intersect_genes) %>%
  arrange(desc(abs(log2FoldChange)))

write_csv(
  intersect_df,
  file.path(OUTPUT_DIR, "DEG_TNBC_CEMiTool_Intersection.csv")
)

# Plain gene list (for next steps)
write_lines(
  intersect_genes,
  file.path(OUTPUT_DIR, "DEG_TNBC_CEMiTool_GeneList.txt")
)

# Summary table (Methods-ready)
summary_tbl <- tibble(
  Metric = c(
    "True DEGs (Up + Down)",
    "TNBC-module genes (CEMiTool)",
    "Intersected genes"
  ),
  Count = c(
    length(deg_genes),
    length(tnbc_module_genes),
    length(intersect_genes)
  )
)

write_csv(
  summary_tbl,
  file.path(OUTPUT_DIR, "Intersection_Summary.csv")
)

## ==============================
## 9. Save R-native objects
## ==============================

# Save intersected DEG table
saveRDS(
  intersect_df,
  file.path(OUTPUT_DIR, "DEG_TNBC_CEMiTool_Intersection.rds")
)

# Save gene vector (very useful for ML / PPI)
saveRDS(
  intersect_genes,
  file.path(OUTPUT_DIR, "DEG_TNBC_CEMiTool_GeneList.rds")
)

# Save summary table
saveRDS(
  summary_tbl,
  file.path(OUTPUT_DIR, "Intersection_Summary.rds")
)


message("==============================================")
message("DEG ∩ TNBC CEMiTool intersection COMPLETED")
message("Results saved in: ", OUTPUT_DIR)
message("==============================================")

############################################################
## END OF SCRIPT
############################################################
