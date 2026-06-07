############################################################
# TNBC Functional Enrichment & High vs Low Risk DE Analysis
# External Validation Biology (No Circular Bias)
# Date: 2026-02-17
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(123)

############################################################
# PACKAGES
############################################################
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(msigdbr)
})

############################################################
# OUTPUT DIRECTORIES
############################################################
base_out <- "functional_enrichment_results"
dir.create(base_out, showWarnings = FALSE)
dir.create(file.path(base_out, "plots"), showWarnings = FALSE)
dir.create(file.path(base_out, "tables"), showWarnings = FALSE)
dir.create(file.path(base_out, "rds_objects"), showWarnings = FALSE)

############################################################
# SECTION 1 — LOAD DATA
############################################################

# 1A — Load trained model (for signature genes only)
model_objects <- readRDS(
  "results_TNBC/LASSO_Cox/08_model_objects/complete_model_objects.rds"
)
final_genes <- model_objects$final_genes

# 1B — Load external validation risk scores
risk_results <- readRDS(
  "external_validation_results/rds_objects/05_risk_scores.rds"
)

# 1C — Load EXTERNAL raw counts (must be integers)
counts_external <- readRDS(
  "external_validation_results/rds_objects/02_counts_filtered.rds"
)

# Integer validation
stopifnot(all(counts_external == round(counts_external)))

# 1D — Load gene mapping
gene_map_objs <- readRDS(
  "external_validation_results/rds_objects/04_gene_mapping.rds"
)
mapping <- gene_map_objs$mapping

############################################################
# SECTION 2 — SIGNATURE GENE ORA (GO + KEGG ONLY)
############################################################

sig_entrez <- bitr(
  final_genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

entrez_ids <- unique(sig_entrez$ENTREZID)

# GO
go_results <- lapply(c("BP","MF","CC"), function(ont) {
  enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
})
names(go_results) <- c("BP","MF","CC")

# KEGG
kegg_res <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa",
  pvalueCutoff = 0.05
)

############################################################
# SECTION 3 — DIFFERENTIAL EXPRESSION (EXTERNAL ONLY)
############################################################

cat("\nRunning DE on external cohorts only — no circular bias\n")

# Safe alignment
common_samples <- intersect(
  colnames(counts_external),
  risk_results$Sample_ID
)

counts_de <- counts_external[, common_samples]
risk_de <- risk_results %>%
  filter(Sample_ID %in% common_samples) %>%
  arrange(match(Sample_ID, common_samples))

stopifnot(identical(colnames(counts_de), risk_de$Sample_ID))

col_data <- data.frame(
  row.names  = risk_de$Sample_ID,
  Risk_Group = factor(risk_de$Risk_Group, levels = c("Low","High")),
  Cohort     = factor(risk_de$Cohort)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_de,
  colData   = col_data,
  design    = ~ Cohort + Risk_Group
)

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("Risk_Group","High","Low")
)

res_df <- as.data.frame(res) %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(mapping, by = "ensembl_gene_id") %>%
  mutate(
    Gene_Symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "",
                         ensembl_gene_id,
                         hgnc_symbol)
  )

############################################################
# SECTION 4 — GSEA ON FULL DE RANKED LIST
############################################################

cat("\nRunning GSEA on full DE-ranked gene list\n")

hallmark_sets <- msigdbr(
  species = "Homo sapiens",
  category = "H"
) %>%
  select(gs_name, entrez_gene)

# Rank using Wald statistic (preferred)
ranked <- res_df %>%
  filter(!is.na(stat)) %>%
  select(Gene_Symbol, stat)

rank_entrez <- bitr(
  ranked$Gene_Symbol,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

ranked2 <- merge(
  ranked,
  rank_entrez,
  by.x = "Gene_Symbol",
  by.y = "SYMBOL"
)

gene_list <- ranked2$stat
names(gene_list) <- ranked2$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

gene_list <- gene_list[!duplicated(names(gene_list))]

gsea_res <- GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark_sets,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

############################################################
# SECTION 5 — SAVE RESULTS
############################################################

saveRDS(
  list(
    go_results = go_results,
    kegg_res   = kegg_res,
    dds        = dds,
    res        = res,
    res_df     = res_df,
    gsea_res   = gsea_res
  ),
  file.path(base_out,"rds_objects/COMPLETE_EXTERNAL_ANALYSIS.rds")
)

write.csv(
  res_df,
  file.path(base_out,"tables/DE_external_HighVsLow.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(gsea_res),
  file.path(base_out,"tables/GSEA_external_Hallmarks.csv"),
  row.names = FALSE
)

cat("\n========================================\n")
cat("ANALYSIS COMPLETE — EXTERNAL COHORT\n")
cat("No circular bias present\n")
cat("========================================\n")
