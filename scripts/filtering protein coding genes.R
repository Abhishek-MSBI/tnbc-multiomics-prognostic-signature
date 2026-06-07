############################################################
# Filter DESeq2 DEGs to Protein-Coding Genes
# (Gene symbols already present — NO ID conversion)
############################################################

library(dplyr)
library(biomaRt)

# -------------------------------
# 1. Load DEG results
# -------------------------------
deg_results <- read.csv(
  "results/DGE_DESeq2/DEG_DESeq2_TNBC_vs_Normal_all_genes.csv",
  stringsAsFactors = FALSE
)

nrow(deg_results)
colnames(deg_results)

# -------------------------------
# 2. Get gene biotype using gene symbols
# -------------------------------
genes <- unique(deg_results$gene_symbol)

mart <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl"
)

gene_biotype <- getBM(
  attributes = c(
    "external_gene_name",
    "gene_biotype"
  ),
  filters = "external_gene_name",
  values = genes,
  mart = mart
)

# -------------------------------
# 3. Merge biotype annotation
# -------------------------------
deg_annotated <- deg_results %>%
  left_join(
    gene_biotype,
    by = c("gene_symbol" = "external_gene_name")
  )

# -------------------------------
# 4. Filter protein-coding genes
# -------------------------------
deg_protein_coding <- deg_annotated %>%
  filter(gene_biotype == "protein_coding")

# -------------------------------
# 5. DEG counts (protein-coding)
# -------------------------------
deg_counts_pc <- table(deg_protein_coding$regulation)

print(deg_counts_pc)

cat("\nProtein-coding DEG summary:\n")
cat("Upregulated      :", deg_counts_pc["Up"], "\n")
cat("Downregulated    :", deg_counts_pc["Down"], "\n")
cat("Not significant  :", deg_counts_pc["NS"], "\n")
cat("Total coding genes:", nrow(deg_protein_coding), "\n")

# -------------------------------
# 6. Save outputs
# -------------------------------
write.csv(
  deg_protein_coding,
  "results/DGE_DESeq2/DEG_DESeq2_protein_coding_all.csv",
  row.names = FALSE
)

write.csv(
  deg_protein_coding %>% filter(regulation != "Not_Significant"),
  "results/DGE_DESeq2/DEG_DESeq2_protein_coding_significant.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(deg_counts_pc),
  "results/DGE_DESeq2/DEG_DESeq2_protein_coding_counts_summary.csv"
)

