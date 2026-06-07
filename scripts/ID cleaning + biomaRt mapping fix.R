## Strip Ensembl version numbers (CRITICAL)
res_df$ensembl_id_clean <- sub("\\..*", "", res_df$ensembl_id)
mart <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  mirror  = "useast"
)
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = unique(res_df$ensembl_id_clean),
  mart       = mart
)
str(gene_map)

res_df <- res_df %>%
  left_join(
    gene_map,
    by = c("ensembl_id_clean" = "ensembl_gene_id")
  )

res_df$label <- ifelse(
  is.na(res_df$hgnc_symbol) | res_df$hgnc_symbol == "",
  res_df$ensembl_id_clean,
  res_df$hgnc_symbol
)
