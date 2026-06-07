
count_matrix <- assay(brca_se, "unstranded")

# Remove Ensembl version numbers
rownames(count_matrix) <- sub("\\..*", "", rownames(count_matrix))

saveRDS(count_matrix, "data_processed/raw_counts_unstranded.rds")



