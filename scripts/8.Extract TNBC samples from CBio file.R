cbio <- read.delim(
  "data_processed/brca_tcga_clinical_data.tsv",
  stringsAsFactors = FALSE
)

colnames(cbio)

tnbc_samples <- cbio |>
  dplyr::filter(
    ER.Status.By.IHC == "Negative",
    PR.status.by.ihc == "Negative",
    IHC.HER2 == "Negative"
  ) |>
  dplyr::pull(Sample.ID) |>
  unique()

length(tnbc_samples)

saveRDS(tnbc_samples, "data_processed/tnbc_sample_ids.rds")



