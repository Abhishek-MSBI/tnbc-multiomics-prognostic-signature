metadata <- colData(brca_se) |>
  as.data.frame() |>
  tibble::rownames_to_column("tcga_barcode") |>
  dplyr::rename(uuid = sample_id) |>
  dplyr::mutate(
    patient = substr(tcga_barcode, 1, 12)
  ) |>
  dplyr::select(
    tcga_barcode,
    uuid,
    patient,
    sample_type,
    vital_status,
    days_to_death,
    days_to_last_follow_up
  )

saveRDS(metadata, "data_processed/metadata_all_samples.rds")

head(metadata$tcga_barcode)
subtype <- TCGAbiolinks::TCGAquery_subtype(tumor = "BRCA")
