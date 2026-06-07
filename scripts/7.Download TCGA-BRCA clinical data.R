clinical <- GDCquery_clinic(
  project = "TCGA-BRCA",
  type    = "clinical"
)

saveRDS(clinical, "data_processed/clinical_BRCA.rds")

grep("ER|PR|HER2|receptor|hormone", colnames(clinical), value = TRUE, ignore.case = TRUE)
