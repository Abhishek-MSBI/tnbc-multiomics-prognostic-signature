query <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query, directory = "data_raw")
brca_se <- GDCprepare(query, directory = "data_raw")

saveRDS(brca_se, "data_processed/brca_se.rds")
