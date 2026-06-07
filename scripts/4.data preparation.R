## ================================
## TCGA-BRCA: Normal tissue extraction
## ================================

rm(list = ls())
options(stringsAsFactors = FALSE)

library(data.table)

## ---- Project root ----
project_dir <- "C:/Users/abhir/Desktop/TNBC/TCGA_TNBC"
setwd(project_dir)

## ---- Paths ----
raw_dir <- file.path(
  "data_raw",
  "TCGA-BRCA",
  "Transcriptome_Profiling",
  "Gene_Expression_Quantification"
)

out_dir <- "data_processed"

dir.create(out_dir, showWarnings = FALSE)

## ---- List all files ----
files <- list.files(
  raw_dir,
  recursive = TRUE,
  full.names = TRUE
)

stopifnot(length(files) > 0)

## ---- Extract TCGA barcodes ----
get_barcode <- function(x) {
  sub(".*(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}).*", "\\1", x)
}

barcodes <- sapply(files, get_barcode)

## ---- Sample type ----
sample_type <- substr(barcodes, 14, 15)

## ---- Solid tissue normal (11) ----
normal_idx <- which(sample_type == "11")

normal_files <- files[normal_idx]
normal_barcodes <- barcodes[normal_idx]

cat("Normal samples identified:", length(normal_files), "\n")

stopifnot(length(normal_files) > 0)

## ---- Read STAR counts ----
read_star <- function(f) {
  fread(f, skip = 4)[, .(gene_id = V1, count = V2)]
}

count_list <- lapply(normal_files, read_star)

genes <- count_list[[1]]$gene_id

normal_counts <- do.call(
  cbind,
  lapply(count_list, function(x) x$count)
)

rownames(normal_counts) <- genes
colnames(normal_counts) <- normal_barcodes

## ---- Save ----
saveRDS(
  normal_counts,
  file.path(out_dir, "counts_normal_raw.rds")
)

cat("✅ Normal count matrix saved successfully\n")
