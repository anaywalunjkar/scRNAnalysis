library(Seurat)
library(openxlsx)
library(ggplot2)

source("config/config.R")
source("scripts/utils.R")

# =========================
# Create directories
# =========================

create_dirs(outdir)

cat("Loading GBM dataset: GSE162631\n")
cat("Samples:", paste(sample_ids, collapse = ", "), "\n")
cat("Mode: Tumour-only | Max", max_cells_per_sample, "cells per sample\n\n")

# =========================
# Load each sample
# =========================

# GSE162631 structure after extraction:
#   raw_counts_matrix/
#     R1_T/  barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
#     R2_T/  ...
#
# Tumour samples only (_T) — drops _N peritumoral samples
# to reduce dataset from 119k to ~20k cells for RAM efficiency.
# Subsampled to max_cells_per_sample (set in config.R) per sample.

seurat_list <- list()

for (sid in sample_ids) {

  sample_path <- file.path(input_dir, sid)

  if (!dir.exists(sample_path)) {
    stop(
      "Sample folder not found: ", sample_path,
      "\nCheck input_dir in config.R and ensure inner zips were extracted."
    )
  }

  cat("Reading:", sid, "\n")

  counts <- Read10X(
    data.dir    = sample_path,
    gene.column = 2
  )

  obj_tmp <- CreateSeuratObject(
    counts       = counts,
    project      = sid,
    min.cells    = min_cells,
    min.features = min_features
  )

  # Tag metadata
  obj_tmp$sample_id  <- sid
  obj_tmp$patient_id <- sub("_[TN]$", "", sid)
  obj_tmp$region     <- ifelse(grepl("_T$", sid), "Tumour", "Non_Tumour")

  # =========================
  # Subsample to max_cells_per_sample
  # =========================

  # Keeps biological representation intact while reducing RAM usage.
  # 119k cells crashes ScaleData on machines with less than 16GB RAM.
  # 5000 per sample = ~20k total — representative, runs on 8GB RAM.

  if (ncol(obj_tmp) > max_cells_per_sample) {
    set.seed(42)
    cells_keep <- sample(colnames(obj_tmp), max_cells_per_sample)
    obj_tmp    <- obj_tmp[, cells_keep]
    cat("  Subsampled to", max_cells_per_sample, "cells\n")
  }

  seurat_list[[sid]] <- obj_tmp

  cat("  Final cells:", ncol(obj_tmp), "| Genes:", nrow(obj_tmp), "\n")
}

# =========================
# Merge all samples
# =========================

obj <- merge(
  x            = seurat_list[[1]],
  y            = seurat_list[-1],
  add.cell.ids = sample_ids,
  project      = project_name
)

cat("\nMerged object created\n")

# =========================
# Add mitochondrial percentage
# =========================

obj[["percent.mt"]] <- PercentageFeatureSet(
  obj,
  pattern = "^MT-"
)

# =========================
# Summary
# =========================

cat("\nPer-sample cell counts:\n")
print(table(obj$sample_id))

cat("\nPer-patient cell counts:\n")
print(table(obj$patient_id))

cat("\nTumour vs Non-Tumour:\n")
print(table(obj$region))

# =========================
# Save metadata
# =========================

write.xlsx(
  obj@meta.data,
  file = paste0(outdir, "/metadata/01_raw_metadata.xlsx"),
  rowNames = TRUE
)

# =========================
# Save object
# =========================

saveRDS(
  obj,
  file = paste0(outdir, "/objects/01_raw_object.rds")
)

cat("\nRaw Seurat object saved\n")
print(obj)