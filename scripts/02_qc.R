library(Seurat)
library(openxlsx)
library(ggplot2)

source("config/config.R")
source("scripts/utils.R")

cat("Running QC analysis...\n")

# =========================
# Load object
# =========================

obj <- readRDS(
  paste0(
    outdir,
    "/objects/01_raw_object.rds"
  )
)

# =========================
# QC violin plot — by patient
# =========================

# Group by patient_id so per-sample outliers are visible
# (single-sample NSCLC used default grouping)

qc_plot <- VlnPlot(
  obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "patient_id",
  ncol = 3,
  pt.size = 0
)

save_plot_jpeg(
  qc_plot,
  paste0(
    outdir,
    "/plots/QC_violin_by_patient.jpeg"
  ),
  width = 2400,
  height = 1200
)

# =========================
# Feature scatter plots
# =========================

scatter1 <- FeatureScatter(
  obj,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  group.by = "patient_id"
)

save_plot_jpeg(
  scatter1,
  paste0(
    outdir,
    "/plots/Scatter_nCount_vs_MT.jpeg"
  )
)

scatter2 <- FeatureScatter(
  obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "patient_id"
)

save_plot_jpeg(
  scatter2,
  paste0(
    outdir,
    "/plots/Scatter_nCount_vs_nFeature.jpeg"
  )
)

# =========================
# QC filtering
# =========================

# GBM thresholds (see config.R):
#   max_features = 8000  (NSCLC was 6000 — glia express more genes)
#   mt_cutoff    = 20    (NSCLC was 10  — CNS cells have higher baseline mt)

cells_before <- ncol(obj)

obj <- subset(
  obj,
  subset =
    nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    percent.mt < mt_cutoff
)

cells_after <- ncol(obj)

cat("QC filtering complete\n")
cat("Cells before:", cells_before, "\n")
cat("Cells after: ", cells_after, "\n")
cat("Removed:     ", cells_before - cells_after, "\n")

print(obj)

cat("Per-patient counts after QC:\n")
print(table(obj$patient_id))

# =========================
# Save metadata
# =========================

write.xlsx(
  obj@meta.data,
  file = paste0(
    outdir,
    "/metadata/02_filtered_metadata.xlsx"
  ),
  rowNames = TRUE
)

# =========================
# Save filtered object
# =========================

saveRDS(
  obj,
  file = paste0(
    outdir,
    "/objects/02_qc_filtered.rds"
  )
)