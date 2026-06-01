library(Seurat)
library(openxlsx)
library(ggplot2)

source("config/config.R")
source("scripts/utils.R")

cat("Running QC analysis...
")

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
# QC violin plot
# =========================

qc_plot <- VlnPlot(
  obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

save_plot_jpeg(
  qc_plot,
  paste0(
    outdir,
    "/plots/QC_violin.jpeg"
  )
)

# =========================
# Feature scatter plots
# =========================

scatter1 <- FeatureScatter(
  obj,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
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
  feature2 = "nFeature_RNA"
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

obj <- subset(
  obj,
  subset =
    nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    percent.mt < mt_cutoff
)

cat("QC filtering complete
")

print(obj)

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
