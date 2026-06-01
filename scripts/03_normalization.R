library(Seurat)
library(ggplot2)

source("config/config.R")
source("scripts/utils.R")

cat("Running normalization and PCA...
")

# =========================
# Load filtered object
# =========================

obj <- readRDS(
  paste0(
    outdir,
    "/objects/02_qc_filtered.rds"
  )
)

# =========================
# Normalize
# =========================

obj <- NormalizeData(obj)

# =========================
# Variable features
# =========================

obj <- FindVariableFeatures(obj)

# =========================
# Top variable genes plot
# =========================

variable_plot <- VariableFeaturePlot(obj)

save_plot_jpeg(
  variable_plot,
  paste0(
    outdir,
    "/plots/Variable_Features.jpeg"
  )
)

# =========================
# Scale data
# =========================

obj <- ScaleData(obj)

# =========================
# Save processed object
# =========================

saveRDS(
  obj,
  file = paste0(
    outdir,
    "/objects/03_normalized.rds"
  )
)