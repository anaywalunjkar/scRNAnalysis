library(Seurat)
library(ggplot2)

source("config/config.R")
source("scripts/utils.R")

# =========================
# Load clustered object
# =========================

obj <- readRDS(
  paste0(
    outdir,
    "/objects/05_clustered.rds"
  )
)

# =========================
# Example marker genes
# =========================

genes <- c(
  "EPCAM",
  "CD3D",
  "MS4A1"
)

# =========================
# Feature plots
# =========================

feature_plot <- FeaturePlot(
  obj,
  features = genes
)

save_plot_jpeg(
  feature_plot,
  paste0(
    outdir,
    "/plots/FeaturePlot.jpeg"
  )
)

# =========================
# Dot plot
# =========================

dot_plot <- DotPlot(
  obj,
  features = genes
) +
  RotatedAxis()

save_plot_jpeg(
  dot_plot,
  paste0(
    outdir,
    "/plots/DotPlot.jpeg"
  )
)
