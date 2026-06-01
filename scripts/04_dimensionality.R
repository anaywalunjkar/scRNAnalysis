library(Seurat)
library(ggplot2)

source("config/config.R")
source("scripts/utils.R")

cat("Running dimensionality reduction...\n")

# =========================
# Load normalized object
# =========================

obj <- readRDS(
  paste0(
    outdir,
    "/objects/03_normalized.rds"
  )
)

# =========================
# PCA
# =========================

obj <- RunPCA(obj)

# =========================
# PCA plot
# =========================

pca_plot <- DimPlot(
  obj,
  reduction = "pca"
)

save_plot_jpeg(
  pca_plot,
  paste0(
    outdir,
    "/plots/PCA.jpeg"
  )
)

# =========================
# Elbow plot
# =========================

elbow_plot <- ElbowPlot(obj)

save_plot_jpeg(
  elbow_plot,
  paste0(
    outdir,
    "/plots/ElbowPlot.jpeg"
  )
)

# =========================
# UMAP
# =========================

obj <- RunUMAP(
  obj,
  dims = 1:num_dims
)

# =========================
# UMAP plot
# =========================

umap_plot <- DimPlot(
  obj,
  reduction = "umap"
)

save_plot_jpeg(
  umap_plot,
  paste0(
    outdir,
    "/plots/UMAP.jpeg"
  )
)

# =========================
# Save dimensionality object
# =========================

saveRDS(
  obj,
  file = paste0(
    outdir,
    "/objects/04_dimensionality.rds"
  )
)

cat("Dimensionality reduction complete\n")