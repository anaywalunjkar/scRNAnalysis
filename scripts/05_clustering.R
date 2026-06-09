library(Seurat)
library(openxlsx)
library(ggplot2)

source("config/config.R")
source("scripts/utils.R")

cat("Running clustering...\n")

# =========================
# Load processed object
# =========================

obj <- readRDS(
  paste0(
    outdir,
    "/objects/04_dimensionality.rds"
  )
)

# =========================
# Neighbors
# =========================

# Use Harmony reduction — not PCA — because UMAP was built on Harmony
# (if you use PCA here the cluster boundaries won't match the UMAP)

obj <- FindNeighbors(
  obj,
  reduction = "harmony",
  dims = 1:num_dims
)

# =========================
# Clustering
# =========================

obj <- FindClusters(
  obj,
  resolution = cluster_resolution
)

cat("Clusters found:", length(unique(obj$seurat_clusters)), "\n")

# =========================
# Cluster UMAP
# =========================

cluster_plot <- DimPlot(
  obj,
  reduction = "umap",
  label = TRUE
) +
  theme_classic(base_size = 16)

save_plot_jpeg(
  cluster_plot,
  paste0(
    outdir,
    "/plots/UMAP_clusters.jpeg"
  )
)

# =========================
# Save metadata
# =========================

write.xlsx(
  obj@meta.data,
  file = paste0(
    outdir,
    "/metadata/cluster_metadata.xlsx"
  ),
  rowNames = TRUE
)

# =========================
# Save clustered object
# =========================

saveRDS(
  obj,
  file = paste0(
    outdir,
    "/objects/05_clustered.rds"
  )
)

cat("Clustering complete\n")