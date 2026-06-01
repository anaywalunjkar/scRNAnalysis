library(Seurat)
library(openxlsx)
library(ggplot2)

source("config/config.R")
source("scripts/utils.R")

cat("Running clustering...
")

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

obj <- FindNeighbors(
  obj,
  dims = 1:num_dims
)

# =========================
# Clustering
# =========================

obj <- FindClusters(
  obj,
  resolution = cluster_resolution
)

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

cat("Clustering complete
")
