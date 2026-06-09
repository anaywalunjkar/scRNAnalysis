library(Seurat)
library(ggplot2)
library(harmony)

source("config/config.R")
source("scripts/utils.R")

cat("Running dimensionality reduction...\n")

# =========================
# Load normalized object
# =========================

obj <- readRDS(
  paste0(outdir, "/objects/03_normalized.rds")
)

# =========================
# PCA
# =========================

obj <- RunPCA(
  obj,
  npcs    = 50,
  verbose = FALSE
)

# =========================
# PCA plot — coloured by patient
# (shows batch effect before Harmony)
# =========================

pca_plot <- DimPlot(
  obj,
  reduction = "pca",
  group.by  = "patient_id"
)

save_plot_jpeg(
  pca_plot,
  paste0(outdir, "/plots/PCA_by_patient_pre_harmony.jpeg")
)

# =========================
# Elbow plot
# =========================

elbow_plot <- ElbowPlot(obj, ndims = 50)

save_plot_jpeg(
  elbow_plot,
  paste0(outdir, "/plots/ElbowPlot.jpeg")
)

# =========================
# Harmony batch correction
# =========================

# Required for GSE162631 (4 patients).
# Corrects patient-level batch effects in PCA space.
# UMAP and clustering use the Harmony embedding — NOT raw PCA.
#
# NOTE: Harmony2 v2.0+ changed argument names:
#   reduction     -> reduction.use
#   assay.use     -> removed (deprecated)
# Using reduction.use to match Harmony2 API.

cat("Running Harmony batch correction on patient_id...\n")

obj <- RunHarmony(
  obj,
  group.by.vars = batch_variable,
  reduction.use = "pca",
  dims.use      = 1:num_dims,
  verbose       = FALSE
)

cat("Harmony complete\n")

# =========================
# Harmony plot
# (patients should now be mixed)
# =========================

harmony_plot <- DimPlot(
  obj,
  reduction = "harmony",
  group.by  = "patient_id"
)

save_plot_jpeg(
  harmony_plot,
  paste0(outdir, "/plots/Harmony_by_patient.jpeg")
)

# =========================
# UMAP on Harmony embedding
# =========================

obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims      = 1:num_dims
)

# =========================
# UMAP plots
# =========================

umap_plot <- DimPlot(
  obj,
  reduction = "umap"
)

save_plot_jpeg(
  umap_plot,
  paste0(outdir, "/plots/UMAP.jpeg")
)

umap_by_patient <- DimPlot(
  obj,
  reduction = "umap",
  group.by  = "patient_id"
)

save_plot_jpeg(
  umap_by_patient,
  paste0(outdir, "/plots/UMAP_by_patient.jpeg")
)

# =========================
# Save object
# =========================

saveRDS(
  obj,
  file = paste0(outdir, "/objects/04_dimensionality.rds")
)

cat("Dimensionality reduction complete\n")