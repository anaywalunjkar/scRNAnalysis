library(Seurat)
library(ggplot2)

source("config/config.R")
source("scripts/utils.R")

cat("Running normalization...\n")

# =========================
# Load filtered object
# =========================

obj <- readRDS(
  paste0(outdir, "/objects/02_qc_filtered.rds")
)

cat("Cells loaded:", ncol(obj), "\n")

# =========================
# Join layers
# =========================

# Seurat v5 stores each sample as a separate layer after merge()
# JoinLayers() collapses them into one before normalizing —
# without this, NormalizeData only runs on the first layer

obj <- JoinLayers(obj)
cat("Layers joined\n")

# =========================
# Normalize
# =========================

# LogNormalize: divides each cell's counts by its total counts,
# multiplies by 10000 (per-cell scaling), then log1p transforms.
# This removes sequencing depth bias across cells.

obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor         = 10000,
  verbose              = FALSE
)

cat("NormalizeData complete\n")

# =========================
# Variable features
# =========================

# Identifies the 3000 genes with the most biological variability
# across cells — ignoring housekeeping genes that are uniformly
# expressed everywhere. PCA will only use these genes.

obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures        = 3000,
  verbose          = FALSE
)

cat("Variable features selected:", length(VariableFeatures(obj)), "\n")

# =========================
# Scale data
# =========================

# Centers each gene to mean=0 and scales to unit variance.
# Regresses out mitochondrial % so it doesn't drive clustering.
# Only run on variable features to save memory.

obj <- ScaleData(
  obj,
  features        = VariableFeatures(obj),
  vars.to.regress = "percent.mt",
  verbose         = FALSE
)

cat("ScaleData complete\n")

# =========================
# Variable features plot
# =========================

variable_plot <- VariableFeaturePlot(obj)

save_plot_jpeg(
  variable_plot,
  paste0(outdir, "/plots/Variable_Features.jpeg")
)

# =========================
# Save processed object
# =========================

saveRDS(
  obj,
  file = paste0(outdir, "/objects/03_normalized.rds")
)

cat("03_normalization.R complete\n")