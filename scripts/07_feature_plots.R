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
# Feature plots per cell type
# =========================

# Loops through gbm_markers list defined in config.R
# One JPEG saved per cell type

for (celltype in names(gbm_markers)) {

  cat("Plotting:", celltype, "\n")

  genes <- filter_genes(gbm_markers[[celltype]], obj)

  if (length(genes) == 0) next

  feature_plot <- FeaturePlot(
    obj,
    features = genes,
    order = TRUE,
    pt.size = 0.2
  )

  save_plot_jpeg(
    feature_plot,
    paste0(
      outdir,
      "/plots/FeaturePlot_",
      celltype,
      ".jpeg"
    ),
    width = 1800,
    height = 1400
  )
}

# =========================
# Key overview panel
# Tumor vs myeloid split
# =========================

overview_genes <- filter_genes(
  c("EGFR", "SOX2", "P2RY12", "CD163"),
  obj
)

overview_plot <- FeaturePlot(
  obj,
  features = overview_genes,
  ncol = 4,
  pt.size = 0.2,
  order = TRUE
)

save_plot_jpeg(
  overview_plot,
  paste0(
    outdir,
    "/plots/FeaturePlot_tumor_vs_myeloid.jpeg"
  ),
  width = 3200,
  height = 900
)

# =========================
# Microglia vs TAM panel
# Most important distinction in GBM
# =========================

mic_vs_tam_genes <- filter_genes(
  c("P2RY12", "TMEM119", "CX3CR1",
    "CD163",  "TREM2",   "SPP1"),
  obj
)

mic_vs_tam_plot <- FeaturePlot(
  obj,
  features = mic_vs_tam_genes,
  ncol = 3,
  pt.size = 0.2,
  order = TRUE
)

save_plot_jpeg(
  mic_vs_tam_plot,
  paste0(
    outdir,
    "/plots/FeaturePlot_microglia_vs_TAM.jpeg"
  ),
  width = 2400,
  height = 1600
)

# =========================
# Dot plot — all cell types
# =========================

# One representative gene per cell type
dot_genes <- filter_genes(
  c(
    "EGFR",   "MKI67", "P2RY12", "CD163",
    "CD3D",   "HAVCR2","GNLY",   "IGHG1",
    "MBP",    "PDGFRA","GFAP",   "PECAM1",
    "PDGFRB"
  ),
  obj
)

dot_plot <- DotPlot(
  obj,
  features = dot_genes
) +
  RotatedAxis()

save_plot_jpeg(
  dot_plot,
  paste0(
    outdir,
    "/plots/DotPlot_all_celltypes.jpeg"
  ),
  width = 2200,
  height = 1400
)

# =========================
# GBM state scoring
# Neftel et al. 2019 model
# =========================

cat("Scoring GBM tumour states...\n")

for (state in names(gbm_states)) {

  genes <- filter_genes(gbm_states[[state]], obj)

  if (length(genes) < 2) next

  score_col <- paste0("state_", state)

  obj <- AddModuleScore(
    obj,
    features = list(genes),
    name = score_col
  )

  # Seurat appends "1" to the name
  actual_col <- paste0(score_col, "1")

  state_plot <- FeaturePlot(
    obj,
    features = actual_col,
    pt.size = 0.2,
    order = TRUE
  ) +
    scale_colour_gradientn(
      colours = c("grey90", "#2166ac")
    ) +
    ggtitle(paste("GBM State:", state))

  save_plot_jpeg(
    state_plot,
    paste0(
      outdir,
      "/plots/GBM_state_",
      state,
      ".jpeg"
    )
  )

  cat("  Saved state plot:", state, "\n")
}

# Save object with module scores added
saveRDS(
  obj,
  file = paste0(
    outdir,
    "/objects/07_feature_plots.rds"
  )
)