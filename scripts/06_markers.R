library(Seurat)
library(openxlsx)

source("config/config.R")

cat("Finding marker genes...\n")

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
# Find markers
# =========================

# Using RNA assay — standard normalization was used in 03_normalization.R
# PrepSCTFindMarkers() is NOT needed here (only required with SCTransform)

markers <- FindAllMarkers(
  obj,
  assay           = "RNA",
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

cat("Total markers found:", nrow(markers), "\n")

# =========================
# Save markers
# =========================

write.xlsx(
  markers,
  file = paste0(
    outdir,
    "/markers/all_markers.xlsx"
  ),
  rowNames = FALSE
)

cat("Marker analysis complete\n")