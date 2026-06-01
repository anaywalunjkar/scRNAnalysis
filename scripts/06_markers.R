library(Seurat)
library(openxlsx)

source("config/config.R")

cat("Finding marker genes...
")

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

markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

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

cat("Marker analysis complete
")
