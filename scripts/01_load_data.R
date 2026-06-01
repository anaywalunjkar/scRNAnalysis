# =========================
# Load libraries
# =========================

library(Seurat)
library(openxlsx)
library(ggplot2)

# =========================
# Load config + utils
# =========================

source("config/config.R")
source("scripts/utils.R")

# =========================
# Create directories
# =========================

create_dirs(outdir)

cat("Loading 10X H5 dataset...
")

# =========================
# Read dataset
# =========================

data <- Read10X_h5(input_file)

counts <- data$`Gene Expression`

# =========================
# Create Seurat object
# =========================

obj <- CreateSeuratObject(
  counts = counts,
  project = project_name,
  min.cells = min_cells,
  min.features = min_features
)

# =========================
# Add mitochondrial percentage
# =========================

obj[["percent.mt"]] <- PercentageFeatureSet(
  obj,
  pattern = "^MT-"
)

# =========================
# Save metadata
# =========================

write.xlsx(
  obj@meta.data,
  file = paste0(
    outdir,
    "/metadata/01_raw_metadata.xlsx"
  ),
  rowNames = TRUE
)

# =========================
# Save object
# =========================

saveRDS(
  obj,
  file = paste0(
    outdir,
    "/objects/01_raw_object.rds"
  )
)

cat("Raw Seurat object created successfully
")

print(obj)
