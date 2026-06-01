library(Seurat)
library(SingleR)
library(celldex)

source("config/config.R")

cat("\n=========================\n")
cat("STEP 08: CELL ANNOTATION\n")
cat("=========================\n")

# Load clustered object

obj <- readRDS(
file.path(outdir, "/objects/05_clustered.rds")
)

cat("Cells loaded:", ncol(obj), "\n")

# Load reference

cat("Loading reference dataset...\n")

ref <- HumanPrimaryCellAtlasData()

cat("Reference loaded\n")

# Run SingleR

cat("Running SingleR annotation...\n")

pred <- SingleR(
test = GetAssayData(
obj,
assay = "RNA",
layer = "data"
),
ref = ref,
labels = ref$label.main
)

cat("Annotation complete\n")

# Add labels

obj$SingleR.labels <- pred$labels

# Print counts

cat("\nCell Type Counts:\n")

print(
table(obj$SingleR.labels)
)

# Save label table

write.csv(
as.data.frame(
table(obj$SingleR.labels)
),
file.path(
outdir,
"SingleR_celltype_counts.csv"
),
row.names = FALSE
)

# UMAP

p <- DimPlot(
obj,
group.by = "SingleR.labels",
label = TRUE,
repel = TRUE
)

jpeg(
annotation_plot,
width = 2200,
height = 1800,
res = 300
)

print(p)

dev.off()

# Save metadata

write.csv(
  obj@meta.data,
  annotation_metadata,
  row.names = TRUE
)

# Save object

saveRDS(
obj,
annotated_rds
)

cat("\nAnnotation outputs saved\n")

