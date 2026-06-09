library(Seurat)
library(openxlsx)
library(ggplot2)
library(dplyr)

source("config/config.R")
source("scripts/utils.R")

cat("\n=========================\n")
cat("STEP 08: CELL ANNOTATION\n")
cat("=========================\n")

# =========================
# Load object
# =========================

obj <- readRDS(
  paste0(outdir, "/objects/07_feature_plots.rds")
)

cat("Cells loaded:", ncol(obj), "\n")
cat("Clusters found:", length(unique(obj$seurat_clusters)), "\n\n")

cat("Cluster sizes:\n")
print(table(obj$seurat_clusters))

# =========================
# Manual cluster labels
# =========================

# !! EDIT THESE after reviewing your markers and feature plots !!
# All 16 clusters (0-15) must have a label or the script crashes.
#
# How to annotate:
#   1. Open output/markers/all_markers.xlsx
#   2. Filter by cluster, sort by avg_log2FC descending
#   3. Top 5-10 genes = the identity of that cluster
#   4. Cross-check with output/plots/FeaturePlot_*.jpeg
#
# Key GBM rules:
#   P2RY12+ TMEM119+ CX3CR1+  → Microglia (resident brain immune)
#   CD163+  TREM2+   SPP1+    → TAM_Macrophage (infiltrating)
#   MKI67+  TOP2A+            → Proliferating (any cell type)
#   MBP+    MOG+     PLP1+    → Oligodendrocytes
#   GFAP+   AQP4+   S100B+   → Astrocytes
#   CD3D+   CD3E+            → T_Cells
#   PDCD1+  HAVCR2+ TIGIT+   → Exhausted_T_Cells

cluster_labels <- c(
  "0"  = "Tumor_GBM",           # placeholder — edit after reviewing markers
  "1"  = "Tumor_GBM",           # placeholder
  "2"  = "Tumor_GBM",           # placeholder
  "3"  = "Microglia",           # placeholder
  "4"  = "TAM_Macrophage",      # placeholder
  "5"  = "Tumor_GBM",           # placeholder
  "6"  = "Oligodendrocytes",    # placeholder
  "7"  = "T_Cells",             # placeholder
  "8"  = "Astrocytes",          # placeholder
  "9"  = "Proliferating_Tumor", # placeholder
  "10" = "Exhausted_T_Cells",   # placeholder
  "11" = "OPC",                 # placeholder
  "12" = "NK_Cells",            # placeholder
  "13" = "B_Plasma_Cells",      # placeholder
  "14" = "Endothelial",         # placeholder
  "15" = "Pericytes"            # placeholder
)

# =========================
# Apply labels
# =========================

# Using AddMetaData with named vector — more robust than direct assignment
# avoids the "No cell overlap" error from direct $ assignment

celltype_vec        <- cluster_labels[as.character(obj$seurat_clusters)]
names(celltype_vec) <- colnames(obj)

obj <- AddMetaData(
  obj,
  metadata = celltype_vec,
  col.name = "celltype"
)

# Check for any unlabelled clusters
unlabelled <- unique(obj$seurat_clusters[is.na(obj$celltype)])
if (length(unlabelled) > 0) {
  cat("WARNING: unlabelled clusters:", paste(unlabelled, collapse = ", "), "\n")
  cat("Add them to cluster_labels in 08_annotation.R\n")
}

cat("\nAnnotation summary:\n")
print(table(obj$celltype))

# =========================
# Colour palette
# =========================

gbm_colours <- c(
  "Tumor_GBM"           = "#e41a1c",
  "Proliferating_Tumor" = "#ff7f00",
  "Microglia"           = "#4daf4a",
  "TAM_Macrophage"      = "#984ea3",
  "T_Cells"             = "#377eb8",
  "Exhausted_T_Cells"   = "#1f78b4",
  "NK_Cells"            = "#a6cee3",
  "B_Plasma_Cells"      = "#fb9a99",
  "Oligodendrocytes"    = "#33a02c",
  "OPC"                 = "#b2df8a",
  "Astrocytes"          = "#fdbf6f",
  "Endothelial"         = "#cab2d6",
  "Pericytes"           = "#6a3d9a"
)

# =========================
# Annotated UMAP
# =========================

annotation_plot <- DimPlot(
  obj,
  reduction = "umap",
  group.by  = "celltype",
  label     = TRUE,
  repel     = TRUE,
  cols      = gbm_colours
) +
  ggtitle("GBM — Annotated Cell Types")

save_plot_jpeg(
  annotation_plot,
  paste0(outdir, "/plots/UMAP_annotated.jpeg"),
  width  = 2200,
  height = 1800
)

# =========================
# UMAP split by patient
# =========================

split_plot <- DimPlot(
  obj,
  reduction = "umap",
  group.by  = "celltype",
  split.by  = "patient_id",
  cols      = gbm_colours,
  pt.size   = 0.2
)

save_plot_jpeg(
  split_plot,
  paste0(outdir, "/plots/UMAP_annotated_by_patient.jpeg"),
  width  = 3600,
  height = 1000
)

# =========================
# Cell type proportions
# per patient bar chart
# =========================

prop_data <- obj@meta.data %>%
  dplyr::count(patient_id, celltype) %>%
  dplyr::group_by(patient_id) %>%
  dplyr::mutate(proportion = n / sum(n)) %>%
  dplyr::ungroup()

prop_plot <- ggplot(
  prop_data,
  aes(x = patient_id, y = proportion, fill = celltype)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = gbm_colours, na.value = "grey70") +
  labs(
    title = "GBM — Cell Type Proportions per Patient",
    x     = "Patient",
    y     = "Proportion",
    fill  = "Cell Type"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot_jpeg(
  prop_plot,
  paste0(outdir, "/plots/Celltype_proportions_by_patient.jpeg"),
  width  = 1800,
  height = 1400
)

# =========================
# Save metadata + object
# =========================

write.xlsx(
  obj@meta.data,
  file     = paste0(outdir, "/metadata/08_annotated_metadata.xlsx"),
  rowNames = TRUE
)

saveRDS(
  obj,
  file = paste0(outdir, "/objects/08_annotated.rds")
)

cat("\nAnnotation outputs saved\n")
cat("NEXT: Review plots and markers, then update cluster_labels above\n")