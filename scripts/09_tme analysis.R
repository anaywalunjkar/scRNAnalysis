library(Seurat)
library(ggplot2)
library(openxlsx)
library(dplyr)

source("config/config.R")
source("scripts/utils.R")

cat("\n=========================\n")
cat("STEP 09: TME ANALYSIS\n")
cat("=========================\n")

# =========================
# Load annotated object
# =========================

obj <- readRDS(
  paste0(
    outdir,
    "/objects/08_annotated.rds"
  )
)

cat("Cells loaded:", ncol(obj), "\n")

# =========================
# Define TME compartments
# =========================

obj$compartment <- dplyr::case_when(
  obj$celltype %in% c("Tumor_GBM", "Proliferating_Tumor") ~ "Tumor",
  obj$celltype %in% c("Microglia", "TAM_Macrophage")       ~ "Myeloid",
  obj$celltype %in% c("T_Cells", "Exhausted_T_Cells",
                      "NK_Cells", "B_Plasma_Cells")        ~ "Lymphoid",
  obj$celltype %in% c("Oligodendrocytes", "OPC",
                      "Astrocytes")                         ~ "Neural",
  obj$celltype %in% c("Endothelial", "Pericytes")          ~ "Vascular",
  TRUE                                                       ~ "Other"
)

cat("\nCompartment breakdown:\n")
print(table(obj$compartment))

# =========================
# UMAP by compartment
# =========================

compartment_cols <- c(
  "Tumor"    = "#e41a1c",
  "Myeloid"  = "#4daf4a",
  "Lymphoid" = "#377eb8",
  "Neural"   = "#ff7f00",
  "Vascular" = "#984ea3",
  "Other"    = "grey70"
)

compartment_umap <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "compartment",
  cols = compartment_cols,
  pt.size = 0.3
) +
  ggtitle("GBM TME — Broad Compartments")

save_plot_jpeg(
  compartment_umap,
  paste0(
    outdir,
    "/plots/TME_compartments_UMAP.jpeg"
  )
)

# =========================
# Compartment proportions
# per patient
# =========================

comp_prop <- obj@meta.data %>%
  dplyr::count(patient_id, compartment) %>%
  dplyr::group_by(patient_id) %>%
  dplyr::mutate(proportion = n / sum(n)) %>%
  dplyr::ungroup()

comp_plot <- ggplot(
  comp_prop,
  aes(x = patient_id, y = proportion, fill = compartment)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = compartment_cols) +
  labs(
    title = "GBM — TME Compartment Proportions per Patient",
    x = "Patient",
    y = "Proportion",
    fill = "Compartment"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot_jpeg(
  comp_plot,
  paste0(
    outdir,
    "/plots/TME_compartments_by_patient.jpeg"
  )
)

# =========================
# T cell exhaustion scoring
# =========================

cat("Scoring T cell exhaustion...\n")

ex_genes <- filter_genes(exhaustion_genes, obj)

cat("Exhaustion genes found:", paste(ex_genes, collapse = ", "), "\n")

obj <- AddModuleScore(
  obj,
  features = list(ex_genes),
  name = "exhaustion_score"
)

# Seurat appends "1"
obj$exhaustion_score <- obj$exhaustion_score1

exhaustion_umap <- FeaturePlot(
  obj,
  features = "exhaustion_score",
  pt.size = 0.2,
  order = TRUE
) +
  scale_colour_gradientn(colours = c("grey90", "#d73027")) +
  ggtitle("GBM — T Cell Exhaustion Score")

save_plot_jpeg(
  exhaustion_umap,
  paste0(
    outdir,
    "/plots/TME_exhaustion_UMAP.jpeg"
  )
)

# =========================
# T cell exhaustion violin
# (T cell clusters only)
# =========================

t_obj <- subset(
  obj,
  subset = celltype %in% c("T_Cells", "Exhausted_T_Cells")
)

if (ncol(t_obj) > 10) {

  exhaustion_vln <- VlnPlot(
    t_obj,
    features = "exhaustion_score",
    group.by = "celltype",
    pt.size = 0.1
  ) +
    ggtitle("GBM — T Cell Exhaustion Score") +
    NoLegend()

  save_plot_jpeg(
    exhaustion_vln,
    paste0(
      outdir,
      "/plots/TME_exhaustion_violin.jpeg"
    )
  )
}

# =========================
# Microglia vs TAM dotplot
# =========================

myeloid_obj <- subset(
  obj,
  subset = celltype %in% c("Microglia", "TAM_Macrophage")
)

mic_tam_genes <- filter_genes(
  c("P2RY12", "TMEM119", "CX3CR1", "HEXB", "SALL1",
    "CD163",  "MRC1",    "TREM2",  "SPP1", "CCL2",
    "CD68",   "AIF1",    "CSF1R"),
  myeloid_obj
)

if (ncol(myeloid_obj) > 50 && length(mic_tam_genes) > 3) {

  mic_tam_dot <- DotPlot(
    myeloid_obj,
    features = mic_tam_genes,
    group.by = "celltype",
    assay = "RNA"
  ) +
    RotatedAxis() +
    coord_flip() +
    ggtitle("GBM — Microglia vs TAM")

  save_plot_jpeg(
    mic_tam_dot,
    paste0(
      outdir,
      "/plots/TME_microglia_vs_TAM_dotplot.jpeg"
    ),
    width = 1600,
    height = 1200
  )
}

# =========================
# Full TME dotplot
# one marker per cell type
# =========================

tme_dot_genes <- filter_genes(
  c("EGFR",   "MKI67",  "P2RY12", "CD163",
    "CD3D",   "HAVCR2", "GNLY",   "IGHG1",
    "MBP",    "PDGFRA", "GFAP",   "PECAM1",
    "PDGFRB"),
  obj
)

tme_dot_plot <- DotPlot(
  obj,
  features = tme_dot_genes,
  group.by = "celltype",
  assay = "RNA"
) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("GBM TME — Representative Markers")

save_plot_jpeg(
  tme_dot_plot,
  paste0(
    outdir,
    "/plots/TME_full_dotplot.jpeg"
  ),
  width = 2200,
  height = 1600
)

# =========================
# Save metadata + object
# =========================

write.xlsx(
  obj@meta.data,
  file = paste0(
    outdir,
    "/metadata/09_TME_metadata.xlsx"
  ),
  rowNames = TRUE
)

saveRDS(
  obj,
  file = paste0(
    outdir,
    "/objects/09_TME.rds"
  )
)

cat("\nTME analysis complete\n")