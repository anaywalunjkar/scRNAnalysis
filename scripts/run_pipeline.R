setwd("C:/anay/scripts/scRNA/GBM")

start_time <- Sys.time()

cat("            .    .                 \n")
cat("           /| _ /|                 \n")
cat("          (  o.o  )                \n")
cat("            > ^ <                  \n")
cat("===================================\n")
cat("  GBM scRNA-seq pipeline by anay   \n")
cat("  Dataset: GSE162631               \n")
cat("===================================\n")

# =========================
# Phase 1: Core pipeline
# Run all steps to clustering
# =========================

source("scripts/01_load_data.R")
source("scripts/02_qc.R")
source("scripts/03_normalization.R")
source("scripts/04_dimensionality.R")
source("scripts/05_clustering.R")
source("scripts/06_markers.R")
source("scripts/07_feature_plots.R")

cat("===================================\n")
cat("Phase 1 complete\n")
cat("NEXT: Open output/markers/all_markers.xlsx\n")
cat("      Review output/plots/FeaturePlot_*.jpeg\n")
cat("      Edit cluster_labels in scripts/08_annotation.R\n")
cat("===================================\n")

# =========================
# Phase 2: Annotation + TME
# Run after editing 08_annotation.R
# =========================

source("scripts/08_annotation.R")
source("scripts/09_TME_analysis.R")

cat("===================================\n")
cat("Pipeline completed successfully\n")
cat("===================================\n")

end_time <- Sys.time()
runtime <- end_time - start_time

cat("\n")
cat("Pipeline runtime:\n")
print(runtime)
cat("\n")

cat("===================================\n")
cat("scRNA-seq PIPELINE COMPLETED SUCCESSFULLY\n")
cat("===================================\n")
cat("\n")

cat("Outputs generated:\n")
cat("- QC plots (per patient)\n")
cat("- Harmony batch correction plots\n")
cat("- PCA + UMAP plots\n")
cat("- Cluster analysis\n")
cat("- Marker genes\n")
cat("- Feature plots (per cell type)\n")
cat("- GBM state scores (MES/AC/OPC/NPC)\n")
cat("- Annotated UMAP\n")
cat("- TME compartment analysis\n")
cat("- T cell exhaustion scores\n")
cat("- Microglia vs TAM dotplot\n")
cat("- Metadata exports (.xlsx)\n")
cat("- Seurat objects (.rds)\n")
cat("\n")

cat("Results saved in:\n")
cat("C:/anay/scripts/scRNA_GBM/output\n")
cat("\n")