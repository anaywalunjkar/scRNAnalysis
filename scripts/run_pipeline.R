setwd("C:/anay/scripts/scRNA")

start_time <- Sys.time()
cat("                 /| _ /|                      \n")
cat("===============  ( o.o ) ====================\n")
cat("                  > ^ <                      \n")
cat("===================================\n")
cat("Starting scRNA-seq pipeline by anay\n")
cat("===================================\n")

source("scripts/01_load_data.R")
source("scripts/02_qc.R")
source("scripts/03_normalization.R")
source("scripts/04_dimensionality.R")
source("scripts/05_clustering.R")
source("scripts/06_markers.R")
source("scripts/07_feature_plots.R")

cat("===================================\n")
cat("Pipeline completed successfully\n")
cat("===================================\n")

end_time <- Sys.time()

runtime <- end_time - start_time

cat("\n")
cat("Pipeline runtime:\n")
print(runtime)
cat("\n")

cat("\n")
cat("===================================\n")
cat("scRNA-seq PIPELINE COMPLETED SUCCESSFULLY\n")
cat("===================================\n")
cat("\n")

cat("Outputs generated:\n")
cat("- QC plots\n")
cat("- PCA plots\n")
cat("- UMAP plots\n")
cat("- Cluster analysis\n")
cat("- Marker genes\n")
cat("- Metadata exports\n")
cat("- Seurat objects (.rds)\n")
cat("\n")

cat("Results saved in:\n")
cat("C:/anay/scripts/scRNA/output\n")
cat("\n")