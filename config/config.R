project_name <- "NSCLC"

# Input file
input_file <- "C:/anay/scripts/scRNA/input/20k_NSCLC_DTC_3p_nextgem_donor_1_count_sample_feature_bc_matrix.h5"

# Output directory
outdir <- "C:/anay/scripts/scRNA/output"

# QC thresholds
min_cells <- 3
min_features <- 200
max_features <- 6000
mt_cutoff <- 10

# PCA / UMAP
num_dims <- 20

# Clustering
cluster_resolution <- 0.5
