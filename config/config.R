# =========================
# GBM Pipeline Config
# Dataset: GSE162631
# =========================

project_name <- "GBM"

# =========================
# Input
# =========================

input_dir <- "C:/anay/scripts/scRNA/GBM/input/raw_counts_matrix"

# Tumour samples only — drops _N peritumoral to reduce dataset size.
# Add "R1_N","R2_N","R3_N","R4_N" back here once you have more RAM
# or are running on a server.

sample_ids <- c(
  "R1_T", "R2_T", "R3_T", "R4_T"
)

# =========================
# Output
# =========================

outdir <- "C:/anay/scripts/scRNA/GBM/output"

# =========================
# Subsampling
# =========================

# Max cells to keep per sample after loading.
# 5000 x 4 tumour samples = ~20000 cells total.
# Runs comfortably on 8GB RAM.
# Set to Inf to disable subsampling (requires 16GB+ RAM).

max_cells_per_sample <- 5000

# =========================
# QC thresholds
# =========================

min_cells    <- 3
min_features <- 200
max_features <- 8000
mt_cutoff    <- 20

# =========================
# Dimensionality + clustering
# =========================

num_dims           <- 30
cluster_resolution <- 0.6

# =========================
# Batch correction
# =========================

batch_variable <- "patient_id"

# =========================
# GBM marker genes
# =========================

gbm_markers <- list(

  Tumor_GBM = c(
    "EGFR", "PTEN", "CDK4", "PDGFRA", "SOX2", "PTPRZ1"
  ),

  Proliferating_Tumor = c(
    "MKI67", "TOP2A", "PCNA", "UBE2C"
  ),

  Microglia = c(
    "P2RY12", "TMEM119", "CX3CR1", "HEXB", "SALL1"
  ),

  TAM_Macrophage = c(
    "CD163", "MRC1", "TREM2", "SPP1", "CCL2", "FCGR1A"
  ),

  T_Cells = c(
    "CD3D", "CD3E", "CD8A", "CD4", "IL7R"
  ),

  Exhausted_T_Cells = c(
    "PDCD1", "HAVCR2", "TIGIT", "LAG3", "CTLA4", "CXCL13"
  ),

  NK_Cells = c(
    "GNLY", "NKG7", "KLRD1", "GZMH", "NCAM1"
  ),

  B_Plasma_Cells = c(
    "CD19", "MS4A1", "IGHG1", "IGKC", "SDC1"
  ),

  Oligodendrocytes = c(
    "MBP", "MOG", "PLP1", "MAG", "MOBP"
  ),

  OPC = c(
    "OLIG1", "PDGFRA", "CSPG4", "SOX10"
  ),

  Astrocytes = c(
    "GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A3"
  ),

  Endothelial = c(
    "PECAM1", "VWF", "CLDN5", "ESAM", "CD34"
  ),

  Pericytes = c(
    "ACTA2", "PDGFRB", "RGS5", "NOTCH3"
  )
)

# =========================
# GBM tumour state genes
# Neftel et al. 2019
# =========================

gbm_states <- list(
  MES = c("CHI3L1", "CD44", "VIM", "FN1", "SERPINE1", "TGFBI"),
  AC  = c("GFAP", "AQP4", "S100B", "AGT", "ALDOC"),
  OPC = c("OLIG1", "OLIG2", "PDGFRA", "CSPG4", "SOX10"),
  NPC = c("SOX2", "SOX4", "DLL3", "ASCL1", "HES6")
)

# =========================
# T cell exhaustion genes
# =========================

exhaustion_genes <- c(
  "PDCD1", "HAVCR2", "TIGIT", "LAG3",
  "CTLA4", "ENTPD1", "TOX", "NR4A1",
  "BATF", "CXCL13"
)