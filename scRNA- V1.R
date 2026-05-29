library(Seurat)
library(ggplot2)
library(openxlsx)

input_file <- "C:/anay/scripts/scRNA/input/20k_NSCLC_DTC_3p_nextgem_donor_1_count_sample_feature_bc_matrix.h5"

outdir <- "C:/anay/scripts/scRNA/output"

if (!dir.exists(outdir)) {
dir.create(outdir, recursive = TRUE)
}

cat("Starting analysis...\n")

data <- Read10X_h5(input_file)

counts <- data$`Gene Expression`

print(dim(counts))

obj <- CreateSeuratObject(
counts = counts,
project = "NSCLC",
min.cells = 3,
min.features = 200
)

print(obj)

obj[["percent.mt"]] <- PercentageFeatureSet(
obj,
pattern = "^MT-"
)

write.xlsx(
obj@meta.data,
file = paste0(outdir, "/metadata_QC.xlsx"),
rowNames = TRUE
)

jpeg(
paste0(outdir, "/QC_violin_plot.jpeg"),
width = 2000,
height = 1200,
res = 200
)

print(
VlnPlot(
obj,
features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
ncol = 3
)
)

dev.off()

jpeg(
paste0(outdir, "/feature_scatter.jpeg"),
width = 1800,
height = 1400,
res = 200
)

print(
FeatureScatter(
obj,
feature1 = "nCount_RNA",
feature2 = "nFeature_RNA"
) + geom_smooth(method = "lm")
)

dev.off()

obj <- subset(
obj,
subset =
nFeature_RNA > 200 &
nFeature_RNA < 2500 &
percent.mt < 5
)

write.xlsx(
obj@meta.data,
file = paste0(outdir, "/filtered_metadata.xlsx"),
rowNames = TRUE
)

obj <- NormalizeData(obj)

obj <- FindVariableFeatures(
obj,
selection.method = "vst",
nfeatures = 2000
)

top10 <- head(VariableFeatures(obj), 10)

write.xlsx(
data.frame(top10),
file = paste0(outdir, "/top10_variable_genes.xlsx"),
rowNames = FALSE
)

jpeg(
paste0(outdir, "/variable_features.jpeg"),
width = 1800,
height = 1400,
res = 200
)

plot1 <- VariableFeaturePlot(obj)

print(
LabelPoints(
plot = plot1,
points = top10,
repel = TRUE
)
)

dev.off()

all.genes <- rownames(obj)

obj <- ScaleData(
obj,
features = all.genes
)

obj <- RunPCA(
obj,
features = VariableFeatures(object = obj)
)

jpeg(
paste0(outdir, "/PCA_heatmap.jpeg"),
width = 1800,
height = 1400,
res = 200
)

DimHeatmap(
obj,
dims = 1,
cells = 500,
balanced = TRUE
)

dev.off()

jpeg(
paste0(outdir, "/elbow_plot.jpeg"),
width = 1600,
height = 1200,
res = 200
)

print(
ElbowPlot(obj)
)

dev.off()

obj <- FindNeighbors(
obj,
dims = 1:15
)

obj <- FindClusters(
obj,
resolution = 0.5
)

jpeg(
paste0(outdir, "/cluster_plot.jpeg"),
width = 1800,
height = 1400,
res = 200
)

print(
DimPlot(
obj,
group.by = "seurat_clusters",
label = TRUE
)
)

dev.off()

obj <- RunUMAP(
obj,
dims = 1:15
)

jpeg(
paste0(outdir, "/UMAP_plot.jpeg"),
width = 1800,
height = 1400,
res = 200
)

print(
DimPlot(
obj,
reduction = "umap",
label = TRUE
)
)

dev.off()

write.xlsx(
obj@meta.data,
file = paste0(outdir, "/final_metadata.xlsx"),
rowNames = TRUE
)

saveRDS(
obj,
file = paste0(outdir, "/nsclc_seurat_object.rds")
)

cat("Analysis completed successfully!\n")
