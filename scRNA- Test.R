library(Seurat)
library(ggplot2)
library(openxlsx)

input_file <- "C:/anay/scripts/scRNA/input/20k_NSCLC_DTC_3p_nextgem_donor_1_count_sample_feature_bc_matrix.h5"

outdir <- "C:/anay/scripts/scRNA/output"

if (!dir.exists(outdir)) {
dir.create(outdir, recursive = TRUE)
}

cat("Starting analysis...\n")

nsclc.sparse.m <- Read10X_h5(
filename = input_file
)

cts <- nsclc.sparse.m$`Gene Expression`

print(dim(cts))

nsclc.seurat.obj <- CreateSeuratObject(
counts = cts,
project = "NSCLC",
min.cells = 3,
min.features = 200
)

print(nsclc.seurat.obj)

nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(
nsclc.seurat.obj,
pattern = "^MT-"
)

write.xlsx(
[nsclc.seurat.obj@meta.data](mailto:nsclc.seurat.obj@meta.data),
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
nsclc.seurat.obj,
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
nsclc.seurat.obj,
feature1 = "nCount_RNA",
feature2 = "nFeature_RNA"
) + geom_smooth(method = "lm")
)

dev.off()

nsclc.seurat.obj <- subset(
nsclc.seurat.obj,
subset =
nFeature_RNA > 200 &
nFeature_RNA < 2500 &
percent.mt < 5
)

write.xlsx(
[nsclc.seurat.obj@meta.data](mailto:nsclc.seurat.obj@meta.data),
file = paste0(outdir, "/filtered_metadata.xlsx"),
rowNames = TRUE
)

nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)

nsclc.seurat.obj <- FindVariableFeatures(
nsclc.seurat.obj,
selection.method = "vst",
nfeatures = 2000
)

top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

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

plot1 <- VariableFeaturePlot(nsclc.seurat.obj)

print(
LabelPoints(
plot = plot1,
points = top10,
repel = TRUE
)
)

dev.off()

all.genes <- rownames(nsclc.seurat.obj)

nsclc.seurat.obj <- ScaleData(
nsclc.seurat.obj,
features = all.genes
)

nsclc.seurat.obj <- RunPCA(
nsclc.seurat.obj,
features = VariableFeatures(object = nsclc.seurat.obj)
)

jpeg(
paste0(outdir, "/PCA_heatmap.jpeg"),
width = 1800,
height = 1400,
res = 200
)

DimHeatmap(
nsclc.seurat.obj,
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
ElbowPlot(nsclc.seurat.obj)
)

dev.off()

nsclc.seurat.obj <- FindNeighbors(
nsclc.seurat.obj,
dims = 1:15
)

nsclc.seurat.obj <- FindClusters(
nsclc.seurat.obj,
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
nsclc.seurat.obj,
group.by = "seurat_clusters",
label = TRUE
)
)

dev.off()

nsclc.seurat.obj <- RunUMAP(
nsclc.seurat.obj,
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
nsclc.seurat.obj,
reduction = "umap",
label = TRUE
)
)

dev.off()

write.xlsx(
[nsclc.seurat.obj@meta.data](mailto:nsclc.seurat.obj@meta.data),
file = paste0(outdir, "/final_metadata.xlsx"),
rowNames = TRUE
)

saveRDS(
nsclc.seurat.obj,
file = paste0(outdir, "/nsclc_seurat_object.rds")
)

cat("Analysis completed successfully!\n")
[]