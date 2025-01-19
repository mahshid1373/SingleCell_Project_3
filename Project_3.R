set.seed(42)
library(dplyr)
# library(spatstat.core)
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(SingleR)
# library(enrichR)
library(CellChat)
library(SingleCellExperiment)
library(SeuratWrappers)
library(tidyverse)
# library(monocle3)
library(celldex)
library(hdf5r)
library(glmGamPoi)


# We wrote down the code in jupyter notebook and then copied it here. The code is written in R language.
################################################################################################
############################################# Week 1 ##########################################
################################################################################################
############################# 1.1 Properties of the Slides (1P)

#define directory with scRNAseq data
setwd("~/Downloads/project_3_dataset/")

# Define file paths for Section 1 and Section 2
section1_path <- "Section_1"

list.files(section1_path)

# Load spatial transcriptomics data
section1_data <- Load10X_Spatial(data.dir = section1_path, 
                                 filename = "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5")

print(section1_data)

### 1.2 Resolution of the spatial transcriptomics technology (1P)
average_cell_size <- 10
resolution_comment <- "The resolution of this technology is approximately 10 times larger than the size of an average eukaryotic cell. This means that each spot represents gene expression averaged over multiple cells."
cat(resolution_comment, "\n")

### 1.3 Output of Space Ranger (1P)
# Extract gene expression data
gene_expression_section1 <- GetAssayData(section1_data, assay = "Spatial", slot = "counts")

print(gene_expression_section1)

# Display metadata for Section 1
metadata_section1 <- section1_data@meta.data
# head(metadata_section1)

## 2 Spatial transcriptomics data in Seurat (3P)
### 2.1 Loading the data
#define directory with scRNAseq data
setwd("~/Downloads/project_3_dataset/")

# Define file paths for Section 1 and Section 2
section1_path <- "Section_1"

list.files(section1_path)

# Load spatial transcriptomics data
section1_data <- Load10X_Spatial(data.dir = section1_path, 
                                 filename = "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5")

print(section1_data)

### 2.2 Inspecting the Seurat-object
### 2.3 Visualization of a feature
DefaultAssay(section1_data) <- "Spatial"
# slotNames(section1_data@assays$Spatial)
section1_data <- NormalizeData(section1_data, assay = "Spatial")
# Visualize the expression of random genes for Section 1
SpatialFeaturePlot(section1_data, features = c("Gm1992", "Rp1"))

VlnPlot(section1_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

# section1_data <- subset(section1_data, subset = nFeature_Spatial > 2000 & nFeature_Spatial < 7500)
section1_data_Filtered <- subset(section1_data, 
                                 subset = nFeature_Spatial > 2500 & nFeature_Spatial < 7500 & 
                                 nCount_Spatial < 60000)

# Visualize gene count and mitochondrial content distributions (Section 1)
VlnPlot(section1_data_Filtered, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

### 3.2 Apply SCTransform
# Apply SCTransform for normalization and variance stabilization
section1_data_Filtered <- SCTransform(section1_data_Filtered, assay = "Spatial", verbose = FALSE)

## 4 Dimensionality reduction, clustering, and visualization (3P) 
### 4.1 Dimensionality reduction
# Perform PCA
section1_data_Filtered <- RunPCA(section1_data_Filtered, assay = "SCT", verbose = FALSE)

# Visualize variance explained by PCs to decide the number of dimensions to retain
ElbowPlot(section1_data_Filtered, ndims = 50)

### 4.2 Clustering
# Choose the number of dimensions based on the Elbow Plot (e.g., 10 PCs)
dims_to_use <- 1:10

# Run UMAP using the first 10 PCs (adjust based on ElbowPlot)
section1_data_Filtered <- RunUMAP(section1_data_Filtered, dims = dims_to_use)

# Find neighbors and clusters
section1_data_Filtered <- FindNeighbors(section1_data_Filtered, dims = dims_to_use, verbose = FALSE)
section1_data_Filtered <- FindClusters(section1_data_Filtered, resolution = 0.5, verbose = FALSE)

# View clusters in UMAP space
DimPlot(section1_data_Filtered, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Clustering on UMAP") +
  theme_minimal()
  
# Plot clusters on the tissue slide
SpatialDimPlot(section1_data_Filtered, group.by = "seurat_clusters", crop = TRUE) +
  ggtitle("Clusters on Tissue Slide") +
  theme_minimal()

# Save normalized data
saveRDS(section1_data_Filtered, file = "section1_normalized_data.rds")


################################################################################################
############################################# Week 2 ##########################################
################################################################################################

## 5 Differential Expression Analysis (8P)
# head(section1_data_Filtered@meta.data)
# Extract the specific cell barcode
specific_sample <- subset(section1_data_Filtered, nCount_Spatial = 9195)

### 5.1 DEG (Differentially Expressed Genes) analysis based on clustering
# Set the active identity to clusters
Idents(section1_data_Filtered) <- "seurat_clusters"

# Perform differential expression analysis between two clusters (e.g., cluster 1 and cluster 2)
cluster_markers <- FindMarkers(section1_data_Filtered, ident.1 = 1, ident.2 = 2, min.pct = 0.25, logfc.threshold = 0.25)

# View the top differentially expressed genes
# head(cluster_markers)

# Save the DEG results to a CSV file for use in the next task
write.csv(cluster_markers, file = "differentially_expressed_genes.csv")

### 5.2 DEG analysis based on the spatial patterning
cluster_markers <- FindAllMarkers(
  section1_data_Filtered,
  assay = "Spatial",
  test.use = "wilcox",        # Wilcoxon rank-sum test for DEG analysis
  min.pct = 0.25,             # Minimum fraction of spots expressing a gene
  logfc.threshold = 0.25      # Minimum log-fold change
)

# View the top markers
# head(cluster_markers)

# Save differentially expressed genes (DEGs) to a CSV file
write.csv(cluster_markers, file = "~/Documents/GitHub/SingleCell_Project_3/DEGs_by_clusters.csv")

# Visualize top DEGs for each cluster
top_markers <- cluster_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

# Heatmap for top markers across clusters
DoHeatmap(section1_data_Filtered, features = top_markers$gene) +
  ggtitle("Top DEGs by Cluster") +
  theme_minimal()


## 6 Merging the data (7P)
### 6.1 Merging without Batch-correction (3P)
# Check the header of the tissue_positions_list.csv file
coordinates <- Seurat::Read10X_Coordinates("~/Downloads/project_3_dataset/Section_1/spatial/tissue_positions_list.csv", filter.matrix = TRUE)

# is.numeric(coordinates)

# head(coordinates)

# Check the header of the tissue_positions_list.csv file
coordinates <- Seurat::Read10X_Coordinates("~/Downloads/project_3_dataset/Section_2/spatial/tissue_positions_list.csv", filter.matrix = TRUE)

# is.numeric(coordinates)

# head(coordinates)

#define directory with scRNAseq data
setwd("~/Downloads/project_3_dataset/")

# Define file paths for Section 1 and Section 2
section1_path <- "Section_1"
section2_path <- "Section_2"

list.files(section1_path)
list.files(section2_path)

# Load spatial transcriptomics data
section1_data <- Load10X_Spatial(data.dir = section1_path, 
                                 filename = "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5")

section2_data <- Load10X_Spatial(data.dir = section2_path, 
                                 filename = "V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5")
                 
# Combine section1_data and section2_data
combined_data <- merge(section1_data, y = section2_data, add.cell.ids = c("Section1", "Section2"))

# Print the combined data
print(combined_data)

# Repeat the necessary pre-processing,
# dimensionality reduction, and clustering steps for the merged dataset

# Preprocessing
combined_data <- SCTransform(combined_data, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction
combined_data <- RunPCA(combined_data, verbose = FALSE)
combined_data <- RunUMAP(combined_data, dims = dims_to_use, verbose = FALSE)

# Clustering
combined_data <- FindNeighbors(combined_data, dims = dims_to_use)
combined_data <- FindClusters(combined_data, resolution = 0.5)

# Load the repr package
library(repr)

# Set the plot width (e.g., 15 inches)
options(repr.plot.width = 15)

# Show the clustering in the UMAP and tissue slides
p1 <- DimPlot(combined_data, reduction = "umap", group.by = "orig.ident", label = TRUE)+ ggtitle("UMAP Ident Combination of Datasets")


# Visualize UMAP
p2 <- DimPlot(combined_data, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("UMAP Clustering Combination of Datasets")

# Display the plots side by side
# p1+p2

SpatialDimPlot(combined_data, images = "slice1", label = TRUE)

section1_data <- SCTransform(section1_data, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction
section1_data <- RunPCA(section1_data, verbose = FALSE)
section1_data <- RunUMAP(section1_data, dims = dims_to_use, verbose = FALSE)

# Clustering
section1_data <- FindNeighbors(section1_data, dims = dims_to_use)
section1_data <- FindClusters(section1_data, resolution = 0.5)

DimPlot(section1_data, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("UMAP Clustering Dataset 1")

section2_data <- SCTransform(section2_data, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction
section2_data <- RunPCA(section2_data, verbose = FALSE)
section2_data <- RunUMAP(section2_data, dims = dims_to_use, verbose = FALSE)

# Clustering
section2_data <- FindNeighbors(section2_data, dims = dims_to_use)
section2_data <- FindClusters(section2_data, resolution = 0.5)

DimPlot(section2_data, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("UMAP Clustering Dataset 2")

# Identify clusters in both samples
clusters_dataset1 <- unique(Idents(section1_data))
clusters_dataset2 <- unique(Idents(section2_data))

clusters_both <- intersect(clusters_dataset1, clusters_dataset2)

clusters_unique_dataset1 <- setdiff(clusters_dataset1, clusters_dataset2)
clusters_unique_dataset2 <- setdiff(clusters_dataset2, clusters_dataset1)

# Print comparison
cat("Clusters present in both samples:", clusters_both, "\n")
cat("Clusters unique to dataset1:", clusters_unique_dataset1, "\n")
cat("Clusters unique to dataset2:", clusters_unique_dataset2, "\n")

### 6.2 Merging with Batch-correction (3P)
#define directory with scRNAseq data
setwd("~/Downloads/project_3_dataset/")

# Define file paths for Section 1 and Section 2
section1_path <- "Section_1"
section2_path <- "Section_2"

list.files(section1_path)
list.files(section2_path)

# Load spatial transcriptomics data
section1_data <- Load10X_Spatial(data.dir = section1_path, 
                                 filename = "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5")

section2_data <- Load10X_Spatial(data.dir = section2_path, 
                                 filename = "V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5")
                 
                 
list <- list(section1 = section1_data, section2 = section2_data)

# Find integration features
features <- SelectIntegrationFeatures(object.list = list)
list <- PrepSCTIntegration(object.list = list, 
                                  anchor.features = features)

# Find anchors and integrate
anchors <- FindIntegrationAnchors(object.list = list, 
                                  normalization.method = "SCT",
                                  anchor.features = features)

integrated_data <- IntegrateData(anchorset = anchors, 
                                 normalization.method = "SCT")

# Process integrated data
integrated_data <- RunPCA(integrated_data) %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.8)



### 6.3 Detection of Batch-effects (1P)
p5 <- DimPlot(combined_data, group.by = "orig.ident") + 
  ggtitle("Without Integration")
# p6 <- DimPlot(combined_data, group.by = "orig.ident") + 
#   ggtitle("With Integration")
p5 + p6

# Calculate silhouette scores for both methods
library(cluster)

# Function to calculate silhouette score
calc_silhouette <- function(seurat_obj) {
  umap_coords <- Embeddings(seurat_obj, "umap")
  batch_labels <- as.numeric(factor(seurat_obj$orig.ident))
  si <- silhouette(batch_labels, dist(umap_coords))
#   mean(si[, "sil_width"])
}

# Calculate scores
sil_score_merged <- calc_silhouette(combined_data)
# sil_score_integrated <- calc_silhouette(integrated_data)

print(paste("Silhouette score without integration:", sil_score_merged))
# print(paste("Silhouette score with integration:", sil_score_integrated))



################################################################################################
############################################# Week 3 ##########################################
################################################################################################
## 7 Cell-type identification (8P)
### 7.1 Automatic Annotation using Data Integration (5P)
# refrence data
cortex.data <- readRDS(file = "~/Downloads/project_3_dataset/allen_cortex.rds")

#define directory with scRNAseq data
setwd("~/Downloads/project_3_dataset/")

# Define file paths for Section 1 and Section 2
section1_path <- "Section_1"
section2_path <- "Section_2"

list.files(section1_path)
list.files(section2_path)

# Load spatial transcriptomics data
section1_data <- Load10X_Spatial(data.dir = section1_path, 
                                 filename = "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5")

section2_data <- Load10X_Spatial(data.dir = section2_path, 
                                 filename = "V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5")
                 
# Combine section1_data and section2_data
combined_data <- merge(section1_data, y = section2_data, add.cell.ids = c("Section1", "Section2"))

# Repeat the necessary pre-processing,
# dimensionality reduction, and clustering steps for the merged dataset

# Preprocessing
combined_data <- SCTransform(combined_data, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction
combined_data <- RunPCA(combined_data, verbose = FALSE)
combined_data <- RunUMAP(combined_data, dims = dims_to_use, verbose = FALSE)

# Clustering
combined_data <- FindNeighbors(combined_data, dims = dims_to_use)
combined_data <- FindClusters(combined_data, resolution = 0.5)


# Check available assays in the cortex.data object
# cortex.data

# Preprocessing the reference dataset
cortex.data <- SCTransform(cortex.data, assay = "RNA", verbose = FALSE)
cortex.data <- RunPCA(cortex.data, verbose = FALSE)
cortex.data <- RunUMAP(cortex.data, dims = dim_to_use, verbose = FALSE)
cortex.data <- FindNeighbors(cortex.data, dims = dim_to_use)
cortex.data <- FindClusters(cortex.data, resolution = 0.5)

# Preprocessing the combined_data transcriptomics dataset
# combined_data <- SCTransform(combined_data, assay = "Spatial", verbose = FALSE)
# combined_data <- RunPCA(combined_data, verbose = FALSE)
# Find anchors for label transfer
# anchors <- FindTransferAnchors(reference = cortex.data, query = combined_data, normalization.method = "SCT")

# # Transfer cell type labels from reference to combined_data dataset
# combined_data <- MapQuery(
#   anchorset = anchors, 
#   query = combined_data, 
#   reference = cortex.data, 
#   refdata = Idents(cortex.data), 
#   reference.reduction = "pca"
# )

# # Run UMAP on combined_data data using reference anchors
# combined_data <- RunUMAP(combined_data, reduction = "ref.pca", dims = 1:30)

# # Plot the UMAP with predicted labels
# DimPlot(combined_data, reduction = "umap", group.by = "predicted.id", label = TRUE) +
#   ggtitle("UMAP with Transferred Cell Type Annotations")









