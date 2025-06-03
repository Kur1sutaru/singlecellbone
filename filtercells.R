library(dplyr)
library(ggplot2)
library(Seurat)
library(readxl)
library(SeuratDisk)


setwd("/Users/cristalvillalba/Desktop/Lee lab bioinfo/emilly/datawithcellqc")

# Convert to Seurat object (automatically handles HDF5)
bl_68774 <- LoadH5Seurat("bl_68774.h5seurat")
bl_68775 <- LoadH5Seurat("bl_68775.h5seurat")
# first using sample bl_68774 - wt sample
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
bl_68774[["percent.mt"]] <- PercentageFeatureSet(bl_68774, pattern = "^MT-")
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(bl_68774, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bl_68774, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## as single nuclei, we need to filter out 5% of mt genes
bl_68774 <- subset(bl_68774, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
bl_68774 <- NormalizeData(bl_68774, normalization.method = "LogNormalize", scale.factor = 33696)
bl_68774 <- FindVariableFeatures(bl_68774, selection.method = "vst", nfeatures = 33696)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(bl_68774), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(bl_68774)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## scaling the data
all.genes <- rownames(bl_68774)
bl_68774 <- ScaleData(bl_68774, features = all.genes)

#pca
bl_68774 <- RunPCA(bl_68774, features = VariableFeatures(object = bl_68774))
VizDimLoadings(bl_68774, dims = 1:2, reduction = "pca")
DimPlot(bl_68774, reduction = "pca") + NoLegend()
### cluster the cells
bl_68774 <- FindNeighbors(bl_68774, dims = 1:30)
bl_68774 <- FindClusters(bl_68774, resolution = 0.5)
bl_68774 <- RunUMAP(bl_68774, dims = 1:30)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(bl_68774, reduction = "umap")



# STEP 1: Define gene list manually (from your message, unique markers)
gene_list <- c(
  "Cd31", "CDH5", "Cdh5", "Emcn", "Erg", "Flt1", "Pecam1",  # Endothelial
  "Kit", "Cd19", "Cd247", "CD43", "CD44", "CD45", "Itgam", "Ptprc",  # Hematopoietic
  "Ccr2", "CD115", "CD11b", "Cd14", "Cd209a", "Cd68", "Cebpe", "Cx3cr1", 
  "F4/80", "Gr-1", "Lgals3", "Ly6c2", "Mpo", "Retnlg"  # Myeloid
)


# STEP 2: Filter to valid gene names in your data
gene_list <- unique(na.omit(gene_list))
genes_in_data <- intersect(gene_list, rownames(bl_68774))

# STEP 3: Identify cells expressing any of those genes (log-normalized expression > 0)
expr_mat <- GetAssayData(bl_68774, slot = "data")
cells_to_remove <- colnames(expr_mat)[colSums(expr_mat[genes_in_data, , drop = FALSE] > 0) > 0]

# STEP 4: Subset object to exclude those cells
seurat_filtered <- subset(bl_68774, cells = setdiff(colnames(bl_68774), cells_to_remove))

# STEP 5: Save the filtered Seurat object
saveRDS(seurat_filtered, "bl_68774_no_filtendohemmyeloid_cells.rds")

#STEP 6: add metadata from cells we want to exclude
# Label cells that express any of the genes
bl_68774$expresses_marker <- ifelse(colnames(bl_68774) %in% cells_to_remove, "Expresses Marker", "Other")

# Plot UMAP
DimPlot(bl_68774, group.by = "expresses_marker", pt.size = 0.5) + ggtitle("Cells expressing marker genes")


library(ggplot2)

# Count cells
cell_counts <- table(bl_68774$expresses_marker) %>% as.data.frame()
colnames(cell_counts) <- c("Category", "CellCount")

# Bar plot
ggplot(cell_counts, aes(x = Category, y = CellCount, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of Cells Before Filtering", x = "", y = "Cell Count") +
  theme(legend.position = "none")


seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:20)

DimPlot(seurat_filtered, pt.size = 0.5) + ggtitle("UMAP After Filtering")

#after filtering
# Count remaining cells
cell_counts_filtered <- data.frame(
  Category = "Remaining After Filtering",
  CellCount = ncol(seurat_filtered)
)

# Bar plot after filtering
ggplot(cell_counts_filtered, aes(x = Category, y = CellCount, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of Cells After Filtering", x = "", y = "Cell Count") +
  theme(legend.position = "none")

saveRDS(seurat_filtered, "bl_68774_no_filtendohemmyeloid_cells.rds")


#### same processing steps in sample bl_68775 Mut
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
bl_68775[["percent.mt"]] <- PercentageFeatureSet(bl_68775, pattern = "^MT-")
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(bl_68775, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bl_68775, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## as single nuclei, we need to filter out 5% of mt genes
bl_68775 <- subset(bl_68775, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
bl_68775 <- NormalizeData(bl_68775, normalization.method = "LogNormalize", scale.factor = 33696)
bl_68775 <- FindVariableFeatures(bl_68775, selection.method = "vst", nfeatures = 33696)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(bl_68775), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(bl_68775)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## scaling the data
all.genes <- rownames(bl_68775)
bl_68775 <- ScaleData(bl_68775, features = all.genes)

#pca
bl_68775 <- RunPCA(bl_68775, features = VariableFeatures(object = bl_68775))
VizDimLoadings(bl_68775, dims = 1:2, reduction = "pca")
DimPlot(bl_68775, reduction = "pca") + NoLegend()
### cluster the cells
bl_68775 <- FindNeighbors(bl_68775, dims = 1:30)
bl_68775 <- FindClusters(bl_68775, resolution = 0.5)
bl_68775 <- RunUMAP(bl_68775, dims = 1:30)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(bl_68775, reduction = "umap")



# STEP 1: Define gene list manually (from your message, unique markers)
gene_list <- c(
  "Cd31", "CDH5", "Cdh5", "Emcn", "Erg", "Flt1", "Pecam1",  # Endothelial
  "Kit", "Cd19", "Cd247", "CD43", "CD44", "CD45", "Itgam", "Ptprc",  # Hematopoietic
  "Ccr2", "CD115", "CD11b", "Cd14", "Cd209a", "Cd68", "Cebpe", "Cx3cr1", 
  "F4/80", "Gr-1", "Lgals3", "Ly6c2", "Mpo", "Retnlg"  # Myeloid
)


# STEP 2: Filter to valid gene names in your data
gene_list <- unique(na.omit(gene_list))
genes_in_data <- intersect(gene_list, rownames(bl_68775))

# STEP 3: Identify cells expressing any of those genes (log-normalized expression > 0)
expr_mat <- GetAssayData(bl_68775, slot = "data")
cells_to_remove <- colnames(expr_mat)[colSums(expr_mat[genes_in_data, , drop = FALSE] > 0) > 0]

# STEP 4: Subset object to exclude those cells
seurat_filtered <- subset(bl_68775, cells = setdiff(colnames(bl_68775), cells_to_remove))

# STEP 5: Save the filtered Seurat object
saveRDS(seurat_filtered, "bl_68775_no_filtendohemmyeloid_cells.rds")

#STEP 6: add metadata from cells we want to exclude
# Label cells that express any of the genes
bl_68775$expresses_marker <- ifelse(colnames(bl_68775) %in% cells_to_remove, "Expresses Marker", "Other")

# Plot UMAP
DimPlot(bl_68775, group.by = "expresses_marker", pt.size = 0.5) + ggtitle("Cells expressing marker genes")


library(ggplot2)

# Count cells
cell_counts <- table(bl_68775$expresses_marker) %>% as.data.frame()
colnames(cell_counts) <- c("Category", "CellCount")

# Bar plot
ggplot(cell_counts, aes(x = Category, y = CellCount, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of Cells Before Filtering", x = "", y = "Cell Count") +
  theme(legend.position = "none")


seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:20)

DimPlot(seurat_filtered, pt.size = 0.5) + ggtitle("UMAP After Filtering")

#after filtering
# Count remaining cells
cell_counts_filtered <- data.frame(
  Category = "Remaining After Filtering",
  CellCount = ncol(seurat_filtered)
)

# Bar plot after filtering
ggplot(cell_counts_filtered, aes(x = Category, y = CellCount, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of Cells After Filtering", x = "", y = "Cell Count") +
  theme(legend.position = "none")

saveRDS(seurat_filtered, "bl_68775_no_filtendohemmyeloid_cells.rds")
