setwd("D:/Baylor/Leelab/singlecell/scprednew")
library("scPred")
library("glue")
library("dplyr")
library("Seurat")
library("magrittr")
library("caret")
options(Seurat.object.assay.version = "v5")
## to avoid compatibility errors, install by devtools::install_github(repo="powellgenomicslab/scPred",  ref="9f407b7436f40d44224a5976a94cc6815c6e837f")
#Training the model
reference <- scPred::longbone_stromal 
query <- scPred::pbmc_2
reference <- longbone_stromal %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)
DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE)
reference <- getFeatureSpace(reference, "cell_type")
reference <- trainModel(reference)
get_probabilities(reference) %>% head()
get_scpred(reference)
plot_probabilities(reference)

### After training the model with the reference, lets predict the query - our data set
# Normalize the data
wtlb66089subset   <- NormalizeData(wtlb66089subset)

# Find variable features
wtlb66089subset   <- FindVariableFeatures(wtlb66089subset)

# Scale the data
wtlb66089subset   <- ScaleData(wtlb66089subset)

# Perform PCA
wtlb66089subset   <- RunPCA(wtlb66089subset)

# Check available assays
Assays(wtlb66089subset)

# Make sure "integrated" is the default assay for predictions
wtlb66089subset   <- scPredict(wtlb66089subset, threshold = 0.08, reference)
DimPlot(wtlb66089subset, group.by = "scpred_prediction", reduction = "scpred")
wtlb66089subset   <- RunUMAP(wtlb66089subset, reduction = "scpred", dims = 1:30)
DimPlot(wtlb66089subset, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
# Store the trained scPred model in the 'misc' slot of the reference Seurat object
wtlb66089subset@misc$scPred_model <- reference@tools$scPred
# Calculate proportions of cell types
cell_proportions <- table(wtlb66089subset@meta.data$scpred_prediction)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Count")
cell_proportions_df$Proportion <- cell_proportions_df$Count / sum(cell_proportions_df$Count)
# Load ggplot2 for plotting
library(ggplot2)

# Create the bar plot with cell counts as labels
ggplot(cell_proportions_df, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +  # Add cell count labels above each bar
  theme_minimal() +
  labs(x = "Predicted Cell Type", y = "Proportion", title = "Cell Type Proportions from scPred") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


saveRDS(wt66089subset, "wtlb66089subset_annotatedscpred.rds")
saveRDS(bl_68774 , "bl_68774_annotated.rds")
saveRDS(wtlb66089subset , "wtlb66089subset_annotated.rds")
saveRDS(wtlb66089subset, "wtlb66089subset_annotated.rds")


### Accessing classifiers
crossTab(wtlb66089subset  , "cell_type", "scpred_prediction")
crossTab(wtlb66089subset  , "cell_type", "scpred_prediction", output = "prop")
get_classifiers(reference)

# Each model can be normally treated using the caret enviroment. 
# For example, we can plot the performance resamples using the plot.train:
caret::plot.train(get_classifiers(reference)[["Chondrocytes"]])


# Calculate proportions of cell types
cell_proportions <- table(reference@meta.data$cell_type)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Count")
cell_proportions_df$Proportion <- cell_proportions_df$Count / sum(cell_proportions_df$Count)
# Load ggplot2 for plotting
library(ggplot2)

# Create the bar plot
ggplot(cell_proportions_df, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Predicted Cell Type", y = "Proportion", title = "Cell Type Proportions from scPred") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability



