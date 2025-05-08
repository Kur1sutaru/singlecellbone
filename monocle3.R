library("Seurat")
library("remotes")
library(dplyr)
library(monocle3)
library("SeuratWrappers")
# first we will try sample bl_68774wtannoscpred.rds
# Convert Seurat object to Monocle3's cell_data_set
cds <- as.cell_data_set(periosteum_seurat)

# Add cluster info from scpred_prediction
cds@colData$cluster <- periosteum_seurat@meta.data$scpred_prediction

# OPTIONAL: Add cell type or condition annotations if desired
# cds@colData$condition <- periosteum_seurat@meta.data$condition_column_name

# Preprocess the data
cds <- preprocess_cds(cds, num_dim = 50)


# 1. Preprocessing (optional if already done)
cds <- preprocess_cds(cds, num_dim = 50)

# 2. Dimensionality reduction
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# ✅ 3. REQUIRED: Run Monocle clustering (so learn_graph will work)
cds <- cluster_cells(cds, reduction_method = "UMAP")

# 4. Overwrite clusters with your scpred annotations
cds@clusters$UMAP$clusters <- cds@colData$cluster

# 5. Now this will work!
cds <- learn_graph(cds)


# Step 4: Choose a root cluster and order cells
table(cds@colData$cluster)  # Check cluster names
root_cells <- colnames(cds)[cds@colData$cluster == "Osteogenic Progenitors"]
cds <- order_cells(cds, root_cells = root_cells[1])



# Plot trajectory colored by your scpred clusters
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE)

plot_cells(cds,
           color_cells_by = "cluster",        # This refers to cds@colData$cluster
           label_groups_by_cluster = TRUE,   # Disable floating cluster labels
           label_leaves = FALSE,
           label_branch_points = FALSE,
           show_trajectory_graph = TRUE)      # Show the trajectory lines


# Optional: Plot by pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE)

# Identify genes that change along the trajectory
deg_pseudotime <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

# Filter significant genes (e.g., q_value < 0.05)
deg_pseudotime_sig <- deg_pseudotime %>% filter(q_value < 0.05)

# Save to CSV
write.csv(deg_pseudotime_sig,
          file = "pseudotime_DEGs_qval0.05.csv",
          row.names = FALSE)
# Assuming your gene names are stored in the row names
rowData(cds)$gene_short_name <- rownames(cds)

# OR, if your gene names are stored in a column called 'id' or 'gene_name'
rowData(cds)$gene_short_name <- rownames(cds)


# Plot selected genes across pseudotime
plot_genes_in_pseudotime(cds[c("Sp7", "Runx2"), ])

# OR top variable genes:
top_genes <- rownames(deg_pseudotime_sig)[1:10]
plot_genes_in_pseudotime(cds[top_genes, ])


earlymesenc_genes <- c("Ly6a",
                    "Pdgfra",
                    "Cd44",
                    "Thy1",
                    "Nes",
                    "Lepr")

plot_cells(cds,
           genes=earlymesenc_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

matureosteo_genes <- c("Bglap",
                       "Spp1",
                       "Dmp1",
                       "Ibsp",
                       "Opn",
                       "Opg")

plot_cells(cds,
           genes=matureosteo_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


chondro_genes <- c("Sox9",
                       "Acan",
                       "Col2a1",
                       "Col10a1",
                       "Col3a1",
                       "Dcn")

plot_cells(cds,
           genes=chondro_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

skelprog_genes <- c("Cd82",
                   "Cd201",
                   "Cd164",
                   "Sox9",
                   "Grem1",
                   "Flk1",
                   "Ska1",
                   "Cxcl12",
                   "Pthlh")

plot_cells(cds,
           genes=skelprog_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)



cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "scpred_prediction",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
plot_cells(cds,
           color_cells_by = "scpred_prediction",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=2.5)



plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=2.5)


skelprog_genes <- c("Cd82",
                    "Cd201",
                    "Cd164",
                    "Sox9",
                    "Grem1",
                    "Flk1",
                    "Ska1",
                    "Cxcl12",
                    "Pthlh")
cds_subset <- cds[rowData(cds)$gene_short_name %in% skelprog_genes,]
gene_fits <- fit_models(cds_subset, model_formula_str = "~scpred_prediction")
fit_coefs <- coefficient_table(gene_fits)
str(fit_coefs)
library(dplyr)

fit_coefs <- coefficient_table(gene_fits)

degs <- fit_coefs %>%
  filter(term != "(Intercept)", !is.na(q_value), q_value < 0.05) %>%
  arrange(q_value)

# Show top results
head(degs)



write.csv(degs, "DEGs_scpred_prediction.csv", row.names = FALSE)

library(dplyr)

top_genes <- degs %>%
  group_by(term) %>%
  slice_min(order_by = q_value, n = 5) %>%
  pull(gene_short_name) %>%
  unique()

library(monocle3)

plot_genes_by_group(
  cds,
  genes = top_genes,
  group_cells_by = "scpred_prediction",
  ordering_type = "cluster_row",
  max.size = 3
)



plot_genes_violin(cds_subset, group_cells_by="scpred_prediction", ncol=2)

library(ggplot2)

plot_genes_violin(
  cds_subset,
  group_cells_by = "scpred_prediction",
  ncol = 2
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels


plot_genes_hybrid(cds_subset, group_cells_by="scpred_prediction", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


# Your gene set
skelprog_genes <- c(
  "Cd82",  "Cd164", "Sox9", "Grem1",
   "Ska1", "Cxcl12", "Pthlh"
)

# Plot genes in pseudotime, colored by scpred_prediction
plot_genes_in_pseudotime(
  cds[skelprog_genes, ],
  color_cells_by = "scpred_prediction",
  min_expr = 0.1,
  label_by_short_name = TRUE,
  ncol = 3
)

library(monocle3)
library(ggplot2)
library(dplyr)

# Ensure gene short names are in place
rowData(cds)$gene_short_name <- rownames(cds)

# Your gene list
skelprog_genes <- c(
  "Cd82", "Cd201", "Cd164", "Sox9", "Grem1",
  "Flk1", "Ska1", "Cxcl12", "Pthlh"
)

# Extract pseudotime and expression data for these genes
plot_df <- plot_genes_in_pseudotime(
  cds[skelprog_genes, ],
  color_cells_by = "scpred_prediction",
  min_expr = 0.1,
  return_data = TRUE  # 👈 IMPORTANT
)

# Now plot with smoothed lines using ggplot2
ggplot(plot_df, aes(x = pseudotime, y = expression, color = color_by)) +
  geom_point(alpha = 0.3, size = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ gene_short_name, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(color = "scpred_prediction", y = "Expression", x = "Pseudotime")



