## match gene ontology terms with tgf beta

library(clusterProfiler)
# Subset significant genes from chondro.de using row names
# Use p val adj or just the p value?
significant_genes <- rownames(subset(WTvsHet, pvalue < 0.05))
# GO enrichment analysis using mouse gene annotations
go_enrichment <- enrichGO(
  gene = WTvsHet,
  OrgDb = org.Mm.eg.db,  # Use the mouse gene annotation database
  keyType = "SYMBOL",    # Assuming your gene identifiers are symbols
  ont = "ALL",            # 'BP' for Biological Process, 'MF' for Molecular Function, 'CC' for Cellular Component
  pAdjustMethod = "none",  # Benjamini-Hochberg correction
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# GO:0141091 tgf beta signalling

# View the top results
head(go_enrichment)
write.csv(go_enrichment, "trabsitioningcellswtxkogsearesults.csv")
barplot(go_enrichment, showCategory = 10, title = "GO Enrichment in Transitioning cells longbone wtxko")



# Subset significant genes from chondro.de with gene names in rownames
significant_genes <- rownames(subset(WTvsHet, padj < 0.05))
# GO enrichment using Ensembl gene IDs
go_enrichment <- enrichGO(
  gene = significant_genes,
  OrgDb = org.Mm.eg.db,  # Change to org.Mm.eg.db for mouse
  keyType = "ENSEMBL",   # Specify Ensembl as the ID type
  ont = "ALL",            # Biological Process; use 'MF' or 'CC' for other categories
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# View results
head(go_enrichment)
write.csv(go_enrichment, "gotermswtxhet.csv")
## Use the specific gene ontology term to see which genes of the DEGs are present in there
# Perform KEGG enrichment to find TGF-beta pathway genes in the significant genes
kegg_enrichment <- enrichKEGG(
  gene = significant_genes,
  organism = "mmu",  # Use 'hsa' for human, 'mmu' for mouse
  keyType = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

barplot(tgfbetafriends, showCategory = 10, title = "TGFbeta related pathways in WTvs Het")

# Extract TGF-beta pathway genes from the enrichment results
# Find the TGF-beta pathway (KEGG pathway ID: hsa04350 for human)
tgf_beta_genes <- kegg_enrichment[kegg_enrichment$ID == "mmu04350", "geneID"]
tgf_beta_gene_list <- unlist(strsplit(as.character(tgf_beta_genes), "/"))  # Split into individual genes

# Display TGF-beta genes
print(tgf_beta_gene_list)


# Define the list of Ensembl gene IDs
ensembl_genes <- row.names(WTvsHet)

# Convert Ensembl IDs to Entrez IDs using the org.Hs.eg.db database
conversion_table <- bitr(
  ensembl_genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db  # Use org.Mm.eg.db for mouse data
)

# View the conversion table
print(conversion_table)
# GO enrichment analysis with Ensembl IDs
go_enrichment <- enrichGO(
  gene = conversion_table$ENTREZID,  # Use Entrez IDs if required
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",  # Use "ENSEMBL" if directly working with Ensembl IDs
  ont = "BP",            # Biological Process; change to "MF" or "CC" if needed
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# View GO results
head(go_enrichment)

# KEGG enrichment analysis with Entrez IDs
kegg_enrichment <- enrichKEGG(
  gene = conversion_table$ENTREZID,
  organism = "mmu",  # "mmu" for mouse in KEGG
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# View KEGG results
head(kegg_enrichment)

write.csv(kegg_enrichment, "keggwtxhet.csv")

library(pathfindR)

# Run pathfindR for KEGG pathway enrichment
pathfindR_results <- run_pathfindR(
  WTvsHet,
  gene_sets = "KEGG"  # You can also specify "Reactome" or "GO-BP" as needed
)

# View the pathfindR results
head(pathfindR_results)

WTvsHet <-na.omit(WTvsHet)



library(ggplot2)

tabelinhaCristal <- read.csv("~/Documents/tabelinhaCristal.csv")
mydata<-tabelinhaCristal

####P da comparacao
#Por tamanho
ggplot(mydata,
       aes_string(
         x = 'Grade',
         y = "Gene",
         size = 'P.value'
       )
) +
  geom_point(colour="red") +
  ylab(NULL) + theme(axis.text.y = element_text(
    size = 12,
    angle = 0,
    hjust = 1,
    vjust = 0,
    face = "plain"
  ))+ scale_size(range = c(3, 8))


#Por cor
ggplot(mydata,
       aes_string(
         x = 'Grade',
         y = "Gene",
         color = 'P.value'
       )
) +
  geom_point(size=5) +
  scale_color_continuous(
    low = "red",
    high = "blue",
    name = 'P-value',
    guide = guide_colorbar(reverse = TRUE)
  ) +
  ylab(NULL) + theme(axis.text.y = element_text(
    size = 12,
    angle = 0,
    hjust = 1,
    vjust = 0,
    face = "plain"
  ))  

#########P do tecido#############
ggplot(mydata,
       aes_string(
         x = 'Grade',
         y = "Gene",
         size = 'P.value.1',
         colour="Histology"
       )
) +
  geom_point() +
  scale_colour_hue(l = 70, c = 150)+
  ylab(NULL) + theme(axis.text.y = element_text(
    size = 12,
    angle = 0,
    hjust = 1,
    vjust = 0,
    face = "plain"
  ))+ scale_size(name="P-value",range = c(3, 8))


dotplot(kegg_enrichment@result, showCategory = 10, title = "KEGG Pathway WTvsHet")
write.csv(kegg_enrichment@result, "wtvshetkeggresults.csv")


gsets_list <- get_gene_sets_list(
  source = "KEGG",
  org_code = "mmu"
)
mmu_kegg_genes <- gsets_list$gene_sets
mmu_kegg_descriptions <- gsets_list$descriptions

## Save both as RDS files for later use
saveRDS(mmu_kegg_genes, "mmu_kegg_genes.RDS")
saveRDS(mmu_kegg_descriptions, "mmu_kegg_descriptions.RDS")
mmu_kegg_genes <- readRDS("mmu_kegg_genes.RDS")
mmu_kegg_descriptions <- readRDS("mmu_kegg_descriptions.RDS")
## Downloading the STRING PIN file to tempdir
url <- "https://stringdb-static.org/download/protein.links.v11.0/10090.protein.links.v11.0.txt.gz"
path2file <- file.path(tempdir(check = TRUE), "STRING.txt.gz")
download.file(url, path2file)

## read STRING pin file
mmu_string_df <- read.table(path2file, header = TRUE)

## filter using combined_score cut-off value of 800
mmu_string_df <- mmu_string_df[mmu_string_df$combined_score >= 800, ]

## fix ids
mmu_string_pin <- data.frame(
  Interactor_A = sub("^10090\\.", "", mmu_string_df$protein1),
  Interactor_B = sub("^10090\\.", "", mmu_string_df$protein2)
)
head(mmu_string_pin, 2)
library(biomaRt)

mmu_ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

converted <- getBM(
  attributes = c("ensembl_peptide_id", "mgi_symbol"),
  filters = "ensembl_peptide_id",
  values = unique(unlist(mmu_string_pin)),
  mart = mmu_ensembl
)
mmu_string_pin$Interactor_A <- converted$mgi_symbol[match(mmu_string_pin$Interactor_A, converted$ensembl_peptide_id)]
mmu_string_pin$Interactor_B <- converted$mgi_symbol[match(mmu_string_pin$Interactor_B, converted$ensembl_peptide_id)]
mmu_string_pin <- mmu_string_pin[!is.na(mmu_string_pin$Interactor_A) & !is.na(mmu_string_pin$Interactor_B), ]
mmu_string_pin <- mmu_string_pin[mmu_string_pin$Interactor_A != "" & mmu_string_pin$Interactor_B != "", ]

head(mmu_string_pin, 2)


# remove self interactions
self_intr_cond <- mmu_string_pin$Interactor_A == mmu_string_pin$Interactor_B
mmu_string_pin <- mmu_string_pin[!self_intr_cond, ]

# remove duplicated inteactions (including symmetric ones)
mmu_string_pin <- unique(t(apply(mmu_string_pin, 1, sort))) # this will return a matrix object

mmu_string_pin <- data.frame(
  A = mmu_string_pin[, 1],
  pp = "pp",
  B = mmu_string_pin[, 2]
)

path2SIF <- file.path(tempdir(), "mmusculusPIN.sif")
write.table(mmu_string_pin,
            file = path2SIF,
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE
)
path2SIF <- normalizePath(path2SIF)

# now convert ensembl in gene symbol to run pathfindR
# Load biomaRt
library(biomaRt)

# Connect to the Ensembl database for mouse (change dataset for human if needed)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")  # "hsapiens_gene_ensembl" for human
# Extract Ensembl gene IDs from the data frame
ensembl_genes <- WTvsHet$gene
# Retrieve gene symbols for the Ensembl IDs
gene_conversion <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),  # Use "mgi_symbol" for mouse, "hgnc_symbol" for human
  filters = "ensembl_gene_id",
  values = ensembl_genes,
  mart = ensembl
)

# View the conversion table
print(gene_conversion)

# Merge the gene symbol data with WTvsHet
WTvsHet_with_symbols <- merge(WTvsHet, gene_conversion, by.x = "gene", by.y = "ensembl_gene_id", all.x = TRUE)

# View the data frame with gene symbols added
print(WTvsHet_with_symbols)

## save and run pathfindR
write.csv(WTvsHet_with_symbols, "wtvshetgenes.csv")
wtvshetgenes <- read.csv("~/Downloads/go_tgfbeta/wtvshetgenes.csv")
View(wtvshetgenes)


## Run pathfindR
example_mmu_output <- run_pathfindR(
  input = wtvshetgenes,
  gene_sets = "KEGG",
  custom_genes = mmu_kegg_genes,
  custom_descriptions = mmu_kegg_descriptions,
  convert2alias = TRUE,
  pin_name_path = path2SIF
)

knitr::kable(example_mmu_output)
# Assuming 'example_mmu_output' is the data frame or table you want to save
write.csv(example_mmu_output, "example_mmu_output.csv", row.names = FALSE)

cluster_enriched_terms(pathwaysresults, use_description = TRUE)

input_results<-input_processing(
  example_mmu_output,
  p_val_threshold = 0.05,
  pin_name_path = "Biogrid",
  convert2alias = TRUE
)

visualize_terms(example_mmu_output,
                input_processed = input_results,
                hsa_KEGG = FALSE,
                pin_name_path = "Biogrid")

# The table must be contain 2 or 3 columns (gene name, p value, and or fold change).
write.csv(pathwaysresults, "c11_resSig_clusters_cutoff0.64kegg.csv")

jpeg("enrichkegg.jpeg")
enrichment_chart(example_mmu_output, top_terms = 20)
dev.off()

term_gene_graph(example_mmu_output, num_terms = 3,use_description = TRUE)

# Plot the enrichment results
enrichment_chart(pathfindR_results)
# Visualize a specific pathway, e.g., "TGF-beta signaling pathway" from KEGG
visualize_hsa_KEGG(
  hsa_kegg_ids <- c("hsa04010", "hsa04310", "hsa04668", "hsa04350"),
  example_mmu_output,
  scale_vals = TRUE,
  node_cols = NULL,
  quiet = TRUE,
  key_gravity = "northeast",
  logo_gravity = "southeast"
)
# Ensure example_mmu_output contains fold change data in a column named "logFC"
# Make sure input_results is the processed output of the pathfindR analysis

# Visualize terms with up/down regulation based on fold change
visualize_terms(
  result_df = example_mmu_output,
  input_processed = input_results,       # Processed input from pathfindR
  hsa_KEGG = FALSE,                      # Set to TRUE if working with human data
  pin_name_path = "Biogrid",             # Interaction network to use (e.g., "Biogrid" for mouse)
  fc_col = "logFC"                       # Specify the fold change column to color code up/down regulation
)


# Create the ggplot volcano plot with gene names
library(ggplot2)

# Add -log10(p-value) to the data in chondro.de
WT$neg_log10_pval <- -log10(chondro.de$p_val_adj)

# Define significance criteria for labeling genes
significant_genes <- WTvsHet[WTvsHet$padj < 0.05 & abs(WTvsHet$log2FoldChange) > 0.25, ]

# Create the ggplot volcano plot with gene names
ggplot(WTvsHet, aes(x = log2FoldChange, y = neg_log10_pval)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "blue") +
  geom_text(
    data = significant_genes,
    aes(label = rownames(significant_genes)),
    size = 3,
    vjust = 1,       # Adjust vertical position
    hjust = 1,       # Adjust horizontal position
    check_overlap = TRUE # Avoid overlapping text
  ) +
  labs(
    title = "Volcano Plot: WT vs Het",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal()



# Load the CSV file
data <- read.csv("/Users/cristalvillalba/Downloads/go_tgfbeta/WTvsHetdegs.csv")  # Replace with your actual file path
# Create Enhanced Volcano plot
EnhancedVolcano(
  data,
  lab = data$X,                      # Use gene names or IDs as labels
  x = 'log2FoldChange',                 # Fold change column
  y = 'padj',                           # Adjusted p-value column
  title = 'Volcano Plot of Differential Expression',
  xlab = 'Log2 Fold Change',
  ylab = '-Log10 Adjusted P-value',
  pCutoff = 0.05,                       # p-value threshold
  FCcutoff = 1.0,                       # Fold change threshold
  pointSize = 2.0,                      # Point size
  labSize = 3.0,                        # Label size
  col = c("grey30", "forestgreen", "royalblue", "red2"), # Colors for significance categories
  legendLabels = c("NS", "Log2 FC", "Adjusted p-value", "Significant")  # Custom legend labels
)

# Connect to the Ensembl database for mouse (use "hsapiens_gene_ensembl" for human data)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")  # "mmusculus_gene_ensembl" for mouse

# Get the Ensembl IDs from the data
ensembl_genes <- data$X  # Replace with the actual column name if different

# Retrieve gene symbols for the Ensembl IDs
gene_conversion <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),  # Use "hgnc_symbol" for human
  filters = "ensembl_gene_id",
  values = ensembl_genes,
  mart = ensembl
)

# Merge the gene symbols with your data
data <- merge(data, gene_conversion, by.x = "X", by.y = "ensembl_gene_id", all.x = TRUE)

# Replace missing symbols with Ensembl IDs if any
data$mgi_symbol[is.na(data$mgi_symbol)] <- data$X[is.na(data$mgi_symbol)]

# Generate the Enhanced Volcano plot using gene symbols as labels
EnhancedVolcano(
  data,
  lab = data$mgi_symbol,                   # Use gene symbols for labels
  x = 'log2FoldChange',                    # Fold change column
  y = 'padj',                              # Adjusted p-value column
  title = 'Volcano Plot of WT vs Het',
  xlab = 'Log2 Fold Change',
  ylab = '-Log10 Adjusted P-value',
  pCutoff = 0.05,                          # p-value threshold
  FCcutoff = 1.0,                          # Fold change threshold
  pointSize = 2.0,                         # Point size
  labSize = 3.0,                           # Label size
  col = c("grey30", "forestgreen", "royalblue", "red2"), # Colors for significance categories
  legendLabels = c("NS", "Log2 FC", "Adjusted p-value", "Significant")  # Custom legend labels
)



