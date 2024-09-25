# Load necessary libraries
library(Seurat)
library(Matrix)
library(hdf5r)
library(harmony)
library(patchwork)
library(dplyr)
library(tibble)
library(jsonlite)
library(ggplot2)


# Function to load and process a single HDF5 file
process_h5_file <- function(file_path, sample_name) {
  h5_file <- H5File$new(file_path, mode = "r")
  matrix_group <- h5_file[["matrix"]]
  
  # Load necessary datasets from the HDF5 file
  barcodes <- matrix_group[["barcodes"]][]
  data <- matrix_group[["data"]][]
  indices <- matrix_group[["indices"]][]
  indptr <- matrix_group[["indptr"]][]
  shape <- matrix_group[["shape"]][]
  
  # Reconstruct the sparse matrix
  expression_matrix <- sparseMatrix(
    i = indices + 1,  # R is 1-based, so we add 1 to the indices
    p = indptr,
    x = data,
    dims = shape
  )
  
  # Create a Seurat object for RNA data
  seurat_obj <- CreateSeuratObject(counts = expression_matrix, project = sample_name)
  
  # Add mitochondrial content percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Basic QC filtering (adjust thresholds based on data specifics)
  seurat_obj <- subset(
    seurat_obj,
    subset = nCount_RNA < 25000 & nCount_RNA > 1000 & percent.mt < 20
  )
  
  # SCTransform for RNA normalization
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
  
  # Run PCA for dimensionality reduction
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  
  # Add sample identifier metadata
  seurat_obj$orig.ident <- sample_name
  
  # Close the HDF5 file
  h5_file$close()
  
  return(seurat_obj)
}

# Step 1: Process each HDF5 file
m1 <- process_h5_file("~/Desktop/m1_cell_feature_matrix.h5", "m1")
m4 <- process_h5_file("~/Desktop/m4_cell_feature_matrix.h5", "m4")
i1 <- process_h5_file("~/Desktop/i1_cell_feature_matrix.h5", "i1")
i5 <- process_h5_file("~/Desktop/i5_cell_feature_matrix.h5", "i5")

# Step 2: Merge the Seurat objects
combined_seurat <- merge(m1, y = list(m4, i1, i5), add.cell.ids = c("m1", "m4", "i1", "i5"))

# Step 3: Run SCTransform on the merged object
combined_seurat <- SCTransform(combined_seurat, verbose = FALSE)

# Step 4: Run PCA on the merged data
combined_seurat <- RunPCA(combined_seurat, verbose = FALSE)

# Step 5: Run Harmony integration
combined_seurat <- RunHarmony(
  object = combined_seurat,
  group.by.vars = "orig.ident",
  dims = 1:30
)

# Step 6: Run UMAP on Harmony integrated dimensions
combined_seurat <- RunUMAP(combined_seurat, reduction = "harmony", dims = 1:30)

# Step 7: Visualize the UMAP with RNA data
DimPlot(combined_seurat, reduction = "umap", group.by = "orig.ident", label = TRUE) + ggtitle("UMAP - RNA (Harmony)")

# Load the JSON file containing the gene mappings (adjust the file path as needed)
gene_info <- fromJSON("~/Library/CloudStorage/Box-Box/Kevin-Engler Project/20240802__164811__Kevinset1-mBrain/output-XETG00154__0015892__M4_0015892__20240802__164829/gene_panel.json")

# Step 8: Create a mapping between feature numbers and gene names
feature_to_gene <- lapply(1:nrow(gene_info$payload$targets), function(i) {
  feature_name <- paste0("Feature", i)  # Create "FeatureXXX"
  
  # Extract the gene name from the nested structure in the 'targets' field
  gene_name <- gene_info$payload$targets$type$data$name[i]  # Access the 'name' field for each gene
  
  return(c(feature_name, gene_name))
})

# Convert the list into a data frame
feature_to_gene_df <- as.data.frame(do.call(rbind, feature_to_gene), stringsAsFactors = FALSE)

# Rename the columns for clarity
colnames(feature_to_gene_df) <- c("Feature", "Gene")

# Step 9: Differential Expression Analysis
next_combined <- combined_seurat

# Step 9.1: Set identifiers for mock and infected samples
next_combined$group <- ifelse(next_combined$orig.ident %in% c("m1", "m4"), "Mock", "Infected")
Idents(next_combined) <- "group"

# Step 9.2: Prepare SCT data for differential expression analysis
next_combined <- PrepSCTFindMarkers(next_combined)

# Step 9.3: Perform differential expression analysis (RNA)
de_genes <- FindMarkers(next_combined, ident.1 = "Mock", ident.2 = "Infected", assay = "SCT")

# Step 9.4: Map feature names to gene names
de_genes_with_names <- de_genes %>%
  rownames_to_column("Feature") %>%
  left_join(feature_to_gene_df, by = "Feature") %>%
  select(Gene, everything())  # Move 'Gene' column to the front

# Step 9.5: Export the results with gene names
write.csv(de_genes_with_names, file = "/Users/stellawroblewski/Desktop/spatial_named_differential_expression_results_mockvsinfected_with_genes.csv", row.names = FALSE)

# Step 10: Volcano Plot 
# Load necessary library

# Step 1: Prepare the data
# Assuming 'de_genes_with_names' is your differential expression data with 'avg_log2FC', 'p_val_adj', and 'Gene' columns

# Step 2: Create a new column for significance based on thresholds
de_genes_with_names <- de_genes_with_names %>%
  mutate(Significance = case_when(
    avg_log2FC > 1 & p_val_adj < 0.05 ~ "Upregulated",
    avg_log2FC < -1 & p_val_adj < 0.05 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Step 3: Create a volcano plot
ggplot(de_genes_with_names, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 2.5) +  # Plot points
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +  # Color points
  theme_minimal() +  # Minimal theme
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano Plot") +  # Axis labels
  theme(legend.position = "right") +  # Legend position
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Fold change threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # p-value threshold line
  geom_text(aes(label = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 1, Gene, "")), 
            check_overlap = FALSE, vjust = 1.5, size = 3)  # Add gene labels to all significant points





# Step 1: Read the spatial and multiome differential expression CSV files
spatial_diffexp <- read.csv("~/Desktop/spatial_named_differential_expression_results_mockvsinfected_with_genes.csv")
multiome_diffexp <- read.csv("~/Desktop/differential_expression_results_mockvsinfected.csv")

# Step 2: Ensure both datasets have a common 'Gene' column and remove any duplicates
# Remove duplicates in case the same gene appears multiple times
spatial_diffexp <- spatial_diffexp %>% distinct(Gene, .keep_all = TRUE)
multiome_diffexp <- multiome_diffexp %>% distinct(Gene, .keep_all = TRUE)

# Step 3: Rename columns in multiome to avoid conflicts when merging
multiome_diffexp <- multiome_diffexp %>%
  rename(
    p_val_multiome = p_val,
    avg_log2FC_multiome = avg_log2FC,
    pct.1_multiome = pct.1,
    pct.2_multiome = pct.2,
    p_val_adj_multiome = p_val_adj
  )

# Step 4: Perform a left join on the 'Gene' column
# This will keep all rows from the spatial dataset and bring in matching multiome data
merged_diffexp <- spatial_diffexp %>%
  left_join(multiome_diffexp, by = "Gene")

# Step 5: Inspect the merged data (optional)
head(merged_diffexp)


# Step 5: Save the merged data to a new CSV file
write.csv(merged_diffexp, file = "/Users/stellawroblewski/Desktop/merged_spatial_multiome_diffexp.csv", row.names = FALSE)