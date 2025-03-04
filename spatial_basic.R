##############################################################################
# 1) LOAD REQUIRED LIBRARIES
##############################################################################
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(tibble)
library(openxlsx)
library(HGNChelper)  # for sc-type annotation

# If any of these are missing, install them first:
# install.packages(c("Seurat", "openxlsx", "HGNChelper", "patchwork"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Mm.eg.db")

##############################################################################
# 2) SOURCE SC-TYPE FUNCTIONS (if using sc-type for annotation)
##############################################################################
# NOTE: Adjust these URLs or use a local copy if needed. 
# The sc-type code below uses an older function signature.
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

##############################################################################
# 3) FUNCTION TO PROCESS ONE XENIUM H5 FILE
##############################################################################
process_xenium_sample <- function(h5_path, sample_id) {
  # 3a) Read the H5 file
  counts_list <- Read10X_h5(h5_path)
  
  # 3b) Extract the "Gene Expression" sub-matrix
  gene_expr <- counts_list[["Gene Expression"]]
  
  # 3c) Create a Seurat object
  so <- CreateSeuratObject(
    counts = gene_expr,
    project = sample_id,
    assay = "RNA"
  )
  
  # 3d) Basic QC
  # Optionally calculate percent.mt if mitochondrial genes follow a pattern (e.g. '^MT-'):
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  
  # Adjust thresholds to suit your data. Example:
  so <- subset(so, subset = nCount_RNA > 500 & percent.mt < 15)
  
  # 3e) SCTransform
  so <- SCTransform(so, verbose = FALSE)
  
  # 3f) PCA + UMAP
  so <- RunPCA(so, verbose = FALSE)
  # Use "umap" as the reduction name to easily reference later (DimPlot(..., reduction="umap"))
  so <- RunUMAP(so, dims = 1:30, reduction.name = "umap", reduction.key = "UMAP_")

  # Keep track of sample identity
  so$orig.ident <- sample_id
  
  return(so)
}

##############################################################################
# 4) PROCESS MULTIPLE SAMPLES (EXAMPLE)
##############################################################################
# Replace these file paths with your own Xenium H5 paths:
# sample1_path <- "path/to/sample1/cell_feature_matrix.h5"
# sample2_path <- "path/to/sample2/cell_feature_matrix.h5"
# etc.

# sample1_obj <- process_xenium_sample(sample1_path, "sample1")
# sample2_obj <- process_xenium_sample(sample2_path, "sample2")

# If you have multiple samples, you can merge them:
# combined <- merge(sample1_obj, y = sample2_obj, add.cell.id = c("S1", "S2"))

# For demonstration, we'll assume we have a single processed object called `combined`.
# If you only have one sample, just rename that object to "combined".

# combined <- sample1_obj  # If you only processed a single sample

##############################################################################
# 5) OPTIONAL: MERGE AND RE-INTEGRATE DATA
##############################################################################
# If you merged multiple samples:
#   1) Set default assay to RNA
#   2) Re-run SCTransform, PCA, UMAP, etc. on the merged object
#
# Example:
# DefaultAssay(combined) <- "RNA"
# combined <- SCTransform(combined, verbose = FALSE)
# combined <- RunPCA(combined, verbose = FALSE)
# combined <- RunUMAP(combined, dims = 1:30, reduction.name = "umap", reduction.key = "UMAP_")
#
# Quick check:
# DimPlot(combined, reduction = "umap", group.by = "orig.ident", label = TRUE)

##############################################################################
# 6) DEFINE GROUPS & RUN DIFFERENTIAL EXPRESSION (EXAMPLE)
##############################################################################
# If you have two experimental groups, define them here. Adjust as needed.
# Example: "GroupA" vs. "GroupB"
# combined$group <- ifelse(combined$orig.ident == "sample1", "GroupA", "GroupB")
#
# Then set the identity to "group" for differential expression:
# Idents(combined) <- "group"
#
# Prepare for DE (required if using SCTransform data):
# combined <- PrepSCTFindMarkers(combined)
#
# Run a Wilcoxon test across all cells in GroupA vs. GroupB:
# de_results <- FindMarkers(
#   object = combined,
#   ident.1 = "GroupA",
#   ident.2 = "GroupB",
#   assay = "SCT",           # or "RNA" if not using SCTransform
#   logfc.threshold = 0,     # no FC threshold, adjust as needed
#   min.pct = 0,             # no pct threshold, adjust as needed
#   test.use = "wilcox"
# )
# de_results <- rownames_to_column(de_results, var = "gene")
#
# Optionally add average expression in each group:
# avg_expr_list <- AverageExpression(combined, group.by = "group", assays = "SCT", slot = "data")
# avg_expr_df <- as.data.frame(avg_expr_list$SCT)
# avg_expr_df <- rownames_to_column(avg_expr_df, var = "gene")
# colnames(avg_expr_df) <- c("gene", "avg_expr_GroupA", "avg_expr_GroupB")
#
# final_df <- left_join(de_results, avg_expr_df, by = "gene") %>% arrange(p_val_adj)
# write.csv(final_df, "GroupA_vs_GroupB_DE_all_genes.csv", row.names = FALSE)

##############################################################################
# 7) CELL-TYPE ANNOTATION WITH SC-TYPE (OLDER VERSION USAGE)
##############################################################################
# sc-type uses a database of marker genes to predict cell types.
# Here is an example pipeline. Adjust "Brain" to another tissue if appropriate.

# Prepare gene sets from sc-type’s database:
# gs_list <- gene_sets_prepare(
#   "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",
#   "Brain" # e.g. "Brain", "Immune", "Pancreas", etc.
# )
#
# Switch to the SCT assay if you used SCTransform, then extract scaled data:
# DefaultAssay(combined) <- "SCT"
# scRNAseqData <- as.matrix(combined[["SCT"]]@scale.data)
#
# Score cell types:
# es.max <- sctype_score(
#   scRNAseqData = scRNAseqData,
#   scaled = TRUE,
#   gs = gs_list$gs_positive,
#   gs2 = gs_list$gs_negative
# )
#
# Assign each cell the highest-scoring cell type:
# combined$cell_type <- apply(es.max, 2, function(x) rownames(es.max)[which.max(x)])
#
# Example UMAP plot colored by sc-type results:
# DimPlot(combined, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE)

##############################################################################
# 8) SUBSETTING A SPECIFIC CELL TYPE OF INTEREST (EXAMPLE)
##############################################################################
# Suppose you want to isolate "Microglial cells" (or any labeled cell type).
# We'll subset, re-run SCTransform, cluster, and identify marker genes.

# microglia_obj <- subset(combined, subset = cell_type == "Microglial cells")
#
# DefaultAssay(microglia_obj) <- "RNA"
# microglia_obj <- SCTransform(microglia_obj, verbose = FALSE)
#
# microglia_obj <- RunPCA(microglia_obj, verbose = FALSE)
# microglia_obj <- FindNeighbors(microglia_obj, dims = 1:20)
# microglia_obj <- FindClusters(microglia_obj, resolution = 0.5)  # adjust as needed
# microglia_obj <- RunUMAP(microglia_obj, dims = 1:20, 
#                          reduction.name = "umap_clusters",
#                          reduction.key = "UMAPcl_")
#
# DimPlot(microglia_obj, reduction = "umap_clusters", group.by = "seurat_clusters", label = TRUE)
#
# # Contingency table of cluster vs. group (if relevant):
# table(microglia_obj$seurat_clusters, microglia_obj$group)
#
# # Identify marker genes for each microglial cluster
# microglia_obj <- PrepSCTFindMarkers(microglia_obj)
# all_markers <- FindAllMarkers(
#   microglia_obj,
#   assay = "SCT",
#   only.pos = TRUE,
#   logfc.threshold = 0.25,
#   min.pct = 0.1
# )
# write.csv(all_markers, "microglia_cluster_markers.csv", row.names = FALSE)

##############################################################################
# 9) GO/KEGG ENRICHMENT (EXAMPLE WITH clusterProfiler)
##############################################################################
# Here’s a brief template for running GO/KEGG using clusterProfiler.
# Adjust p-value/FC thresholds and the organism DB as needed.

# library(clusterProfiler)
# library(org.Mm.eg.db)  # for mouse, use org.Hs.eg.db for human
#
# # Example: Suppose 'final_df' is a table of DE results with columns:
# #          "gene", "p_val_adj", "avg_log2FC"
# sig_genes <- final_df %>%
#   filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
#   pull(gene)
#
# # Convert SYMBOL -> ENTREZ ID
# gene_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# entrez_ids <- unique(gene_ids$ENTREZID)
#
# # GO enrichment (Biological Process)
# go_results <- enrichGO(
#   gene = entrez_ids,
#   OrgDb = org.Mm.eg.db,
#   keyType = "ENTREZID",
#   ont = "BP",
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.05,
#   qvalueCutoff = 0.05
# )
# go_results_df <- as.data.frame(go_results)
# write.csv(go_results_df, "GO_Enrichment.csv", row.names = FALSE)
#
# # KEGG enrichment
# kegg_results <- enrichKEGG(
#   gene = entrez_ids,
#   organism = "mmu",  # "mmu" for mouse, "hsa" for human
#   pvalueCutoff = 0.05
# )
# kegg_results_df <- as.data.frame(kegg_results)
# write.csv(kegg_results_df, "KEGG_Enrichment.csv", row.names = FALSE)
#
# # Optional quick plots:
# # barplot(go_results, showCategory = 10)
# # dotplot(kegg_results, showCategory = 10)

##############################################################################
# DONE. ADAPT THIS PIPELINE TO YOUR NEEDS.
##############################################################################
# - Replace sample paths and IDs with your own.
# - Adjust QC thresholds (UMI counts, mitochondrial content).
# - Modify DE thresholds (log2FC, p-values).
# - Tweak resolution for clustering to get more or fewer clusters.
# - Update references (e.g., sc-type tissue, or use a local marker set).
# - Use appropriate organism DB for GO/KEGG (e.g., org.Hs.eg.db for human).
#
# Happy analyzing!
##############################################################################
