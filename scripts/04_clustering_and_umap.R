################################################################################
# 04_clustering_and_umap.R
#
# Purpose:
#   Graph construction, clustering, and UMAP embedding.
#
# Source:
#   Refactored from Report_from_data_analysisV4.Rmd
#
# IMPORTANT:
#   - No parameters are invented
#   - Seurat defaults are used where the Rmd did not specify values
#   - Order of operations matches the Rmd logic
################################################################################

source("scripts/01_setup_environment.R")

# ------------------------------------------------------------------------------
# Load input object
# ------------------------------------------------------------------------------

input_file <- file.path(
  "data", "processed",
  "GSE181064_seurat_preprocessed_integrated.rds"
)

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nRun 03_preprocessing_and_integration.R first."
  )
}

GSE181064.seurat <- readRDS(input_file)

# ------------------------------------------------------------------------------
# Use Harmony reduction for downstream steps 
# ------------------------------------------------------------------------------

reduction_to_use <- "harmony"

# ------------------------------------------------------------------------------
# Neighborhood graph 
# ------------------------------------------------------------------------------

GSE181064.seurat <- FindNeighbors(
  object    = GSE181064.seurat,
  reduction = reduction_to_use,
  dims = 1:10
)

# ------------------------------------------------------------------------------
# Clustering (defaults, as in Rmd)
# ------------------------------------------------------------------------------

GSE181064.seurat <- FindClusters(
  object = GSE181064.seurat, 
  resolution = 1.2
)

# ------------------------------------------------------------------------------
# UMAP 
# ------------------------------------------------------------------------------

GSE181064.seurat <- RunUMAP(
  object    = GSE181064.seurat,
  reduction = reduction_to_use,
  dims = 1:10,
  n.neighbors = 50, 
  min.dist = .4
)

# ------------------------------------------------------------------------------
# Minimal validation
# ------------------------------------------------------------------------------

stopifnot(
  "umap" %in% names(GSE181064.seurat@reductions),
  "seurat_clusters" %in% colnames(GSE181064.seurat@meta.data)
)

# ------------------------------------------------------------------------------
# Save output
# ------------------------------------------------------------------------------

output_file <- file.path(
  "data", "processed",
  "GSE181064_seurat_clustered_umap.rds"
)

saveRDS(GSE181064.seurat, file = output_file)

message("Saved clustered Seurat object with UMAP to: ", output_file)

# ------------------------------------------------------------------------------
# Save UMAP plots
# ------------------------------------------------------------------------------
pdf(file = file.path(paths$umap, "post_clustering_UMAPs.pdf"))
UMAPPlot(GSE181064.seurat, label = T) + NoLegend() # general clusters
UMAPPlot(GSE181064.seurat, label = T, group.by = "tissue") + NoLegend() # sepparation according to cell origin
UMAPPlot(GSE181064.seurat, label = T, group.by = "orig.ident") + NoLegend() # showing that there are no batch effects between patients
FeaturePlot(GSE181064.seurat, features = "LILRB1", order = T) # ILT2 expression
dev.off()
