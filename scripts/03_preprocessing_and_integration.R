################################################################################
# 03_preprocessing_and_integration.R
#
# Purpose:
#   RNA preprocessing and batch integration using Harmony.
#
# Source:
#   Refactored from Report_from_data_analysisV4.Rmd
#
# IMPORTANT:
#   - Parameters and logic are identical to the original Rmd
#   - Only code style and structure are adjusted
################################################################################

source("scripts/01_setup_environment.R")

# ------------------------------------------------------------------------------
# Load input object
# ------------------------------------------------------------------------------

input_file <- file.path("data", "processed", "GSE181064_seurat_raw.rds")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

GSE181064.seurat <- readRDS(input_file)

# ------------------------------------------------------------------------------
# Set assay
# ------------------------------------------------------------------------------

DefaultAssay(GSE181064.seurat) <- "RNA"

# ------------------------------------------------------------------------------
# filter out Peritumoral compartment and split the dataset
# ------------------------------------------------------------------------------
GSE181064.seurat <- subset(GSE181064.seurat, tissue %in% c("pbmc",  "tumor"))
GSE181064.seurat <- split(GSE181064.seurat, GSE181064.seurat$orig.ident)

# ------------------------------------------------------------------------------
# Normalization 
# ------------------------------------------------------------------------------

GSE181064.seurat <- NormalizeData(
  GSE181064.seurat,
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# Variable feature selection 
# ------------------------------------------------------------------------------

GSE181064.seurat <- FindVariableFeatures(
  GSE181064.seurat,
  verbose = FALSE
)

features <- VariableFeatures(
  object    = GSE181064.seurat,
  nfeatures = 3000
)

# ------------------------------------------------------------------------------
# Remove TCR genes from feature selection
# ------------------------------------------------------------------------------

features <- features[
  !stringr::str_detect(features, "^TR[ABDG][VJDC]")
]

VariableFeatures(GSE181064.seurat) <- features

# ------------------------------------------------------------------------------
# Scaling 
# ------------------------------------------------------------------------------

GSE181064.seurat <- ScaleData(
  GSE181064.seurat,
  features = features,
  verbose  = FALSE
)

# ------------------------------------------------------------------------------
# PCA 
# ------------------------------------------------------------------------------

GSE181064.seurat <- RunPCA(
  GSE181064.seurat,
  features = features,
  verbose  = FALSE
)

# ------------------------------------------------------------------------------
# Harmony integration 
# ------------------------------------------------------------------------------

GSE181064.seurat <- IntegrateLayers(
  object          = GSE181064.seurat,
  method          = HarmonyIntegration,
  orig.reduction  = "pca",
  new.reduction   = "harmony",
  verbose         = FALSE
)

# ------------------------------------------------------------------------------
# Save output
# ------------------------------------------------------------------------------

output_file <- file.path(
  "data", "processed", "GSE181064_seurat_preprocessed_integrated.rds"
)

saveRDS(GSE181064.seurat, file = output_file)

message("Saved preprocessed and Harmony-integrated object to: ", output_file)
