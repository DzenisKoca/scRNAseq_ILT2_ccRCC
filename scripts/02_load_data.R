################################################################################
# 02_load_data.R
#
# Purpose:
#   Load raw scRNA-seq data and create the initial Seurat object.
#
# Source:
#   Refactored from Report_from_data_analysisV4.Rmd (section 1.2)
#
# Notes:
#   - All biological and technical parameters are unchanged
#   - Paths are aligned to project directory structure
################################################################################

source("scripts/01_setup_environment.R")

# ------------------------------------------------------------------------------
# Parameters 
# ------------------------------------------------------------------------------

raw_data_dir <- file.path("data", "raw", "scRNA_seq")
project_name <- "ccRCC_T_Cells"
min_cells    <- 3
min_features <- 200

raw_tcr_data_dir <- file.path("data", "raw", "TCRseq")

# ------------------------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------------------------

if (!dir.exists(raw_data_dir)) {
  stop(
    "Raw data directory does not exist: ", raw_data_dir,
    "\nExpected raw 10X data under data/raw/"
  )
}

# ------------------------------------------------------------------------------
# Load raw counts
# ------------------------------------------------------------------------------

message("Reading raw scRNA-seq data from: ", raw_data_dir)

raw_counts <- Read10X(data.dir = raw_data_dir)

# ------------------------------------------------------------------------------
# Create Seurat object 
# ------------------------------------------------------------------------------

GSE181064.seurat <- CreateSeuratObject(
  counts       = raw_counts,
  project      = project_name,
  min.cells    = min_cells,
  min.features = min_features, 
)

colnames(GSE181064.seurat) <- paste0(colnames(GSE181064.seurat), "-1") 

# ------------------------------------------------------------------------------
# Minimal validation
# ------------------------------------------------------------------------------

stopifnot(
  inherits(GSE181064.seurat, "Seurat"),
  ncol(GSE181064.seurat) > 0,
  nrow(GSE181064.seurat) > 0
)

message(
  "Created Seurat object with ",
  ncol(GSE181064.seurat), " cells and ",
  nrow(GSE181064.seurat), " features."
)
# ------------------------------------------------------------------------------
# Create basic metadata
# ------------------------------------------------------------------------------
# Basic metadata are contained within cell barcodes. These will be extracted 
#using a following script

scRNAseq_barcode_data <- as_tibble(
  str_split_fixed(colnames(GSE181064.seurat),
                  pattern = "_", n = Inf)[, c(1, 2, 3, 4, 5, 6)])

colnames(scRNAseq_barcode_data) <- c("patient", "sample", "sample_id", "NA", "tissue", "cell_types")

GSE181064.seurat[["tissue"]] <- scRNAseq_barcode_data$tissue
GSE181064.seurat[["patient"]] <- scRNAseq_barcode_data$patient
GSE181064.seurat[["sample"]] <- scRNAseq_barcode_data$sample
GSE181064.seurat[["sample_id"]] <- scRNAseq_barcode_data$sample_id

# ------------------------------------------------------------------------------
# Save object
# ------------------------------------------------------------------------------

output_file <- file.path("data", "processed", "GSE181064_seurat_raw.rds")

saveRDS(GSE181064.seurat, file = output_file)

message("Saved raw Seurat object to: ", output_file)
