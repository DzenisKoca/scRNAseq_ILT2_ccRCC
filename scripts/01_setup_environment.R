################################################################################
# 01_setup_environment.R
#
# Purpose:
#   Centralized setup of the R environment for the project.
#   This script must be sourced by all downstream scripts.
#
# Responsibilities:
#   - Load required libraries
#   - Define global paths
#   - Set reproducibility options
#
# Non-responsibilities:
#   - Data loading
#   - Analysis
#   - Plotting
################################################################################

# Reproducibility ------------------------------------------------------------

set.seed(42)

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

# Package loading ------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "Seurat",
  "patchwork",
  "scRepertoire"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages:\n",
    paste(missing_packages, collapse = "\n"),
    "\n\nInstall them before proceeding."
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(patchwork)
  library(scRepertoire)
})

# Project paths --------------------------------------------------------------

paths <- list(
  data_raw       = "data/raw",
  data_processed = "data/processed",
  scripts        = "scripts",
  functions      = "functions",
  rmd            = "rmd",
  figures        = "results/figures",
  umap           = "results/figures/UMAP",
  cl_anot        = "results/figures/cluster_annotation",
  diff_expr      = "results/figures/differential_expression",
  tables         = "results/tables",
  output         = "output"
)

# Validate paths -------------------------------------------------------------

missing_dirs <- names(paths)[!dir.exists(unlist(paths))]

if (length(missing_dirs) > 0) {
  stop(
    "The following required directories do not exist:\n",
    paste(missing_dirs, collapse = "\n"),
    "\n\nRun 00_create_project_structure.R first."
  )
}

# Final message --------------------------------------------------------------

message("Environment setup complete.")
