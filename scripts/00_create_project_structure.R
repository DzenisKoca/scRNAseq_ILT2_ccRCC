################################################################################
# 00_create_project_structure.R
#
# Purpose:
#   Create and validate the directory structure for the project repository.
#   This script is idempotent and safe to re-run.
#
# Usage:
#   source("scripts/00_create_project_structure.R")
#
# Notes:
#   - Uses base R only
#   - Does NOT load any libraries
#   - Does NOT create files, only directories
################################################################################

# Define required directories -----------------------------------------------

project_dirs <- c(
  "data/raw",
  "data/raw/scRNA_seq",
  "data/raw/TCRseq",
  "data/raw/metadata",
  "data/metadata",
  "data/processed",
  "scripts",
  "functions",
  "rmd",
  "results/figures",
  "results/figures/UMAP",
  "results/figures/cluster_annotation",
  "results/figures/differential_expression",
  "results/tables",
  "output"
)

# Create directories --------------------------------------------------------

for (dir_path in project_dirs) {
  if (!dir.exists(dir_path)) {
    dir.create(
      path = dir_path,
      recursive = TRUE,
      showWarnings = FALSE
    )
    message("Created directory: ", dir_path)
  } else {
    message("Directory already exists: ", dir_path)
  }
}

# Validation ----------------------------------------------------------------

missing_dirs <- project_dirs[!dir.exists(project_dirs)]

if (length(missing_dirs) > 0) {
  stop(
    "Project structure validation failed. Missing directories:\n",
    paste(missing_dirs, collapse = "\n")
  )
}

message("Project directory structure successfully created and validated.")
