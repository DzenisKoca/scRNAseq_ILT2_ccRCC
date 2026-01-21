################################################################################
# 05_cluster_annotation.R
#
# Purpose:
#   Inspect marker expression and assign final cluster names.
#
# Source:
#   Refactored from Report_from_data_analysisV4.Rmd
#   Sections:
#     - 2.2 Cluster naming
#     - 2.3 Assigning cluster names
#
# IMPORTANT:
#   - Marker panels and cluster-to-label mapping are taken verbatim
#   - No new biological assumptions are introduced
################################################################################

source("scripts/01_setup_environment.R")

# ------------------------------------------------------------------------------
# Load clustered object
# ------------------------------------------------------------------------------

input_file <- file.path(
  "data", "processed",
  "GSE181064_seurat_clustered_umap.rds"
)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

GSE181064.seurat <- readRDS(input_file)

# ------------------------------------------------------------------------------
# Set identities to clusters (verbatim)
# ------------------------------------------------------------------------------

Idents(GSE181064.seurat) <- "seurat_clusters"

# ------------------------------------------------------------------------------
# Cluster inspection DotPlots (verbatim marker panels)
# ------------------------------------------------------------------------------

marker_panels <- list(
  
  # Lineage
  lineage_cd8 = c(
    "CD8A", "CD8B"
  ),
  # Cluster of interest (ILT2)
  ilt2 = c(
    "LILRB1"
  ),
  # Exhausted T cells
  exhaustion = c(
    "PDCD1", "LAG3", "HAVCR2", "TOX2"
  ),
  # Effector memory / cytotoxic T cells
  cytotoxic = c(
    "GNLY", "GZMB", "KLRD1", "GZMM"
  ),
  # Effector transitional CD8+ T cells (per publication)
  effector_transitional = c(
    "HIST1H1E", "WDR74", "AC245014.3"
  ),
  # Tissue-resident memory (RM)
  residency = c(
    "CD69", "IL7R", "ANXA1"
  ),
  # RM CD4+ T cells (high expression)
  rm_cd4 = c(
    "CDKN1A", "EGR1"
  ),
  # Cycling / proliferating CD8+ T cells
  cycling = c(
    "TOP2A", "MKI67"
  ),
  # Naive Treg
  naive_treg = c(
    "SELL", "RTKN2"
  ),
  # Activated Treg
  activated_treg = c(
    "FOXP3", "TNFRSF4"
  )
)

pdf(file.path(paths$cl_anot, "DotPlot_expression_markers.pdf"))
for (panel in names(marker_panels)) {
  print(
    DotPlot(
      GSE181064.seurat,
      features = marker_panels[[panel]],
      scale    = FALSE
    ) +
      ggtitle(panel) +
      theme(
        axis.text.x.bottom = element_text(
          angle = 45, vjust = 0.9, hjust = 0.9
        )
      )
  )
}
dev.off()

# ------------------------------------------------------------------------------
# Assign final cluster names
# ------------------------------------------------------------------------------

new.cluster.ids <- c("0"="1-Exhausted CD8+ T cells (Tumor)",
                     "1"="1-Exhausted CD8+ T cells (Tumor)",
                     "2"="2-RM CD8+ T cells (Tumor)",
                     "3"="5-RM CD4+ T cells (Tumor)",
                     "4"="6-EM CD8+ ILT2- T cells (PBMC)",
                     "5"="4-Naive CD4+ T cells (PBMC)",
                     "6"="4-Naive CD4+ T cells (PBMC)",
                     "7"="7-EM CD8+ ILT2+ T cells (PBMC)",
                     "8"="1-Exhausted CD8+ T cells (Tumor)",
                     "9"="8-EM CD8+ ILT2+ T cells (Tumor)",
                     "10"="1-Exhausted CD8+ T cells (Tumor)",
                     "11"="3-Effector transitional CD8+ T cells (PBMC)",
                     "12"="1-Exhausted CD8+ T cells (Tumor)",
                     "13"="4-Naive CD4+ T cells (PBMC)",
                     "14"="9-Activated Treg CD4+ T cells (Tumor)",
                     "15"="10-Cycling CD8+ T cells (Tumor)",
                     "16"="Drop",
                     "17"="11-Naive Treg CD4+ T cells (PBMC)",
                     "18"="Drop"
)


GSE181064.seurat <- RenameIdents(
  GSE181064.seurat,
  new.cluster.ids
)

GSE181064.seurat$annotated_cluster <- Idents(GSE181064.seurat)

# ------------------------------------------------------------------------------
# Remove artefactual cluster (verbatim)
# ------------------------------------------------------------------------------

GSE181064.seurat <- subset(
  GSE181064.seurat,
  idents = "Drop",
  invert = TRUE
)

# ------------------------------------------------------------------------------
# Save annotated object
# ------------------------------------------------------------------------------

output_file <- file.path(
  "data", "processed",
  "GSE181064_seurat_annotated.rds"
)

saveRDS(GSE181064.seurat, file = output_file)

message("Saved annotated Seurat object to: ", output_file)

# ------------------------------------------------------------------------------
# Save UMAP plots
# ------------------------------------------------------------------------------
pdf(file = file.path(paths$umap, "annotated_UMAPs.pdf"), height = 8, width = 12)
UMAPPlot(GSE181064.seurat, 
         label = T,
         repel = FALSE)
dev.off()
