################################################################################
# 07_tcr_analysis.R
#
# Purpose:
#   Load TCR V(D)J data, integrate it with the annotated Seurat object,
#   and compute clonotype-level metadata using scRepertoire.
#
# Source:
#   Refactored from Report_from_data_analysisV4.Rmd
#
# IMPORTANT:
#   - Uses scRepertoire defaults unless explicitly specified in the Rmd
#   - No visualization is performed here
################################################################################

source("scripts/01_setup_environment.R")

# ------------------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------------------

seurat_file <- file.path(
  "data", "processed",
  "GSE181064_seurat_annotated.rds"
)

tcr_data_dir <- file.path(
  "data", "raw", "TCRseq"
)

if (!file.exists(seurat_file)) {
  stop(
    "Annotated Seurat object not found: ", seurat_file,
    "\nRun 05_cluster_annotation.R first."
  )
}


if (!dir.exists(tcr_data_dir)) {
  stop(
    "TCR data directory not found: ", tcr_data_dir,
    "\nExpected 10X V(D)J output under data/tcr/."
  )
}

# ------------------------------------------------------------------------------
# Load Seurat object
# ------------------------------------------------------------------------------

GSE181064.seurat <- readRDS(seurat_file)

GSE181064.seurat$annotated_cluster <- 
  factor(GSE181064.seurat$annotated_cluster, 
         levels = c(
           "2-RM CD8+ T cells (Tumor)",
           "3-Effector transitional CD8+ T cells (PBMC)",
           "6-EM CD8+ ILT2- T cells (PBMC)",
           "7-EM CD8+ ILT2+ T cells (PBMC)",
           "8-EM CD8+ ILT2+ T cells (Tumor)",
           "1-Exhausted CD8+ T cells (Tumor)",
           "10-Cycling CD8+ T cells (Tumor)"
         ))

GSE181064.seurat.cd8 <- subset(GSE181064.seurat, !is.na(GSE181064.seurat$annotated_cluster))
GSE181064.seurat.cd8$annotated_cluster <- as.vector(GSE181064.seurat.cd8$annotated_cluster)

scRNAseq_barcode_data <- as_tibble(
  str_split_fixed(colnames(GSE181064.seurat.cd8),
                  pattern = "_", n = Inf)[, c(1, 2, 3, 4, 5, 6)])

colnames(scRNAseq_barcode_data) <- c("patient", "sample", "sample_id", "NA", "tissue", "cell_types")


# ------------------------------------------------------------------------------
# TCRseq data import
# ------------------------------------------------------------------------------
# Basic metadata are contained within filenames. These will be extracted 
#using a following script, in a manner which allows for combination
#with scRNAseq data

TCR_seq_file_name_data <- as_tibble(
  str_split_fixed(list.files(tcr_data_dir),
                  pattern = "_", 5)[, c(3, 4)])

colnames(TCR_seq_file_name_data) <- c("sample", "sample_id")

# Adding the rest of information needed to reconstruct the cell barodes
TCR_seq_file_name_data <- TCR_seq_file_name_data %>% 
  left_join(distinct(scRNAseq_barcode_data))


# Loading TCRseq data
GSE181064.tcr <- loadContigs(input = tcr_data_dir )

GSE181064.tcr <-
  combineTCR(
    input.data = GSE181064.tcr,
    filterNonproductive = TRUE, 
    samples = paste(TCR_seq_file_name_data$patient,
                    TCR_seq_file_name_data$sample,
                    TCR_seq_file_name_data$sample_id,
                    TCR_seq_file_name_data$`NA`,
                    TCR_seq_file_name_data$tissue,
                    TCR_seq_file_name_data$cell_types,sep = "_")
    
  )

GSE181064.tcr <- addVariable(
  input.data = GSE181064.tcr,
  variable.name = "patient",
  TCR_seq_file_name_data$patient)

GSE181064.seurat.cd8 <- combineExpression(
  input.data = GSE181064.tcr,
  proportion = TRUE,
  sc.data = GSE181064.seurat.cd8,
  cloneCall = "aa",
  group.by = "patient",
  addLabel = TRUE
)

# ------------------------------------------------------------------------------
# Compute clonal expansion categories
# ------------------------------------------------------------------------------
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

UMAP_clonal_expansion <- DimPlot(GSE181064.seurat.cd8, 
        group.by = "cloneSize", 
        pt.size = 1, order = F) +
  theme(plot.title = element_blank())+
  scale_color_manual(values=rev(colorblind_vector[c(1,2,3,5,7)]))

pdf(file = file.path(paths$umap, "UMAP_clonal_expansion.pdf"), height = 8, width = 12)
UMAP_clonal_expansion 
dev.off()

# ------------------------------------------------------------------------------
# Compute clonal overlaps
# ------------------------------------------------------------------------------
clonaloverlap_jaccard_plot <- clonalOverlap(
  input.data = GSE181064.seurat.cd8,
  method = "jaccard",
  cloneCall = "aa", 
  group.by = "annotated_cluster", 
  order.by = c(
    "2-RM CD8+ T cells (Tumor)",
    "3-Effector transitional CD8+ T cells (PBMC)",
    "6-EM CD8+ ILT2- T cells (PBMC)",
    "7-EM CD8+ ILT2+ T cells (PBMC)",
    "8-EM CD8+ ILT2+ T cells (Tumor)",
    "1-Exhausted CD8+ T cells (Tumor)",
    "10-Cycling CD8+ T cells (Tumor)"
  )
) + 
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

pdf(file = "results/figures/clonaloverlap_jaccard_plot.pdf", height = 8, width = 12)
clonaloverlap_jaccard_plot 
dev.off()

# ------------------------------------------------------------------------------
# Compute Startrac diversity
# ------------------------------------------------------------------------------
GSE181064.seurat.cd8$annotated_cluster <- 
  factor(GSE181064.seurat.cd8$annotated_cluster, 
         levels = c(
           "2-RM CD8+ T cells (Tumor)",
           "3-Effector transitional CD8+ T cells (PBMC)",
           "6-EM CD8+ ILT2- T cells (PBMC)",
           "7-EM CD8+ ILT2+ T cells (PBMC)",
           "8-EM CD8+ ILT2+ T cells (Tumor)",
           "1-Exhausted CD8+ T cells (Tumor)",
           "10-Cycling CD8+ T cells (Tumor)"
         ))

startrac_div_plot <- StartracDiversity(
  sc.data = GSE181064.seurat.cd8, 
  cloneCall = "aa",
  index = c("migr", "tran", "expa"),
  type = "tissue",
  group.by = "patient") + 
  scale_x_discrete(limits = c(
    "2-RM CD8+ T cells (Tumor)",
    "3-Effector transitional CD8+ T cells (PBMC)",
    "6-EM CD8+ ILT2- T cells (PBMC)",
    "7-EM CD8+ ILT2+ T cells (PBMC)",
    "8-EM CD8+ ILT2+ T cells (Tumor)",
    "1-Exhausted CD8+ T cells (Tumor)",
    "10-Cycling CD8+ T cells (Tumor)"
  )) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))


pdf(file = "results/figures/startrac_div_plott.pdf", height = 8, width = 12)
startrac_div_plot 
dev.off()

# ------------------------------------------------------------------------------
# Aluvial plot of shared clonotypes
# ------------------------------------------------------------------------------
GSE181064.seurat.cd8.pt18 <- subset(GSE181064.seurat.cd8, orig.ident == "Pt18")

clonalCompare.pt18.cd8 <- clonalCompare(GSE181064.seurat.cd8.pt18,
              cloneCall = "aa",
              top.clones = 1000) + 
  NoLegend() + 
  scale_x_discrete(limits = c(
    "2-RM CD8+ T cells (Tumor)",
    "3-Effector transitional CD8+ T cells (PBMC)",
    "6-EM CD8+ ILT2- T cells (PBMC)",
    "7-EM CD8+ ILT2+ T cells (PBMC)",
    "8-EM CD8+ ILT2+ T cells (Tumor)",
    "1-Exhausted CD8+ T cells (Tumor)",
    "10-Cycling CD8+ T cells (Tumor)"
  )) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

pdf(file = "results/figures/clonalCompare.pt18.cd8.pdf", height = 8, width = 12)
clonalCompare.pt18.cd8 
dev.off()
# ------------------------------------------------------------------------------
# Aluvial plot of shared clonotypes in ILT2+ cells
# ------------------------------------------------------------------------------

GSE181064.seurat.cd8.pt18.ilt2pos <- subset(GSE181064.seurat.cd8.pt18, 
                                            cells = WhichCells(
                                              GSE181064.seurat.cd8.pt18, 
                                              expression = LILRB1 > 0))

Idents(GSE181064.seurat.cd8.pt18.ilt2pos) <- "tissue"

clonalCompare.pt18.ilt2.pos <- clonalCompare(GSE181064.seurat.cd8.pt18.ilt2pos,
              cloneCall = "nt",
              chain = "both",
              top.clones = 100) + 
  NoLegend() + 
  scale_x_discrete(limits = c("tumor", "pbmc" ))+
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

pdf(file = "results/figures/clonalCompare.pt18.ilt2.pos.pdf", height = 8, width = 12)
clonalCompare.pt18.ilt2.pos 
dev.off()

message("Analysis is complete.")










