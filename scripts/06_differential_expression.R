################################################################################
# 06_differential_expression.R
#
# Purpose:
#   Perform differential expression analysis between selected annotated clusters
#   and visualize top markers using DotPlot.
#
# Source:
#   Refactored from Report_from_data_analysisV4.Rmd
#
# IMPORTANT:
#   - Contrasts, parameters, and logic are identical to the Rmd
#   - This script is exploratory / descriptive, not a full DE pipeline
################################################################################

source("scripts/01_setup_environment.R")

# ------------------------------------------------------------------------------
# Load annotated object
# ------------------------------------------------------------------------------

input_file <- file.path(
  "data", "processed",
  "GSE181064_seurat_annotated.rds"
)

if (!file.exists(input_file)) {
  stop(
    "Annotated Seurat object not found: ", input_file,
    "\nRun 05_cluster_annotation.R first."
  )
}

GSE181064.seurat <- readRDS(input_file)
GSE181064.seurat <- JoinLayers(GSE181064.seurat)
# ------------------------------------------------------------------------------
# Set identities to annotated clusters
# ------------------------------------------------------------------------------

Idents(GSE181064.seurat) <- "annotated_cluster"

# ------------------------------------------------------------------------------
# Differential expression
# ------------------------------------------------------------------------------


diff.expression <- FindMarkers(
  object   = GSE181064.seurat,
  random.seed = 42,
  ident.1 = "8-EM CD8+ ILT2+ T cells (Tumor)",
  ident.2 = "1-Exhausted CD8+ T cells (Tumor)",
  min.pct  = 0.05,
  logfc.threshold = 0,
  test.use = "MAST",
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# Visualization of top markers 
# ------------------------------------------------------------------------------
pdf(file.path(paths$diff_expr, "DotPlot_upregulated_ILT2+.pdf"), height = 10, width = 15)

diff.expression <- dplyr::arrange(
  diff.expression,
  avg_log2FC
)

DotPlot(
  GSE181064.seurat,
  features = rownames(diff.expression)[1:20]
) +
  theme(
    axis.text.x.bottom = element_text(
      angle = 45,
      vjust = 0.9,
      hjust = 0.9
    )
  )
dev.off()
# ------------------------------------------------------------------------------
# Visualization of bottom markers 
# ------------------------------------------------------------------------------
pdf(file.path(paths$diff_expr, "DotPlot_downregulated_ILT2+.pdf"), height = 10, width = 15)

diff.expression <- dplyr::arrange(
  diff.expression,
  desc(avg_log2FC)
)

DotPlot(
  GSE181064.seurat,
  features = rownames(diff.expression)[1:20]
) +
  theme(
    axis.text.x.bottom = element_text(
      angle = 45,
      vjust = 0.9,
      hjust = 0.9
    )
  )
dev.off()
# ------------------------------------------------------------------------------
# Volcano plot
# ------------------------------------------------------------------------------
diff.expression.volcano <- as_tibble(diff.expression, rownames = "gene")

GSE181064.seurat.mean <-
  AverageExpression(GSE181064.seurat,
                    group.by = "annotated_cluster",
                    return.seurat = TRUE, ) 

GSE181064.seurat.mean <- as_tibble(Seurat::GetAssayData(GSE181064.seurat.mean), 
                                   rownames = "gene") %>% 
  select(all_of(c("gene", "g1-Exhausted CD8+ T cells (Tumor)", "g8-EM CD8+ ILT2+ T cells (Tumor)")))



diff.expression.volcano <- diff.expression.volcano %>% 
  left_join(GSE181064.seurat.mean) %>% 
  mutate(mean_expression = `g1-Exhausted CD8+ T cells (Tumor)` +  `g8-EM CD8+ ILT2+ T cells (Tumor)`) %>% 
  mutate(mean_expression = log2(mean_expression + 1))

genes_label <- c(
  "LILRB1",
  "KLRF1",
  "ITGAM",
  "CX3CR1",
  "FCRL6",
  "NCR3",
  "KLRG1",
  "FCGR3A",
  "ADGRG1",
  "IL7R",
  "CD69",
  "CD4",
  "CD8A",
  "CD8B",
  "CD3D",
  "CD3E",
  "CD3G",
  "CXCR3",
  "CXCR4",
  "PDCD1",
  "TIGIT",
  "HAVCR2",
  "LAG3",
  "CTLA4",
  "CD27",
  "HLA-DRA",
  "HLA-DRB1",
  "HLA-DRB5"
)

diff.expression.volcano <- diff.expression.volcano %>% 
  mutate(Significance = 
           case_when(avg_log2FC > 0 & p_val_adj < 0.05 ~ "EM CD8+ILT2+ T cells (Tumor)",
                     avg_log2FC < 0 & p_val_adj < 0.05 ~ "Exhausted CD8+ T cells (Tumor)", 
                     .default = "Non-significant")) %>% 
  arrange(Significance)

volcano.plot <- diff.expression.volcano %>%
  ggplot(aes(x = mean_expression, y = avg_log2FC, colour = Significance)) +
  geom_point() +
  scale_color_manual(values =  c("EM CD8+ILT2+ T cells (Tumor)" = "skyblue", 
                                "Exhausted CD8+ T cells (Tumor)" = "red",
                                "Non-significant" = "gray"
                                )) +
  ggrepel::geom_text_repel(
    data = diff.expression.volcano %>%  filter(gene %in% genes_label),
    aes(x = mean_expression, y = avg_log2FC, label = gene),
    color = "black"
  ) +
  geom_point(
    data = diff.expression.volcano %>%  filter(gene %in% genes_label),
    aes(x = mean_expression, y = avg_log2FC),
    shape = 3,
    color = "black"
  ) +
  geom_hline(yintercept = 0) +
  ggpubr::theme_pubr()

pdf(
  file.path(paths$diff_expr, "volcano_plot.pdf"),
  height = 10,
  width = 10
)
volcano.plot
dev.off()


# ------------------------------------------------------------------------------
# Save DE table
# ------------------------------------------------------------------------------

output_file <- file.path(
  "results", "tables",
  "DE_EM_CD8_ILT2pos_vs_Exhausted_CD8_Tumor.tsv"
)

write.table(
  diff.expression,
  file      = output_file,
  sep       = "\t",
  quote     = FALSE,
  col.names = NA
)

message("Saved differential expression table to: ", output_file)
