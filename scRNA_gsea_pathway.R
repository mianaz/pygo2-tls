cat("\014")
graphics.off()
rm(list=ls())

library(Seurat)
library(dplyr)
library(ggplot2)
library(decoupleR)
library(patchwork)
library(RColorBrewer)
library(tibble)
library(tidyr)
library(pheatmap)

this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)
set.seed(2025)

# Specify data and output path
data.path="scRNA_raw/"
output.path="output/"
table.path="output/Table/"
plot.path="output/Plots/"

# Load data
mPCa <- readRDS(paste0(output.path,"mPCa_annotated.rds"))

# get pathway-target interactions from PROGENy
net <- decoupleR::get_progeny(organism = 'mouse', top=1000)
# Extract the normalized counts
mat <- mPCa[["SCT"]]$data %>% as.matrix()

# Run mlm
acts <- decoupleR::run_mlm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor = 'weight', 
                           minsize = 5)

write.csv(acts, paste0(table.path, "scRNA_pathway_activities_mlm.csv"), row.names = FALSE)
# Extract mlm and store it in pathwaysmlm in data
mPCa[['pathwaysmlm']] <- acts %>%
  tidyr::pivot_wider(id_cols = 'source', 
                     names_from = 'condition',
                     values_from = 'score') %>%
  tibble::column_to_rownames(var = 'source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(mPCa) <- "pathwaysmlm"

# Scale the data
mPCa <- ScaleData(mPCa)
mPCa[["pathwaysmlm"]]$data <- mPCa[["pathwaysmlm"]]$scale.data

p1 <- Seurat::DimPlot(mPCa, 
                      reduction = "tsne", 
                      label = TRUE, 
                      pt.size = 0.5) + 
  Seurat::NoLegend() + 
  ggplot2::ggtitle('Cell types')

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p2 <- FeaturePlot(mPCa, features = c("Androgen")) + 
  ggplot2::scale_colour_gradient2(low = colors[1], mid = 'white', high = colors[2]) +
  ggplot2::ggtitle('Androgen activity')

p <- p1 | p2
p

# Extract activities from object as a long dataframe
mPCa$celltype_group <- paste0(mPCa$celltype_final, "_", mPCa$group)
df <- t(as.matrix(mPCa@assays$pathwaysmlm$data)) %>%
  as.data.frame() %>%
  dplyr::mutate(cluster = mPCa$celltype_group) %>%
  tidyr::pivot_longer(cols = -cluster, 
                      names_to = "source", 
                      values_to = "score") %>%
  dplyr::group_by(cluster, source) %>%
  dplyr::summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  tidyr::pivot_wider(id_cols = 'cluster', 
                     names_from = 'source',
                     values_from = 'mean') %>%
  tibble::column_to_rownames(var = 'cluster') %>%
  as.matrix()

# Color scale
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-1.25, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 1.25, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 20,
                   cellheight = 20,
                   treeheight_row = 20,
                   treeheight_col = 20,
                   filename = paste0(plot.path, "scRNA_pathway_activity_heatmap.png"),
                   width = 8,
                   height = 12,
                   dpi="300") 


saveRDS(mPCa, paste0(output.path,"mPCa_annotated_pathwaysmlm.rds"))
