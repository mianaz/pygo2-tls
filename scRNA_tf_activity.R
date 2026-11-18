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
net <- decoupleR::get_collectri(organism = 'mouse', 
                                split_complexes = FALSE)

# Extract the normalized counts
mat <- mPCa[["SCT"]]$data %>% as.matrix()

# Run ulm
acts <- decoupleR::run_ulm(mat = mat, 
                    net = net, 
                    .source = 'source', 
                    .target = 'target',
                    .mor='mor', 
                    minsize = 5)

write.csv(acts, paste0(table.path, "scRNA_tf_activities_ulm.csv"), row.names = FALSE)
# Extract mlm and store it in pathwaysmlm in data
mPCa[['tfsulm']] <- acts %>%
  tidyr::pivot_wider(id_cols = 'source', 
                     names_from = 'condition',
                     values_from = 'score') %>%
  tibble::column_to_rownames(var = 'source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(mPCa) <- "tfsulm"

# Scale the data
mPCa <- ScaleData(mPCa)
mPCa[["tfsulm"]]$data <- mPCa[["tfsulm"]]$scale.data

p1 <- Seurat::DimPlot(mPCa, 
                      reduction = "tsne", 
                      label = TRUE, 
                      pt.size = 0.5) + 
  Seurat::NoLegend() + 
  ggplot2::ggtitle('Cell types')

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p2 <- FeaturePlot(mPCa, features = c("Nfkb1")) + 
  ggplot2::scale_colour_gradient2(low = colors[1], mid = 'white', high = colors[2]) +
  ggplot2::ggtitle('NF-kB activity')

p <- p1 | p2
p

# Extract activities from object as a long dataframe
mPCa$celltype_group <- paste0(mPCa$celltype_final, "_", mPCa$group)
n_tfs <- 100

# Extract activities from object as a long dataframe
df <- t(as.matrix(mPCa@assays$tfsulm$data)) %>%
  as.data.frame() %>%
  dplyr::mutate(cluster = mPCa$celltype_group) %>%
  tidyr::pivot_longer(cols = -cluster, 
                      names_to = "source", 
                      values_to = "score") %>%
  dplyr::group_by(cluster, source) %>%
  dplyr::summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(mean)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  dplyr::filter(source %in% tfs) %>%
  tidyr::pivot_wider(id_cols = 'cluster', 
                     names_from = 'source',
                     values_from = 'mean') %>%
  tibble::column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 20,
                   treeheight_col = 20,
                   filename = paste0(plot.path, "scRNA_tf_activities_ulm_heatmap.png"),
                   width = 25,
                   height = 10,
                   dpi=300) 


# by group only
# Extract activities from object as a long dataframe
n_tfs <- 25
df <- t(as.matrix(mPCa@assays$tfsulm$data)) %>%
  as.data.frame() %>%
  dplyr::mutate(cluster = mPCa$group) %>%
  tidyr::pivot_longer(cols = -cluster, 
                      names_to = "source", 
                      values_to = "score") %>%
  dplyr::group_by(cluster, source) %>%
  dplyr::summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(mean)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  dplyr::filter(source %in% tfs) %>%
  tidyr::pivot_wider(id_cols = 'cluster', 
                     names_from = 'source',
                     values_from = 'mean') %>%
  tibble::column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 20,
                   treeheight_col = 20,
                   filename = paste0(plot.path, "scRNA_tf_activities_ulm_heatmap_group.png"),
                   width = 8,
                   height = 3,
                   dpi=300) 

