library(tidyverse)
library(cowplot)
library(SCP)
library(Seurat)
library(ggh4x)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(tidyplots)
library(patchwork)

##### Fig 3B UMAP/TSNE #####
mPCa <- readRDS("./output/mPCa_annotated.rds")
CellDimPlot(mPCa, group.by = "celltypes",label=TRUE, theme_use="theme_blank",label_insitu = TRUE)
ggsave("./output/Plots/Fig3B_tsne_celltype.pdf",height=5, width=6,units = "in", device = "pdf")
ggsave("./output/Plots/Fig3B_tsne_celltype.png",height=5, width=6,units = "in", device = "png", dpi=600)

CellDimPlot(mPCa, group.by = "celltypes",theme_use="theme_blank", split.by = "group")+ plot_layout(guides = "collect")
ggsave("./output/Plots/Fig3B_tsne_celltype_split.pdf",height=5, width=12,units = "in", device = "pdf")
ggsave("./output/Plots/Fig3B_tsne_celltype_split.png",height=5, width=12,units = "in", device = "png", dpi=600)

p1 <- CellStatPlot(mPCa, group.by = "group", stat.by = "celltypes")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
p2 <- CellStatPlot(mPCa, group.by = "sample_id", stat.by = "celltypes")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
p <- p1+p2+plot_layout(guides = "collect")
ggsave("./output/Plots/Fig3B_celltype_proportion.png",height=5, width=8,units = "in", device = "png")

##### Fig 3C marker dotplot#####
cano <- data.frame(
  Epi=c("Epcam","Krt8","Krt18"),
  Fibro=c("Dcn","Sfrp4","Col3a1"),
  T_reg=c("Foxp3","Itgb8", "Neb"),
  T_gd=c("Rorc","Il23r","Trdc"),
  T_CD4=c("Cd4","Cd40lg","Cd3d"),
  T_CD8=c("Cd8a","Cd8b1","Gzmk"),
  NK=c("Itga2","Ncr1","Klrb1c"),
  Neutrophil=c("Csf3r","S100a8", "Siglecf"),
  Macro_Mono=c("Csf1r","Vcan","Chil3"),
  ILC2=c("Il5","Areg","Il1rl1"),
  DC=c("Clec9a","Xcr1","Batf3"),
  B=c("Cd19","Ms4a1","Ighm")
)

cano <- cano %>%
  pivot_longer(cols = everything(), names_to = "cluster", values_to = "gene")

p <- DotPlot(mPCa,
             features = unique(cano$gene),
             cols = c("#ffffff", "firebrick3"),
             group.by = "celltypes")
gene_cluster_map <- cano %>%
  dplyr::select(gene, cluster) %>%
  distinct()

p$data <- p$data %>%
  left_join(gene_cluster_map, by = c("features.plot" = "gene"))

p$data$cluster <- factor(p$data$cluster,
                         levels = c("Epi", "Fibro", "Macro_Mono", "DC", "NK",
                                    "ILC2", "T_CD4", "T_CD8", "T_gd",
                                    "T_reg", "B", "Neutrophil"))

p$data$id <- factor(p$data$id, levels = rev(levels(p$data$id)))
p$data <- p$data %>%
  dplyr::group_by(features.plot) %>%
  dplyr::filter(avg.exp.scaled > 0)
cols <- get_palette("Paired", k = length(unique(p$data$cluster)))

strip <- strip_themed(
  background_x = elem_list_rect(fill = cols),
  text_x = elem_list_text(color = "white")
)

p1 <- p$data %>%
  ggplot(aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  facet_wrap2(~cluster,
              scales = "free_x",
              strip = strip,
              nrow = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title = element_blank(),
        strip.text.x = element_text(size = 14),
        strip.background = element_rect(color = "white"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  scale_color_gradient(low = "#ffffff",
                       high = "firebrick3",
                       name = "avg.exp") +
  scale_size(name = "% exp")

p1
ggsave("./output/Plots/Fig3C_marker_dotplot.pdf", p1, height = 4, width = 12, units = "in", device = "pdf")
ggsave("./output/Plots/Fig3C_marker_dotplot.png", p1, height = 4, width = 12, units = "in", device = "png",dpi=600)

##### Fig 3D CD8/Treg Ratio Bar Plot #####
mPCa$sample_id_fixed <- gsub("_", "-", mPCa$sample_id)
cell_counts <- mPCa@meta.data %>%
  mutate(sample_id = sample_id_fixed) %>%
  filter(celltypes %in% c("T_CD8", "T_reg")) %>%
  group_by(sample_id, celltypes, group) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = celltypes, values_from = n, values_fill = 0) %>%
  mutate(cd8_treg_ratio = T_CD8 / (T_reg + 1))

cell_counts$sample_id <- factor(cell_counts$sample_id,
                                levels = c("DKO1", "DKO2", "TKO1", "TKO2", "TKO-CSI1", "TKO-CSI2"))
cell_counts$group <- factor(cell_counts$group, levels = c("DKO", "TKO", "TKO_CSI"))

p <- ggplot(cell_counts, aes(x = sample_id, y = cd8_treg_ratio, fill = group)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = c("DKO" = "#E41A1C", "TKO" = "#377EB8", "TKO_CSI" = "#4DAF4A")) +
  labs(x = "", y = "CD8/Treg Ratio", fill = "Group") +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

print(p)
ggsave("./output/Plots/Fig3D_CD8_Treg_ratio.pdf", p, width = 4, height = 4, units = "in", device = "pdf")
ggsave("./output/Plots/Fig3D_CD8_Treg_ratio.png", p, width = 4, height = 4, units = "in", device = "png", dpi = 600)

##### Fig 3E T cell proportions#####
t <- subset(mPCa, celltypes %in% c("T_naive", "T_reg", "T_gd", "T_CD4", "T_CD8"))
p1 <- CellStatPlot(t, group.by = "group", stat.by = "celltypes", palcolor=get_palette("Paired",k=12)[7:10])
p2 <- CellStatPlot(t, group.by = "sample_id", stat.by = "celltypes",palcolor=get_palette("Paired",k=12)[7:10])
p <- p1+p2+plot_layout(guides = "collect")
p
ggsave("./output/Plots/Fig3E_tcell_proportion_plot.png", p,
       height = 6, width =6, units = "in", device = "png", dpi=600)

##### Fig 3E Chemokine/checkpoint heatmap #####
t_cells <- subset(mPCa, celltypes %in% c("T_CD8", "T_CD4"))
DefaultAssay(t_cells) <- "RNA"

genes_plot <- c("Pdcd1", "Gzmb", "Ifng", "Tnf", "Ccr7", "Cxcl5")

expr_data <- GetAssayData(t_cells, layer = "data")[genes_plot, ] %>%
  as.matrix() %>%
  t() %>%
  as.data.frame()

expr_data$celltype <- t_cells$celltypes

expr_long <- expr_data %>%
  pivot_longer(
    cols = all_of(genes_plot),
    names_to = "gene",
    values_to = "expression"
  )

plot_list <- list()

for (celltype in c("T_CD4", "T_CD8")) {
  cells_subset <- subset(t_cells, celltypes == celltype)

  for (gene in genes_plot) {
    expr_data <- data.frame(
      expression = GetAssayData(cells_subset, slot = "data")[gene, ],
      group = cells_subset$group
    )

    expr_data$group <- as.factor(expr_data$group)

    plot_list[[paste(celltype, gene, sep = "_")]] <- expr_data %>%
      tidyplot(x = group, y = expression, color = group) %>%
      add_violin() %>%
      add_data_points_jitter(size = 0.5, alpha = 0.3) %>%
      add_mean_dash(width = 0.3, linewidth = 0.5) %>%
      add_test_pvalue(method = "wilcox_test", hide_info = TRUE, ref.group = "DKO") %>%
      adjust_y_axis(title = paste(gene, "Expression")) %>%
      adjust_x_axis(title = "", rotate_labels = TRUE) %>%
      add_title(paste(celltype, "-", gene)) %>%
      adjust_font(fontsize=14)
  }
}

wrap_plots(plot_list, ncol = 4, guides = "collect")
ggsave(
  "./output/Plots/Fig3E_chemokine_checkpoint_violin.png",
  wrap_plots(plot_list, ncol = 4, guides = "collect"),
  height = 8,
  width = 12,
  units = "in",
  device = "png",
  dpi = 600
)

genes_plot <- genes_plot[genes_plot %in% rownames(t_cells)]
genes_of_interest <- genes_plot

avg_expr_cd8 <- AverageExpression(
  subset(t_cells, celltypes == "T_CD8"),
  features = genes_of_interest, group.by = "sample_id", assays = "RNA", slot = "data"
)$RNA

avg_expr_cd4 <- AverageExpression(
  subset(t_cells, celltypes == "T_CD4"),
  features = genes_of_interest, group.by = "sample_id", assays = "RNA", slot = "data"
)$RNA

common_samples <- intersect(colnames(avg_expr_cd8), colnames(avg_expr_cd4))
avg_expr_cd8 <- avg_expr_cd8[, common_samples, drop = FALSE]
avg_expr_cd4 <- avg_expr_cd4[, common_samples, drop = FALSE]

sample_order <- colnames(avg_expr_cd8)

zscore <- function(m) t(scale(t(m)))
cd8_z <- zscore(avg_expr_cd8)
cd4_z <- zscore(avg_expr_cd4)

keep_genes <- rownames(cd8_z)[apply(cd8_z, 1, sd, na.rm = TRUE) > 0 | apply(cd4_z, 1, sd, na.rm = TRUE) > 0]
cd8_z <- cd8_z[keep_genes, , drop = FALSE]
cd4_z <- cd4_z[keep_genes, , drop = FALSE]

colnames_cd8 <- paste0(rownames(cd8_z), " (CD8)")
colnames_cd4 <- paste0(rownames(cd4_z), " (CD4)")
combined <- cbind(t(cd8_z), t(cd4_z))
colnames(combined) <- c(colnames_cd8, colnames_cd4)
combined <- combined[sample_order, , drop = FALSE]

col_block <- factor(c(rep("CD8 T cells", length(colnames_cd8)),
                      rep("CD4 T cells", length(colnames_cd4))),
                    levels = c("CD8 T cells", "CD4 T cells"))

sample_groups <- data.frame(
  sample_id = sample_order,
  group = c("DKO","DKO","TKO","TKO","TKO_CSI","TKO_CSI")[seq_along(sample_order)]
)

mPCa$sample_id_fixed <- gsub("_", "-", mPCa$sample_id)
cell_counts <- mPCa@meta.data %>%
  mutate(sample_id = sample_id_fixed) %>%
  filter(celltypes %in% c("T_CD8", "T_reg")) %>%
  group_by(sample_id, celltypes) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = celltypes, values_from = n, values_fill = 0) %>%
  mutate(cd8_treg_ratio = T_CD8 / (T_reg + 1))
cell_counts <- cell_counts[match(sample_order, cell_counts$sample_id), ]

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)
my_breaks <- c(seq(-2, 0, length.out = 50), seq(0.01, 2, length.out = 50))
col_fun <- colorRamp2(breaks = my_breaks, colors = colors.use)

row_ha <- rowAnnotation(
  Group = sample_groups$group,
  col = list(Group = c("DKO"="#E41A1C","TKO"="#377EB8","TKO_CSI"="#4DAF4A")),
  width = unit(5, "mm"),
  show_annotation_name = FALSE
)
ratio_bar <- rowAnnotation(
  "CD8/Treg\nratio" = anno_barplot(cell_counts$cd8_treg_ratio,
                                   bar_width = 0.8,
                                   gp = gpar(fill = "#7B1FA2", col = "white"),
                                   width = unit(2.5, "cm")),
  annotation_name_side = "top",
  annotation_name_rot = 0
)

col_titles <- c("CD8 T cells", "CD4 T cells")

ht_main <- Heatmap(
  as.matrix(combined),
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 14),
  column_names_gp = gpar(fontsize = 14),
  column_title_gp = gpar(fontsize = 14),
  column_split = col_block,
  column_title = col_titles,
  column_gap = unit(3, "mm"),
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 14),
    title_gp = gpar(fontsize = 14)
  )
)

ht_list <- row_ha + ht_main + ratio_bar
print(ht_list)

png("./output/Plots/Fig3E_Tcells_CD8_CD4_core_plus.png",
    width = 10, height = 3.2, units = "in", res = 300)
draw(ht_list)
dev.off()

##### Fig 3F TLS Imprint score #####
tls.data <- mPCa@meta.data[,c("TLS-imprint", "group")]
p <- tls.data |>
  tidyplot(x=group, y=`TLS-imprint`) |>
  add_violin() |>
  add_data_points(size=0.5, alpha=0.1) |>
  add_test_pvalue(method="wilcox_test", hide_info = TRUE) |>
  adjust_y_axis(title = "TLS Imprint Score") |>
  adjust_x_axis(title = "", rotate_labels = TRUE) |>
  adjust_font(fontsize=14)
p
ggsave("./output/Plots/Fig3F_TLS_imprint_score.png", p, height = 4, width = 4, units = "in", device = "png", dpi=600)

##### Fig 3G B/T cells dotplot #####
b_cells <- readRDS("./output/mPCa_b_cells.rds")
p1 <- DotPlot(b_cells, features = c("Cxcr5", "Ccr7"), cols = c("#ffffff", "firebrick3"), group.by = "group")+
  coord_flip()+
  ggtitle("B cells")+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=14),
        axis.text.y=element_text(size=14),
        axis.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 14))
p1

mhc1 <- c("H2-D1", "H2-K1", "H2-L")
mhc2 <- c("H2-Aa", "H2-Ab1", "H2-Eb1")
p1 <- DotPlot(b_cells, features = c("MHC-I", "MHC-II",mhc1, mhc2), cols = c("#ffffff", "firebrick3"), group.by = "group")+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=14),
        axis.text.y=element_text(size=14),
        axis.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
p1

t_cells_subset <- subset(mPCa, celltypes %in% c("T_CD4", "T_CD8", "T_reg", "T_gd"))
p2 <- DotPlot(t_cells_subset, features = c("Cxcr5", "Ccr7"), cols = c("#ffffff", "firebrick3"), group.by = "group")+
  coord_flip()+
  ggtitle("T cells")+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=14),
        axis.text.y=element_text(size=14),
        axis.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 14))
p2

p1+p2
ggsave("./output/Plots/Fig3G_H_BT_dotplot.png", p1+p2+plot_layout(nrow=1),
       height = 4, width =8, units = "in", device = "png", dpi=600)

##### Fig 3H Ig dotplot #####
b <- subset(mPCa, celltype_final %in% c("B"))
p <- DotPlot(b, features = c("Ighg1","Ighg3", "Igha", "Ighm", "Ighd"), cols = c("#ffffff", "firebrick3"), group.by = "sample_id",assay="RNA")+
  coord_flip()+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=14),
        axis.text.y=element_text(size=14),
        axis.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
p
ggsave("./output/Plots/Fig3H_Ig_dotplot.png", p, height = 4, width = 4, units = "in", device = "png", dpi=600)

##### Fig 3I MHC scores #####
p <- DotPlot(mPCa, features = c("MHC-I", "MHC-II", "MHC"), cols = c("#ffffff", "firebrick3"), group.by = "group", assay = "RNA") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
p
ggsave("./output/Plots/Fig3I_MHC_score.png", p, height = 4, width = 4, units = "in", device = "png", dpi = 600)

mhc2.data <- b_cells@meta.data[,c("MHC-II", "group")]
p <- mhc2.data |>
  tidyplot(x=group, y=`MHC-II`) |>
  add_violin() |>
  add_data_points(size=0.5, alpha=0.1) |>
  add_test_pvalue(method="wilcox_test", hide_info = TRUE) |>
  adjust_y_axis(title = "MHC-II Score") |>
  adjust_x_axis(title = "", rotate_labels = TRUE) |>
  adjust_font(fontsize=14)
p
ggsave("./output/Plots/Fig3J_MHC2_score.png", p, height = 4, width = 4, units = "in", device = "png", dpi=600)
