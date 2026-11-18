cat("\014")
graphics.off()
rm(list=ls())

# data manipulation
library(tidyverse)
library(data.table)
library(openxlsx)

# visualization
library(cowplot)
library(ggpubr)
library(ggsci)
library(patchwork)
library(ggrepel)
library(ggplotify)
library(pheatmap)

# DEG
library(DESeq2)

# Enrichment
library(clusterProfiler)
library(fgsea)
library(GSVA)
library(ensembldb)
library(biomaRt)

this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)
set.seed(2025)

# Specify data and output path
data.path="RNA_raw/"
output.path="output/"
table.path="output/Table/"
plot.path="output/Plots/"

##### DEG #####
counts <- read.csv(paste0(data.path,"raw_counts.csv"), row.names=1)
pattern <- "^((?:DKO|TKO)(?:_CSI)?)"

# by genotype+treatment
condition <- str_match(colnames(counts), pattern)[,2]
condition <- factor(condition, levels=c("DKO", "TKO", "DKO_CSI", "TKO_CSI"))
coldata <- data.frame(row.names = colnames(counts), condition = condition)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
df <- data.frame()
res <- as.data.frame(results(dds,contrast=c("condition","TKO","DKO")))
wb <- createWorkbook()
addWorksheet(wb, "TKO_VS_DKO")
writeData(wb, sheet = "TKO_VS_DKO", res, rowNames=TRUE)
res$comparison="TKO_VS_DKO"
df <- rbind(df, res)

res <- as.data.frame(results(dds,contrast=c("condition","DKO_CSI","DKO")))
addWorksheet(wb, "DKO_CSI_VS_DKO")
writeData(wb, sheet = "DKO_CSI_VS_DKO", res, rowNames=TRUE)
res$comparison="DKO_CSI_VS_DKO"
df <- rbind(df, res)

res <- as.data.frame(results(dds,contrast=c("condition","TKO_CSI","TKO")))
addWorksheet(wb, "TKO_CSI_VS_TKO")
writeData(wb, sheet = "TKO_CSI_VS_TKO", res, rowNames=TRUE)
res$comparison="TKO_CSI_VS_TKO"
df <- rbind(df, res)


res <- as.data.frame(results(dds,contrast=c("condition","TKO_CSI","DKO_CSI")))
addWorksheet(wb, "TKO_CSI_VS_DKO_CSI")
writeData(wb, sheet = "TKO_CSI_VS_DKO_CSI", res, rowNames=TRUE)
res$comparison="TKO_CSI_VS_DKO_CSI"
df <- rbind(df, res)

saveWorkbook(wb, file = "./output/Table/DESeq2_DEGs.xlsx", overwrite = TRUE)
write.xlsx(df, file = "./output/Table/DESeq2_DEGs_all.xlsx")
saveRDS(dds, paste0(output.path, "DESeq.rds"))

##### PCA, heatmap #####
dds <- estimateSizeFactors(dds)
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1) # LRT to find DE across groups
res_lrt <- results(dds_lrt)
write.xlsx(data.frame(gene=rownames(res_lrt),as.data.frame(res_lrt)), file = "./output/Table/DESeq2_DE_results_LRT.xlsx")

# rlog normalization (for visualization only)
rlog_data <- rlog(dds, blind=TRUE)
rlog_counts <- assay(rlog_data)

norm_counts <- rlog_counts
pca_result <- prcomp(t(norm_counts))
pc_scores <- as.data.frame(pca_result$x)
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

pca_data <- cbind(pc_scores, coldata)
# pca plot
cols <- setNames(get_palette("Paired", k=length(unique(pca_data$condition))), levels(pca_data$condition))
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = rownames(pca_data))) +
  geom_point(size = 3) +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf) +
  labs(x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)"),
       color = "condition",
       fill = "condition") +
  scale_color_manual(values=cols) +
  theme_cowplot() +
  theme(legend.position = "right")
p
ggsave(paste0(plot.path, "PCA_plot.png"), plot = p, width = 6, height = 4, dpi = 300)

# heatmap (DE genes)
norm_counts <- counts(dds, normalized=TRUE)
sig_genes <- subset(res_lrt, padj<0.05) %>% rownames()
norm_counts <- norm_counts[which(rownames(norm_counts) %in% sig_genes),]
scaled_log_counts <- t(scale(t(norm_counts)))
# plot heatmap
cols <- get_palette("Paired", k=9)
coldata <- data.frame(condition = factor(condition))
rownames(coldata) <- colnames(scaled_log_counts)
colcolor <- setNames(cols[1:4], levels(coldata$condition))
anno_colors <- list(condition = colcolor, cluster = rowcolor)
col_cont <- colorRampPalette(c("#2857A3", "white", "#E41A1C"))(100)
heatmap_output <- pheatmap(scaled_log_counts,
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           color = col_cont,
                           annotation_col = coldata,
                           annotation_colors = anno_colors,
                           cluster_cols = TRUE,
                           angle_col = 45,
                           fontsize_col = 8,
                           scale = "none")

ggsave(paste0(plot.path, "Heatmap_DEGs.png"), plot = heatmap_output, width = 6, height = 6, dpi = 300)

##### GSVA #####
# For consistency, use normalized counts from deseq2
exp <- rlog_counts
gene_sets <- read.csv("./external data/Taavitsainen_gsva_gene_sets.csv")
gene_sets <- list(
  AR_Activity = na.omit(gene_sets$AR_Activity),
  AR_Full_Length = na.omit(gene_sets$AR_Full_Length_ARFL),
  AR_Variant = na.omit(gene_sets$AR_Variant_ARV),
  NEPC_Upregulated = na.omit(gene_sets$NEPC_Upregulated),
  NEPC_Downregulated = na.omit(gene_sets$NEPC_Downregulated),
  MTORC1_Signaling = na.omit(gene_sets$MTORC1_Signaling),
  MYC_Targets_1=na.omit(gene_sets$MYC_Targets_V1),
  MYC_Targets_2=na.omit(gene_sets$MYC_Targets_V2),
  Persistent_Clusters=na.omit(gene_sets$Persistent_Clusters),
  Persist=na.omit(gene_sets$Persist),
  Initial_Clusters=na.omit(gene_sets$Initial_Clusters),
  Enz_Induced_Clusters=na.omit(gene_sets$Enz_Induced_Clusters),
  BRCAness=na.omit(gene_sets$BRCAness),
  PROS_Genesis=na.omit(gene_sets$PROSGenesis)
)
gene_sets <- lapply(gene_sets, function(x) x[x != ""])
# run GSVA
gsvaPar <- gsvaParam(exp, gene_sets)
gsva.es <- gsva(gsvaPar, verbose=TRUE)
# by group
gsva_df <- as.data.frame(t(gsva.es))
gsva_df$Group <- dds$condition
gsva_df$Group <- factor(gsva_df$Group, levels=c("DKO", "TKO", "DKO_CSI", "TKO_CSI"))
gsva_df <- gsva_df %>%
  group_by(Group) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("Group")
gsva_df <- t(gsva_df)
pheatmap(gsva_df,
         color = col_cont,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         cellwidth = 15,
         cellheight = 15,
         fontsize = 10,
         angle_col = 45,
         main = "GSVA",
         breaks = seq(-0.5, 0.5, length.out = 101),
         legend_breaks = seq(-0.5, 0.5, length.out = 5),  # Custom legend range
         legend_labels = sprintf("%.2f", seq(-0.5, 0.5, length.out = 5)),
         border_color = NA,
         treeheight_col = 0,
         treeheight_row = 30)
# by sample
gsva_df <- as.data.frame(gsva.es)
pheatmap(gsva_df,
         color = col_cont,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         cellwidth = 15,
         cellheight = 15,
         fontsize = 10,
         angle_col = 45,
         main = "GSVA",
         breaks = seq(-0.5, 0.5, length.out = 101),
         legend_breaks = seq(-0.5, 0.5, length.out = 5),  # Custom legend range
         legend_labels = sprintf("%.2f", seq(-0.5, 0.5, length.out = 5)),
         border_color = NA,
         treeheight_col = 0,
         treeheight_row = 30)
ggsave(paste0(plot.path, "Heatmap_GSVA.png"), width = 6, height = 6, dpi = 300)


# Save the main DESeq object for volcano plots and other analyses
saveRDS(dds, paste0(output.path, "DESeq.rds"))
saveRDS(rlog_data, paste0(output.path, "rlog_data.rds"))

sessionInfo()
