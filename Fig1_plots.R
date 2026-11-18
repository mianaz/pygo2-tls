library(DESeq2)
library(tidyplots)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(tidyplots)

##### Fig 1D volcano plot #####
dds <- readRDS("./output/DESeq.rds")
df <- results(dds, contrast=c("condition", "TKO", "DKO")) %>%
  as.data.frame() %>%
  mutate(
    neg_log10_p = -log10(pvalue),
    direction = if_else(log2FoldChange > 0, "up", "down", NA),
    candidate = abs(log2FoldChange) >= 1 & pvalue < 0.05,
    gene_id = rownames(.)
  )

genes_label <- c("Myc", "Tgfb1", "Tmprss2", "Slc45a3", "Eno2", "Ncam1", "Syp", "Pate4", "Svs3a", "Svs3b", "Svs5")

df |>
  tidyplot(x = log2FoldChange, y = neg_log10_p) |>
  add_data_points(data = filter_rows(!candidate),
                  color = "lightgrey", rasterize = FALSE) |>
  add_data_points(data = filter_rows(candidate, direction == "up"),
                  color = "#FF7777", alpha = 0.3) |>
  add_data_points(data = filter_rows(candidate, direction == "down"),
                  color = "#7DA8E6", alpha = 0.3) |>
  add_data_points(data = filter_rows(gene_id %in% genes_label, candidate, direction == "up"),
                  color = "red", size = 1.2, alpha=1) |>
  add_data_points(data = filter_rows(gene_id %in% genes_label, candidate, direction == "down"),
                  color = "blue", size = 1.2, alpha=1) |>
  add_reference_lines(x = c(-1, 1), y = -log10(0.05)) |>
  adjust_x_axis(limits=c(-10, 10)) |>
  adjust_y_axis(limits=c(0, 20)) |>
  add_data_labels_repel(data = filter_rows(gene_id %in% genes_label), label = gene_id,
                        color = "#000000", background = TRUE, fontsize = 12) |>
  adjust_x_axis_title("$log[2]~FC$") |>
  adjust_y_axis_title("$-log[10]~p-value$") |>
  add_title("TKO vs DKO") |>
  adjust_font(fontsize=14) |>
  adjust_size(width = 6, height = 6, unit="cm") |>
  tidyplots::save_plot(filename="./output/Plots/Fig1D_TKO_vs_DKO_volcano.pdf", dpi=300) |>
  tidyplots::save_plot(filename="./output/Plots/Fig1D_TKO_vs_DKO_volcano.png", dpi=300)

##### Fig 1E scatter plot #####
plot_data <- readRDS("./output/ar_ne_signature_data.rds")

plot_data$Classification <- factor(plot_data$Classification,
  levels = c("ARPC", "NEPC", "DPPC", "DNPC", "Cell line", "DKO", "TKO"))

color_values <- c(
  "ARPC" = "pink",
  "NEPC" = "#fdae61",
  "DPPC" = "#2c7bb6",
  "DNPC" = "darkgreen",
  "Cell line" = "grey30",
  "DKO" = "blue",
  "TKO" = "red"
)

shape_values <- c(
  "ARPC" = 16, "NEPC" = 16, "DPPC" = 16, "DNPC" = 16,
  "Cell line" = 0,"DKO" = 17, "TKO" = 17)

plot_data$size <- ifelse(plot_data$Dataset == "Kumar", 2, 3)
plot_data$alpha <- ifelse(plot_data$Dataset == "Kumar", 0.4, 1)

p <- ggplot(plot_data,
            aes(x = NE_score, y = AR_score,
                color = Classification, shape = Classification)) +
  geom_point(aes(alpha = alpha), size = 3) +
  scale_alpha_identity() +
  scale_color_manual(values = color_values) +
  scale_shape_manual(values = shape_values) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "NE signature", y = "AR signature") +
  theme_cowplot() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = "right"
  )

label_data <- rbind(
  subset(plot_data, Classification == "Cell line") %>%
    mutate(label = gsub("_PROSTATE", "", Sample)),
  subset(plot_data, Classification %in% c("DKO", "TKO")) %>%
    mutate(label = Classification)
)

if(nrow(label_data) > 0) {
  p <- p + geom_text_repel(
    data = label_data,
    aes(label = label),
    size = 5,
    show.legend = FALSE,
    box.padding = 0.5,
    force = 2,
    max.overlaps = Inf
  )
}

print(p)
ggsave("./output/Plots/Fig1E_Scatter_AR_NE.png", p, width = 5, height = 4, dpi = 300)
ggsave("./output/Plots/Fig1E_Scatter_AR_NE.pdf", p, width = 5, height = 4)