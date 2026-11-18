library(tidyverse)
library(cowplot)
library(ggrepel)
library(data.table)
library(tidyplots)
library(ggh4x)

##### Fig 2A TCGA pan cancer TLS signature #####
rank_data <- readRDS("./output/tcga_tls_rank_data.rds")
# Filter out PYGO2 high/low for non-target signatures
rank_data <- rank_data %>%
  filter(!(tcga_abbrev %in% c("PRAD_PYGO2_high", "PRAD_PYGO2_low") & 
             !signature_clean %in% c("TLS Imprint", "TLS 12-Chemokine")))

# Recalculate ranks after filtering
rank_data <- rank_data %>%
  group_by(signature) %>%
  arrange(desc(mean_score)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Define signature-specific cancer types to label
signature_labels <- list(
  "JAVELIN Renal 101 T effector" = list(
    top = c("DLBC", "SKCM"),
    bottom = c("GBM", "LGG")
  ),
  "IMmotion150 T effector" = list(
    top = c("DLBC", "SKCM"),
    bottom = c("GBM", "LGG")
  ),
  "TLS 12-Chemokine" = list(
    top = c("DLBC", "LUAD"),
    bottom = c("GBM", "LGG")
  ),
  "TLS Imprint" = list(
    top = c("DLBC", "LUAD"),
    bottom = c("GBM", "LGG")
  )
)

# Create label data with signature-specific cancer types
label_data <- rank_data %>%
  mutate(
    should_label = FALSE,
    label_color = case_when(
      tcga_abbrev == "PRAD" ~ "green3",
      tcga_abbrev == "PRAD_PYGO2_high" ~ "red",
      tcga_abbrev == "PRAD_PYGO2_low" ~ "blue",
      TRUE ~ "black"
    ),
    display_label = case_when(
      tcga_abbrev == "PRAD_PYGO2_high" ~ "PRAD (PYGO2_hi)",
      tcga_abbrev == "PRAD_PYGO2_low" ~ "PRAD (PYGO2_low)",
      TRUE ~ tcga_abbrev
    )
  )

# Mark cancer types that should be labeled for each signature
for (sig in names(signature_labels)) {
  top_types <- signature_labels[[sig]]$top
  bottom_types <- signature_labels[[sig]]$bottom

  label_data <- label_data %>%
    mutate(
      should_label = case_when(
        signature_clean == sig & tcga_abbrev %in% c(top_types, bottom_types) ~ TRUE,
        grepl("PRAD", tcga_abbrev) ~ TRUE,  # Always label PRAD samples
        TRUE ~ should_label
      )
    )
}

# Filter to only labeled samples
label_data <- label_data %>%
  filter(should_label) %>%
  select(-should_label)

# Add highlighting flags for signature-specific cancer types
rank_data <- rank_data %>%
  mutate(
    is_highlighted = FALSE,
    highlight_type = NA_character_
  )

# Mark cancer types for highlighting based on signature-specific labels
for (sig in names(signature_labels)) {
  top_types <- signature_labels[[sig]]$top
  bottom_types <- signature_labels[[sig]]$bottom

  rank_data <- rank_data %>%
    mutate(
      is_highlighted = case_when(
        signature_clean == sig & tcga_abbrev %in% top_types ~ TRUE,
        signature_clean == sig & tcga_abbrev %in% bottom_types ~ TRUE,
        TRUE ~ is_highlighted
      ),
      highlight_type = case_when(
        signature_clean == sig & tcga_abbrev %in% top_types ~ "top",
        signature_clean == sig & tcga_abbrev %in% bottom_types ~ "bottom",
        TRUE ~ highlight_type
      )
    )
}

# Create the faceted plot
p_faceted <- ggplot(rank_data, aes(x = mean_score, y = rank)) +
  geom_point(
    aes(
      fill = case_when(
        tcga_abbrev == "PRAD" ~ "green3",
        tcga_abbrev == "PRAD_PYGO2_high" ~ "red",
        tcga_abbrev == "PRAD_PYGO2_low" ~ "blue",
        highlight_type == "top" ~ "darkred",
        highlight_type == "bottom" ~ "navy",
        TRUE ~ "grey70"
      ),
      alpha = case_when(
        grepl("PRAD", tcga_abbrev) ~ 0.9,
        is_highlighted ~ 0.9,
        TRUE ~ 0.4
      )
    ),
    shape = 21, color = "black", size = 2.8
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = display_label, color = label_color),
    size = 4, 
    fontface = "bold"
  ) +
  scale_fill_identity() +
  scale_alpha_identity() +
  scale_color_identity() +
  scale_y_reverse() +
  facet_wrap(~ signature_clean, scales = "free_x", ncol = 2) +
  theme_cowplot() +
  theme(
    strip.background = element_rect(fill = "white", color = "grey85"),
    strip.text = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14)
  ) +
  labs(
    title = "",
    x = "Mean ssGSEA score",
    y = "Rank of TCGA cancer types"
  )
print(p_faceted)
ggsave("./output/Plots/Fig2A_TCGA_TLS_signature_rank.pdf", p_faceted, width = 7, height = 6)
ggsave("./output/Plots/Fig2A_TCGA_TLS_signature_rank.png", p_faceted, width = 7, height = 6, dpi = 300)


##### Fig 2B PCTA, SU2C PYGO2 vs TLS Imprint Violin #####
plot_data <- readRDS("pcta_tls_signature_data.rds")

# Create violin plot with p-value instead of asterisks
p<- plot_data |>
  tidyplot(x = group, y = ssGSEA_score, color = group) |>
  add_violin() |>
  add_data_points_jitter(size = 1.5, alpha=0.05) |>
  add_mean_dash(width = 0.2, linewidth = 1.5) |>
  add_test_pvalue(method = "wilcox.test", hide_info = TRUE, hide.ns = FALSE, bracket.nudge.y = 0.3, label.size=4) |>
  adjust_title("PCTA (n=1073)") |>
  adjust_y_axis_title("TLS Imprint score") |>
  adjust_x_axis_title("") |>
  adjust_x_axis(rotate_labels = TRUE) |>
  adjust_font(fontsize = 14) |>
  adjust_legend_position(position = "none") |>
  adjust_size(width = 4, height = 4, unit = "cm")

p
ggsave("./output/Plots/Fig2B_PCTA_TLS_violin.pdf", p, width = 4, height = 4)
ggsave("./output/Plots/Fig2B_PCTA_TLS_violin.png", p, width = 4, height = 4, dpi = 300)

# su2c
plot_data <- readRDS("pcta_tls_signature_data.rds")

# Create violin plot with p-value instead of asterisks
p<- plot_data |>
  tidyplot(x = group, y = ssGSEA_score, color = group) |>
  add_violin() |>
  add_data_points_jitter(size = 1.5, alpha=0.25) |>
  add_mean_dash(width = 0.2, linewidth = 1.5) |>
  add_test_pvalue(method = "wilcox.test", hide_info = TRUE, hide.ns = FALSE, bracket.nudge.y = 0.3, label.size=4) |>
  adjust_title("SU2C (n=134)") |>
  adjust_y_axis_title("TLS Imprint score") |>
  adjust_x_axis_title("") |>
  adjust_x_axis(rotate_labels = TRUE) |>
  adjust_font(fontsize = 14) |>
  adjust_legend_position(position = "none") |>
  adjust_size(width = 4, height = 4, unit = "cm")

p
ggsave("./output/Plots/Fig2B_SU2C_TLS_violin.pdf", p, width = 4, height = 4)
ggsave("./output/Plots/Fig2B_SU2C_TLS_violin.png", p, width = 4, height = 4, dpi = 300)

##### Fig 2C GSEA Dotplot #####
gsea_combined <- fread("./output/Table/GSEA_all_results.csv")

gsea_combined <- gsea_combined %>%
  mutate(
    source = factor(source, levels = c("GOBP","GOCC","GOMF","C7")),
    condition = factor(condition)
  )
selected_pathways <- c("GOMF_MHC_PROTEIN_COMPLEX_BINDING",
                       "GOCC_IMMUNOGLOBULIN_COMPLEX",
                       "GOMF_ANTIGEN_BINDING",
                       "GOMF_IMMUNE_RECEPTOR_ACTIVITY",
                       "GOBP_B_CELL_MEDIATED_IMMUNITY",
                       "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY",
                       "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN",
                       "GOBP_POSITIVE_REGULATION_OF_RESPONSE_TO_CYTOKINE_STIMULUS")

pretty_pathway <- function(x) {
  x %>%
    str_remove("^GO[:_](BP|CC|MF)[:_]") %>%
    str_remove("^(GOBP|GOCC|GOMF)[:_]") %>%
    str_replace_all("[_]+", " ") %>%
    str_squish() %>%
    str_to_title() %>%
    str_replace_all(c(
      "\\bMhc\\b" = "MHC", "\\bHla\\b" = "HLA", "\\bIfn\\b" = "IFN",
      "\\bTcr\\b" = "TCR", "\\bNk\\b" = "NK", "\\bCd(\\d+)\\b" = "CD\\1",
      "\\bDna\\b" = "DNA", "\\bRna\\b" = "RNA", "\\bAtp\\b" = "ATP",
      "\\bGtp\\b" = "GTP", "\\bMhc I\\b" = "MHC I", "\\bMhc Ii\\b"= "MHC II"
    ))
}

plot_df <- gsea_combined %>%
  filter(pathway %in% selected_pathways) %>%
  mutate(
    pathway_clean = pretty_pathway(pathway),
    pathway_clean = str_wrap(pathway_clean, width = 40)
  )

pathway_order <- plot_df %>%
  group_by(source, pathway_clean) %>%
  summarize(order_val = median(NES, na.rm = TRUE), .groups = "drop") %>%
  arrange(source, desc(order_val)) %>%
  pull(pathway_clean)

plot_df <- plot_df %>%
  mutate(
    pathway_clean = factor(pathway_clean, levels = unique(pathway_order)),
    source = fct_relevel(source, c("GOBP","GOCC","GOMF"))
  )

nes_lim <- max(2, ceiling(max(abs(plot_df$NES), na.rm = TRUE)))

strip_cols <- c(
  "GOBP" = "#E74C3C",
  "GOCC" = "#3498DB",
  "GOMF" = "#2ECC71"
)

strip_theme <- strip_themed(
  background_y = elem_list_rect(fill = strip_cols[levels(plot_df$source)]),
  text_y       = elem_list_text(color = "white", face = "bold", angle = 90)
)

p <- ggplot(plot_df, aes(x = condition, y = pathway_clean)) +
  geom_point(aes(size = -log10(padj), color = NES)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, limits = c(-nes_lim, nes_lim), name = "NES") +
  scale_size_continuous(range = c(2, 8), name = expression(-log[10]~P[adj])) +
  facet_grid2(source ~ ., scales = "free_y", space = "free_y", strip = strip_theme) +
  theme_cowplot() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    strip.text.y.right = element_text(angle = 0, hjust = 0.5, size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

print(p)
ggsave("./output/Plots/Fig2B_gsea_dotplot.pdf", p, width = 7, height = 6)
ggsave("./output/Plots/Fig2B_gsea_dotplot.png", p, width = 7, height = 6, dpi = 300)
