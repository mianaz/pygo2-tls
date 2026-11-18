library(tidyverse)
library(data.table)
library(DESeq2)
library(msigdbr)
library(fgsea)
library(GSVA)
library(HGNChelper)
library(sva)

set.seed(2025)

##### Fig 1E: AR/NE signature analysis #####
ccle_tpm <- read.csv("./external data/ccle_tpm_symbol.csv", row.names=1)
log2_ccle <- log(ccle_tpm + 1)
log2_kumar <- read.csv("./external data/gse77930.csv", row.names=1)

pygo2 <- read.csv("./external data/pygo2_tpm_symbol.csv", row.names=1)
samples_keep <- c("DKO_1", "DKO_2", "TKO_1", "TKO_2")
pygo2 <- pygo2[, colnames(pygo2) %in% samples_keep]
log2_pygo2 <- log(pygo2 + 1)

datasets <- list(Pygo2 = log2_pygo2, CCLE = log2_ccle, Kumar = log2_kumar)

ar_signature <- c("KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1",
                  "PLPP1", "PMEPA1", "PART1", "ALDH1A3", "STEAP4")
ne_signature <- c("SYP", "CHGA", "CHGB", "ENO2", "CHRNB2",
                  "SCG3", "SCN3A", "PCSK1", "ELAVL4", "NKX2-1")

common_genes <- Reduce(intersect, lapply(datasets, rownames))
for(i in names(datasets)) datasets[[i]] <- as.matrix(datasets[[i]][common_genes,])

combined_data <- do.call(cbind, datasets)
batch <- factor(rep(names(datasets), sapply(datasets, ncol)))
combined_data <- combined_data[apply(combined_data, 1, var) > 0,]

corrected_data <- tryCatch({
  ComBat(combined_data, batch, mean.only = TRUE, par.prior = TRUE, prior.plots = FALSE)
}, error = function(e) combined_data)

dataset_lengths <- sapply(datasets, ncol)
start_idx <- c(1, cumsum(dataset_lengths)[-length(dataset_lengths)] + 1)
end_idx <- cumsum(dataset_lengths)

corrected_datasets <- list()
for(i in seq_along(datasets)) {
  corrected_datasets[[names(datasets)[i]]] <- corrected_data[, start_idx[i]:end_idx[i]]
}

all_results <- list()
for(dataset_name in names(corrected_datasets)) {
  expr_mat <- corrected_datasets[[dataset_name]]
  z_scores <- t(scale(t(expr_mat)))

  ar_genes_present <- intersect(ar_signature, rownames(z_scores))
  ar_score <- if(length(ar_genes_present) > 0) {
    colMeans(z_scores[ar_genes_present,, drop=FALSE], na.rm=TRUE)
  } else rep(0, ncol(expr_mat))

  ne_genes_present <- intersect(ne_signature, rownames(z_scores))
  ne_score <- if(length(ne_genes_present) > 0) {
    colMeans(z_scores[ne_genes_present,, drop=FALSE], na.rm=TRUE)
  } else rep(0, ncol(expr_mat))

  ar_score_norm <- (ar_score - min(ar_score)) / (max(ar_score) - min(ar_score))
  ne_score_norm <- (ne_score - min(ne_score)) / (max(ne_score) - min(ne_score))

  all_results[[dataset_name]] <- data.frame(
    Sample = colnames(expr_mat), Dataset = dataset_name,
    AR_score = ar_score_norm, NE_score = ne_score_norm,
    stringsAsFactors = FALSE
  )
}

plot_data <- do.call(rbind, all_results)

cell_lines <- c("LNCAPCLONEFGC", "MDAPCA2B", "VCAP", "22RV1", "NCIH660",
                "PC3", "DU145", "PRECLH")

plot_data$Type <- NA
for(i in 1:nrow(plot_data)) {
  sample <- plot_data$Sample[i]
  is_cell_line <- any(sapply(cell_lines, function(cl) grepl(cl, sample, ignore.case = TRUE)))

  if(is_cell_line) {
    plot_data$Type[i] <- "Cell line"
  } else if(grepl("^X.*CRPC", sample)) {
    plot_data$Type[i] <- "Kumar et al CRPC"
  } else if(grepl("DKO_[0-9]", sample)) {
    plot_data$Type[i] <- "DKO"
  } else if(grepl("TKO_[0-9]", sample)) {
    plot_data$Type[i] <- "TKO"
  } else {
    plot_data$Type[i] <- "Other"
  }
}

plot_data$AR_NE_status <- case_when(
  is.na(plot_data$AR_score) | is.na(plot_data$NE_score) ~ NA_character_,
  plot_data$AR_score > 0.5 & plot_data$NE_score <= 0.5 ~ "ARPC",
  plot_data$AR_score <= 0.5 & plot_data$NE_score > 0.5 ~ "NEPC",
  plot_data$AR_score > 0.5 & plot_data$NE_score > 0.5 ~ "DPPC",
  TRUE ~ "DNPC"
)

plot_data$Classification <- ifelse(plot_data$Dataset=="Kumar",
                                   plot_data$AR_NE_status, plot_data$Type)
plot_data$Sample <- gsub("\\.", "_", plot_data$Sample)
plot_data$Sample <- gsub("^X", "", plot_data$Sample)

saveRDS(plot_data, "./output/ar_ne_signature_data.rds")

##### Fig 2A: TCGA TLS signature analysis #####
exp <- fread("./external data/TCGA_rsem_tpm.txt")
meta <- fread("./external data/TCGA_metadata.txt")

exp %>% column_to_rownames("V1") -> exp

# Load gene sets
gene_sets <- read.csv("./external data/hugaboom_TLS_signature.csv")
gene_sets <- list(
  Meylan_TLS = na.omit(gene_sets$Meylan_TLS),
  Xu_TLS = na.omit(gene_sets$Xu_TLS),
  Immotion_150_Teff = na.omit(gene_sets$Immotion_150_Teff),
  JAVELIN_RENAL_101_Teff = na.omit(gene_sets$JAVELIN_RENAL_101_Teff)
)

gene_sets <- lapply(gene_sets, function(x) {
  check <- checkGeneSymbols(x)
  check$Suggested.Symbol <- sapply(strsplit(check$Suggested.Symbol, " /// "), `[`, 1)
  x <- check$Suggested.Symbol[!is.na(check$Suggested.Symbol)]
  return(x)
})

gene_sets_present <- lapply(gene_sets, function(genes) {
  genes[genes %in% rownames(exp)]
})
gene_sets_present <- gene_sets_present[sapply(gene_sets_present, length) > 0]

# ssGSEA scores
param <- ssgseaParam(expr = as.matrix(exp), geneSets = gene_sets_present)
gsva_results <- gsva(param)

score_df <- data.frame(sample = gsub("\\.", "-", colnames(exp)))
for(gene_set in rownames(gsva_results)) {
  score_df[[gene_set]] <- gsva_results[gene_set, ]
}

# cancer type
score_df$cancer_type <- meta$`primary disease or tissue`[match(score_df$sample, meta$sample)]
score_df <- score_df[!is.na(score_df$cancer_type), ]

# TCGA abbreviation mapping
tcga_abbrev <- c(
  "Diffuse Large B-Cell Lymphoma" = "DLBC",
  "Lung Adenocarcinoma" = "LUAD",
  "Lung Squamous Cell Carcinoma" = "LUSC",
  "Thymoma" = "THYM",
  "Head & Neck Squamous Cell Carcinoma" = "HNSC",
  "Testicular Germ Cell Tumor" = "TGCT",
  "Stomach Adenocarcinoma" = "STAD",
  "Pancreatic Adenocarcinoma" = "PAAD",
  "Cervical & Endocervical Cancer" = "CESC",
  "Mesothelioma" = "MESO",
  "Bladder Urothelial Carcinoma" = "BLCA",
  "Skin Cutaneous Melanoma" = "SKCM",
  "Esophageal Carcinoma" = "ESCA",
  "Breast Invasive Carcinoma" = "BRCA",
  "Cholangiocarcinoma" = "CHOL",
  "Prostate Adenocarcinoma" = "PRAD",
  "Rectum Adenocarcinoma" = "READ",
  "Colon Adenocarcinoma" = "COAD",
  "Liver Hepatocellular Carcinoma" = "LIHC",
  "Thyroid Carcinoma" = "THCA",
  "Uterine Corpus Endometrioid Carcinoma" = "UCEC",
  "Kidney Clear Cell Carcinoma" = "KIRC",
  "Ovarian Serous Cystadenocarcinoma" = "OV",
  "Sarcoma" = "SARC",
  "Uterine Carcinosarcoma" = "UCS",
  "Kidney Chromophobe" = "KICH",
  "Kidney Papillary Cell Carcinoma" = "KIRP",
  "Acute Myeloid Leukemia" = "LAML",
  "Pheochromocytoma & Paraganglioma" = "PCPG",
  "Adrenocortical Cancer" = "ACC",
  "Brain Lower Grade Glioma" = "LGG",
  "Glioblastoma Multiforme" = "GBM",
  "Uveal Melanoma" = "UVM"
)
score_df$tcga_abbrev <- tcga_abbrev[score_df$cancer_type]
score_df <- score_df[!is.na(score_df$tcga_abbrev), ]

# mean scores
mean_scores <- score_df %>%
  group_by(tcga_abbrev, cancer_type) %>%
  summarise(across(all_of(names(gene_sets_present)), mean, na.rm = TRUE))

# Stratify TCGA-PRAD by PYGO2 expression (top 25% vs bottom 25%)
prad_samples <- score_df$sample[score_df$tcga_abbrev == "PRAD"]

pygo2_expr <- exp["PYGO2", prad_samples]
pygo2_vec <- as.numeric(pygo2_expr)
names(pygo2_vec) <- prad_samples
pygo2_vec <- pygo2_vec[!is.na(pygo2_vec)]

q25 <- quantile(pygo2_vec, 0.25)
q75 <- quantile(pygo2_vec, 0.75)

pygo2_high_samples <- names(pygo2_vec)[pygo2_vec >= q75]
pygo2_low_samples <- names(pygo2_vec)[pygo2_vec <= q25]

prad_pygo2_high <- score_df %>%
  filter(sample %in% pygo2_high_samples) %>%
  summarise(across(all_of(names(gene_sets_present)), mean, na.rm = TRUE)) %>%
  mutate(tcga_abbrev = "PRAD_PYGO2_high", cancer_type = "Prostate Adenocarcinoma (PYGO2 high)")

prad_pygo2_low <- score_df %>%
  filter(sample %in% pygo2_low_samples) %>%
  summarise(across(all_of(names(gene_sets_present)), mean, na.rm = TRUE)) %>%
  mutate(tcga_abbrev = "PRAD_PYGO2_low", cancer_type = "Prostate Adenocarcinoma (PYGO2 low)")

mean_scores_with_prad_strat <- bind_rows(mean_scores, prad_pygo2_high, prad_pygo2_low)

rank_data <- list()
for(gs in names(gene_sets_present)) {
  temp <- mean_scores_with_prad_strat[, c("tcga_abbrev", "cancer_type", gs)]
  names(temp)[3] <- "mean_score"
  temp <- temp[order(temp$mean_score, decreasing = TRUE), ]
  temp$rank <- 1:nrow(temp)
  temp$signature <- gs
  temp$is_PRAD <- grepl("PRAD", temp$tcga_abbrev)
  rank_data[[gs]] <- temp
}
rank_data <- do.call(rbind, rank_data)

rank_data <- rank_data %>%
  filter(!is.na(tcga_abbrev) & trimws(tcga_abbrev) != "")
rank_data$signature_clean <- gsub("_", " ", rank_data$signature)

sig_type <- c(
  Meylan_TLS = "TLS",
  Xu_TLS = "TLS",
  Immotion_150_Teff = "T effector",
  JAVELIN_RENAL_101_Teff = "T effector"
)

pretty_names <- c(
  Meylan_TLS = "TLS Imprint",
  Xu_TLS = "TLS 12-Chemokine",
  Immotion_150_Teff = "IMmotion150 T effector",
  JAVELIN_RENAL_101_Teff = "JAVELIN Renal 101 T effector"
)

rank_data$type <- sig_type[rank_data$signature]
rank_data$signature_clean <- pretty_names[rank_data$signature]

# Save for plotting
saveRDS(rank_data, paste0(output.path, "tcga_tls_rank_data.rds"))

##### Fig 2B PYGO2 vs TLS Imprint Correlation Analysis #####
expr <- read.delim("./external data/su2c_data_mrna_seq_fpkm_polya.txt")
expr <- expr %>%
  group_by(Hugo_Symbol) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("Hugo_Symbol")
# Stratify by PYGO2 expression
pygo2_expr <- expr["PYGO2", ]
pygo2_vec <- as.numeric(pygo2_expr)
names(pygo2_vec) <- names(pygo2_expr)

# Filter NAs
pygo2_vec <- pygo2_vec[!is.na(pygo2_vec)]

q25 <- quantile(pygo2_vec, 0.25)
q75 <- quantile(pygo2_vec, 0.75)

selected_samples <- names(pygo2_vec)[pygo2_vec <= q25 | pygo2_vec >= q75]

metadata <- data.frame(
  sample = selected_samples,
  group = ifelse(pygo2_vec[selected_samples] >= q75, "PYGO2_High", "PYGO2_Low"),
  stringsAsFactors = FALSE
)

# Subset expression matrix for selected samples only
expr_subset <- as.matrix(expr[, selected_samples])

# Run ssGSEA with new GSVA API
param <- ssgseaParam(expr = expr_subset, geneSets = tls_genes)
ssgsea_scores <- gsva(param)

# Prepare plot data
plot_data <- ssgsea_scores |>
  as.data.frame() |>
  tibble::rownames_to_column("gene_set") |>
  tidyr::pivot_longer(-gene_set, names_to = "sample", values_to = "ssGSEA_score") |>
  dplyr::left_join(metadata, by = "sample")

plot_data$group <- factor(plot_data$group, levels = c("PYGO2_Low", "PYGO2_High"))
write.csv(plot_data, "./output/Table/su2c_tls_signature_data.csv")

### su2c
gene_sets <- read.csv("./external data/Hugaboom_TLS_signature.csv")
expr <- read.csv("./external data/pcta_expression_data.csv")
expr <- expr %>%
  subset(select = -c(Entrez)) %>%
  group_by(Symbol) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("Symbol")
meta <- read.csv("./external data/pcta_metadata.csv")
colnames(expr) <- meta$Expr_id
# Stratify by PYGO2 expression
pygo2_expr <- expr["PYGO2", ]
pygo2_vec <- as.numeric(pygo2_expr)
names(pygo2_vec) <- names(pygo2_expr)

# Filter NAs
pygo2_vec <- pygo2_vec[!is.na(pygo2_vec)]

q25 <- quantile(pygo2_vec, 0.25)
q75 <- quantile(pygo2_vec, 0.75)

selected_samples <- names(pygo2_vec)[pygo2_vec <= q25 | pygo2_vec >= q75]

metadata <- data.frame(
  sample = selected_samples,
  group = ifelse(pygo2_vec[selected_samples] >= q75, "PYGO2_High", "PYGO2_Low"),
  stringsAsFactors = FALSE
)

# Subset expression matrix for selected samples only
expr_subset <- as.matrix(expr[, selected_samples])

# Define TLS gene set
tls_genes <- list(TLS_signature = gene_sets$Meylan_TLS)

# Run ssGSEA with new GSVA API
param <- ssgseaParam(expr = expr_subset, geneSets = tls_genes)
ssgsea_scores <- gsva(param)

# Prepare plot data
plot_data <- ssgsea_scores |>
  as.data.frame() |>
  tibble::rownames_to_column("gene_set") |>
  tidyr::pivot_longer(-gene_set, names_to = "sample", values_to = "ssGSEA_score") |>
  dplyr::left_join(metadata, by = "sample")

plot_data$group <- factor(plot_data$group, levels = c("PYGO2_Low", "PYGO2_High"))
saveRDS(plot_data, "./output/pcta_tls_signature_data.rds")


##### Fig 2B: GSEA analysis #####

dds <- readRDS("./output/DESeq.rds")

gs_gobp <- msigdbr(db_species = "MM", species = "Mus musculus",
                   collection = "M5", subcollection = "GO:BP") %>%
  transmute(gs_name, gene_symbol, source = "GOBP")

gs_gocc <- msigdbr(db_species = "MM", species = "Mus musculus",
                   collection = "M5", subcollection = "GO:CC") %>%
  transmute(gs_name, gene_symbol, source = "GOCC")

gs_gomf <- msigdbr(db_species = "MM", species = "Mus musculus",
                   collection = "M5", subcollection = "GO:MF") %>%
  transmute(gs_name, gene_symbol, source = "GOMF")

gs_c7   <- msigdbr(db_species = "MM", species = "Mus musculus",
                   collection = "M7") %>%                       # Immunologic signatures
  transmute(gs_name, gene_symbol, source = "C7")

genesets_df <- bind_rows(gs_gobp, gs_gocc, gs_gomf, gs_c7)

genesets_list <- split(genesets_df$gene_symbol, genesets_df$gs_name)
gs_source_lut <- genesets_df %>% distinct(gs_name, source) %>%
  mutate(source = factor(source, levels = c("GOBP","GOCC","GOMF","C7")))

# Ensure DKO is reference; contrast each non-DKO vs DKO
dds$condition <- relevel(dds$condition, ref = "DKO")
conditions <- levels(dds$condition)
non_ref_conditions <- setdiff(conditions, "DKO")

all_gsea <- vector("list", length(non_ref_conditions))
names(all_gsea) <- non_ref_conditions

for (cond in non_ref_conditions) {
  cat("  Running GSEA for", cond, "vs DKO...\n")
  res <- results(dds, contrast = c("condition", cond, "DKO"))
  ranked <- na.omit(setNames(res$stat, rownames(res))) %>% sort(decreasing = TRUE)
  gsea  <- fgsea(pathways = genesets_list, stats = ranked, minSize = 15, maxSize = 500)
  if (nrow(gsea) == 0) next
  gsea$condition <- cond
  all_gsea[[cond]] <- gsea
}

gsea_combined <- bind_rows(all_gsea)
stopifnot(nrow(gsea_combined) > 0)

gsea_combined <- gsea_combined %>%
  rename(pathway = pathway, NES = NES, padj = padj) %>%
  left_join(gs_source_lut, by = c("pathway" = "gs_name")) %>%
  mutate(
    source = fct_drop(source),
    condition = factor(condition, levels = c("DKO_CSI", "TKO", "TKO_CSI") %>% intersect(non_ref_conditions)),
  ) %>%
  dplyr::filter(!is.na(source), !is.na(condition))


fwrite(gsea_combined, "./output/Table/GSEA_all_results.csv", row.names = FALSE)

## unused ## 


sessionInfo()
