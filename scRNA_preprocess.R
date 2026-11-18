cat("\014")
graphics.off()
rm(list=ls())

# General
library(dplyr)
library(ggplot2)
library(stringr)
library(Matrix)

# Analysis
library(Seurat)
library(scDblFinder)
library(harmony)
library(HGNChelper)
library(SingleCellExperiment)
library(decontX)
library(gprofiler2)
library(SingleR)
library(celldex)
library(scrapper)
library(data.table)

# Plotting
library(patchwork)
library(RColorBrewer)

this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)
set.seed(2025)

source("./custom_annotation.R")

# Specify data and output path
data.path="scRNA_raw/"
output.path="output/"
table.path="output/Table/"
plot.path="output/Plots/"

# Obtain list of data files
if(dir.exists(data.path)) files <- list.files(path=data.path) else print("Invalid Data Path!")
# Create output path
if(!dir.exists(output.path)) dir.create(output.path, recursive = T)
if(!dir.exists(table.path)) dir.create(table.path, recursive = T)
if(!dir.exists(plot.path)) dir.create(plot.path, recursive = T)

files <- files[grep("_filtered_feature_bc_matrix.h5", files)]
sample_names <- gsub("_sample_filtered_feature_bc_matrix.h5", "", 
                     grep("sample_filtered_feature_bc_matrix.h5", files, value=T))


genotype <- str_match(sample_names, pattern <- "^(TKO|DKO)")[,2]
treatment <- str_match(sample_names, pattern <- "(CSI)")[,2]
treatment[is.na(treatment)] <- "None"

samples <- list()
for(i in 1:length(sample_names)){
  message("Reading ", sample_names[i], "...")
  dat <- Read10X_h5(paste0(data.path, files[i]))
  dat <- dat[["Gene Expression"]]
  
  message("Checking gene symbols...")
  check <- suppressMessages(checkGeneSymbols(rownames(dat), species = "mouse"))
  check$Suggested.Symbol <- sapply(strsplit(check$Suggested.Symbol, " /// "), `[`, 1)
  dat <- dat[!is.na(check$Suggested.Symbol), ]
  check <- check[!is.na(check$Suggested.Symbol), ]
  rownames(dat) <- check$Suggested.Symbol[rownames(dat)%in%check$x]
  
  message("Combining counts...")
  duplicated_genes <- unique(rownames(dat)[duplicated(rownames(dat))])
  for (gene in duplicated_genes) {
    gene_indices <- which(rownames(dat) == gene)
    dat[gene_indices[1], ] <- colSums(dat[gene_indices, ])
    dat <- dat[-gene_indices[-1], , drop = FALSE]
  }
  
  message("Denoising and detecting doublets...")
  sce <- SingleCellExperiment(list(counts=dat))
  # remove empty cells
  lib.sizes <- colSums(counts(sce))
  if(sum(lib.sizes == 0) > 0) message(sum(lib.sizes == 0), " empty droplets detected.")
  sce <- sce[, lib.sizes > 0]
  sce <- decontX(sce)
  # doublet finder
  sce <- scDblFinder(sce)
  
  message("Creating Seurat Object...")
  srt <- CreateSeuratObject(counts = round(decontXcounts(sce)), 
                            project = sample_names[i], 
                            min.cells = 1, 
                            min.features = 1)
  
  # Add doublet metadata
  srt <- AddMetaData(srt, sce$scDblFinder.class, col.name="scDblFinder.class")
  srt <- AddMetaData(srt, sce$scDblFinder.score, col.name="scDblFinder.score")
  srt <- AddMetaData(srt, sce$scDblFinder.weighted, col.name="scDblFinder.weighted")
  srt <- AddMetaData(srt, sce$scDblFinder.cxds_score, col.name="scDblFinder.cxds_score")

  
  cat("Creating Seurat Object for", sample_names[i], "...\n")
  samples[[i]] <- srt
  samples[[i]]$genotype <- genotype[i]
  samples[[i]]$treatment <- treatment[i]
  samples[[i]]$sample_id <- sample_names[i]
  samples[[i]]$group <- gsub("1|2", "", sample_names[i])
  rm(dat, check, genes)
}

names(samples) <- sample_names
# save list of raw seurat objects
saveRDS(samples, paste0(output.path,"mPCa_list.rds"))

##### QC #####
samples <- lapply(samples, function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^(Mt|mt)-", assay="RNA")
  x[["percent.rb"]] <- PercentageFeatureSet(x, pattern = "^(Rpl|Rps)", assay="RNA")
  x[["percent.hb"]] <- PercentageFeatureSet(x, pattern = "^(Hbb|Hba)-", assay="RNA")
  return(x)
})

# merge object and save raw
mPCa <- merge(samples[[1]], samples[2:length(samples)])
mPCa <- JoinLayers(mPCa)
mPCa$group<- factor(mPCa$group, levels=c("DKO", "TKO", "TKO_CSI"))
mPCa$sample_id <- factor(mPCa$sample_id, levels=c("DKO1", "DKO2", "TKO1", "TKO2", "TKO_CSI1", "TKO_CSI2"))
saveRDS(mPCa, paste0(output.path,"mPCa_raw.rds"))

# filtering
cells <- mPCa@meta.data
cells$sample_id <- factor(cells$sample_id, levels=c("DKO1", "DKO2", "TKO1", "TKO2", "TKO_CSI1", "TKO_CSI2"))
cells$scDblFinder.class <- factor(cells$scDblFinder.class, levels = c("singlet", "doublet"))
summary_preQC <- cells %>% group_by(sample_id) %>% summarize(
  ncells = n(),
  mean_nCount = mean(nCount_RNA),
  median_nCount = median(nCount_RNA),
  mean_nFeature = mean(nFeature_RNA),
  median_nFeature = median(nFeature_RNA),
  mean_percent.mt = mean(percent.mt),
  median_percent.mt = median(percent.mt),
  doublet_rate = sum(scDblFinder.class == "doublet")/ncells,
)

# Filter doublets
cells <- filter(cells, scDblFinder.class == 'singlet')
# Filter cells
median_nCount <- median(cells$nCount_RNA) # 7936
mad_nCount <- mad(cells$nCount_RNA) # 4280
median_nFeature <- median(cells$nFeature_RNA) # 2544
mad_nFeature <- mad(cells$nFeature_RNA) # 1300
median_mt <- median(cells$percent.mt) # 1.79
mad_mt <- mad(cells$percent.mt) #2.40
# thresholds
thresholds_nCount <- c(0, median_nCount + 5*mad_nCount) # 29337
thresholds_nFeature <- c(100, median_nFeature + 5*mad_nFeature) #9034 # remove cells with less than 100 genes
thresholds_mt <- c(0, 25) # 25

qc_metrics <- data.frame(nCount=c(median_nCount, mad_nCount, thresholds_nCount), nFeature=c(median_nFeature, mad_nFeature, thresholds_nFeature), percent.mt=c(median_mt, mad_mt, thresholds_mt))
rownames(qc_metrics) <- c("Median", "MAD", "Lower Threshold", "Upper Threshold")

# filter cells
cells_filtered <- cells %>%
  dplyr::filter(
    nCount_RNA >= thresholds_nCount[1],
    nCount_RNA <= thresholds_nCount[2],
    nFeature_RNA >= thresholds_nFeature[1],
    nFeature_RNA <= thresholds_nFeature[2],
    percent.mt >= thresholds_mt[1],
    percent.mt <= thresholds_mt[2]
  )

cat("\nCells remaining after filtering:\n")
table(cells_filtered$sample_id)

# Create summary_postQC with proper handling of empty groups
summary_postQC <- cells_filtered %>% 
  group_by(sample_id) %>% 
  summarize(
    ncells = n(),
    mean_nCount = ifelse(ncells > 0, mean(nCount_RNA), NA),
    median_nCount = ifelse(ncells > 0, median(nCount_RNA), NA),
    mean_nFeature = ifelse(ncells > 0, mean(nFeature_RNA), NA),
    median_nFeature = ifelse(ncells > 0, median(nFeature_RNA), NA),
    mean_percent.mt = ifelse(ncells > 0, mean(percent.mt), NA),
    median_percent.mt = ifelse(ncells > 0, median(percent.mt), NA),
    doublet_rate = ifelse(ncells > 0, sum(scDblFinder.class == "doublet")/ncells, NA),
  )

fwrite(summary_preQC, paste0(table.path, "summary_preQC.csv"), quote=F)
fwrite(summary_postQC, paste0(table.path, "summary_postQC.csv"), quote=F)
# filter seurat object
mPCa <- subset(mPCa, cells=rownames(cells_filtered))
saveRDS(mPCa, paste0(output.path,"mPCa_filtered.rds"))

##### Clustering #####
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_s =checkGeneSymbols(mmus_s, species = "mouse")$Suggested.Symbol
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m=checkGeneSymbols(mmus_g2m, species = "mouse")$Suggested.Symbol
mPCa <- NormalizeData(mPCa) %>% CellCycleScoring(s.features = mmus_s, g2m.features = mmus_g2m)

mPCa <- SCTransform(mPCa, vars.to.regress = c("percent.mt", "G2M.Score", "S.Score")) %>% 
  RunPCA(dims=1:50)%>% 
  FindNeighbors(dims=1:50) %>%
  FindClusters(resolution=1, group.singletons = FALSE) %>%
  RunTSNE(dims = 1:50)

# no significant batch effects
DimPlot(mPCa, reduction="tsne")+DimPlot(mPCa, reduction="tsne", group.by = "sample_id")

markers <- FindAllMarkers(mPCa, min.pct = 0.25, logfc.threshold = 0.25, only.pos = F, return.thresh=0.01, test.use = "wilcox")
top10 <- markers %>% filter(!grepl("Rik$|^Gm", gene)) %>%  group_by(cluster) %>% slice_max(order_by = avg_log2FC, n=10)
print(top10, n=300)

fwrite(markers, paste0(table.path, "cluster_markers.csv"), row.names=FALSE)
saveRDS(mPCa, paste0(output.path,"mPCa_clustered.rds"))

##### Annotation #####
# singleR
immgen <- celldex::ImmGenData() # get immgen data from celldex
# get annotations
annot_singler <- SingleR(mPCa[["RNA"]]$data, 
                         ref=immgen,
                         labels=immgen$label.fine,
                         cluster=mPCa$seurat_clusters)
# transfer labels to original object
cluster<-annot_singler@rownames
class<-annot_singler$labels
names(class)<-cluster
mPCa[['celltype_singler']]<-unname(class[mPCa@meta.data$seurat_clusters])
table(mPCa$celltype_singler)
DimPlot(mPCa, group.by=c("celltype_singler", "seurat_clusters"),  label=T)&NoLegend()

# Azimuth
# predictions were run using online version of azimuth. but would have produced identical results with:
#Azimuth::AzimuthApp(reference="mouse_pansci")
# web app: https://app.azimuth.hubmapconsortium.org/app/mouse-pansci
pred <- read.delim("external data/azimuth_pred.tsv",header = T)
head(pred)
mPCa$celltype_azimuth <- pred$predicted.Main_cell_type[match(Cells(mPCa), pred$cell)]
DimPlot(mPCa, group.by=c("celltype_azimuth", "seurat_clusters"), label=T)&NoLegend()

# LLM-based
gsea <- run_gsea(markers, species = "mouse", gene_set = "M8") # M8, cell type gene sets

# GPT-4o predictions
res <- run_aiCelltype(mPCa,
  species = "mouse",
  tissue = "prostate",
  gsea_data = gsea,
  de_markers = markers,
  api_key = Sys.getenv("OPENAI_API_KEY"),
  top_n = 20)
fwrite(res, paste0(table.path, "gpt4o_celltype_predictions.csv"), row.names=TRUE)
labels <- res$predicted_label
names(labels) <- as.character(res$cluster)
labels <- labels[as.character(mPCa$seurat_clusters)]
names(labels) <- NULL
mPCa$celltype_gpt <- labels
DimPlot(mPCa, group.by = c("celltype_gpt","seurat_clusters"), label=T, repel=T)

# manual annotation
# T cell
DotPlot(mPCa, features=c("Cd3d", "Cd3e"))
treg <- c("Il2ra", "Foxp3", "Ctla4", "Il7r") # Cd25=Il2ra, Cd127=Il7r
cd8t <- c("Cd8a", "Ccr7", "Cd28", "Ifng", "Gzmb")
gdt <- c("Trgv2", "Il23r", "Il17a", "Rorc")
cd4t <- c("Cd4", "Cxcr3", "Ccr4", "Ccr3", "Ccr6", "Ccr10", "Cxcr5", "Icos")
nkt <- c("Nkg7", "Klrb1c", "Il4")
DotPlot(mPCa, features=c("Cd3d", cd8t, nkt, gdt, treg, cd4t))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# NK cell
DotPlot(mPCa, c("Ptprc","Itgam", "Il2rb", "Klrb1c", "Klrk1", "Ncr1")) # Cd11b=Itgam, Cd122=Il2rb, Nk1.1=Klrb1c, Nkg2d=Klrk1

# B cell
DotPlot(mPCa, c("Ptprc","Cd19", "Ms4a1")) # cd20=ms4a1
naive.b <- c("Ighm","Ighd", "Ccr7", "Sell")
DotPlot(mPCa, naive.b)
plasma <- c("Jchain", "Ighg2c", "Ighg2b", "Ighg1", "Ighg3")
DotPlot(mPCa, plasma)

# ILC (non-B, non-T lymphoid)
ilc1 <- c("Tbx21", "Ncr1", "Il7r", "Il2rb", "Il2rg", "Cd160")
DotPlot(mPCa, ilc1)
ilc2 <- c("Gata3", "Rora", "Klrg1", "Il17rb", "Il1rl1", "Icos", "Ptgdr2")
DotPlot(mPCa, ilc2)
ilc3 <- c("Rorc","Ncr1", "Il7r","Csf2ra", "Ltbr")
DotPlot(mPCa, ilc3)

# granolucytes
neu <- c("Ptprc", "Ly6g", "Mpo", "S100a8", "Csf3r", "Cd177")
DotPlot(mPCa, neu)
eos<- c("Ptprc", "Siglecf", "Ccr3", "Il5ra", "Prg2", "Epx")
DotPlot(mPCa, eos)
bas <- c("Ptprc", "Mcpt8", "Fcer1a", "Il4", "Ccr3", "Cd200r3")
DotPlot(mPCa, bas)

# myeloid
macro <- c("Ptprc", "Cd68", "Adgre1", "Mertk", "Cd163", "Mrc1")  #f4/80=adgre1
DotPlot(mPCa, macro)
mono <- c("Ptprc", "Ly6c2", "Cd14", "Itgam", "Fcgr1", "Ccr2", "Fcgr3")  #fcgr1=cd64
DotPlot(mPCa, mono)
DotPlot(mPCa, c("Ptprc", "H2-Ab1", "H2-Aa", "Itgax", "Cd209a", "Xcr1", "Zbtb46"))
pdc <- c("Clec4a4","Il3ra", "Tcf4", "Irf7", "Nrp1", "Lilrb4a")
DotPlot(mPCa, pdc)
cdc1 <- c("Batf3", "Tnfrsf11a", "Ccr7")
DotPlot(mPCa, features=cdc1)
cdc2 <- c("Itgam","Sirpa", "Clec10a", "Irf4", "Ccr7")
DotPlot(mPCa, cdc2)

# stromal, epithelial
DotPlot(mPCa, c("Epcam", "Ptprc", "Vwf", "Acta2", "Dcn", "Col1a1"))

celltypes <- c("Neutrophil", "Neutrophil", "T_CD8", "Neutrophil", "Epi", "T_gd", # 0-5
        "B", "T_gd", "T_CD4", "NK", "Macro/Mono", # 6-10
        "Neutrophil", "T_reg", "DC", "T_CD8", "NK", # 11-15
        "Fibro", "ILC", "DC", "DC") # 16-19

mPCa$celltype_final <- "Others"
for(i in seq_along(unique(mPCa$seurat_clusters))){
  mPCa$celltype_final[mPCa$seurat_clusters==as.character(i-1)] <- celltypes[as.numeric(i)]
}
DimPlot(mPCa, group.by=c("celltype_final", "celltype_gpt", "celltype_singler", "celltype_azimuth", "seurat_clusters"), label=T)&NoLegend()
mPCa$celltype_final <- factor(mPCa$celltype_final, levels=c("Epi", "Fibro", "Macro/Mono", "DC", "NK", "ILC", "T_CD4", "T_CD8", "T_gd", "T_reg", "B", "Neutrophil"))

Idents(mPCa) <- "celltype_final"
markers <- FindAllMarkers(mPCa, min.pct = 0.25, logfc.threshold = 0.25, only.pos = F, return.thresh=0.01, test.use = "wilcox")
fwrite(markers, paste0(table.path, "celltype_markers.csv"), row.names=FALSE)

##### Module Scores #####
chemo12 <- c("Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl8", "Ccl19", "Ccl21a", "Cxcl9", "Cxcl10", "Cxcl11", "Cxcl13", 
             "Ccl17") # ccl18 does not have analog, functional analog could be Ccl8 and Ccl17
mPCa <- AddModuleScore(mPCa, features=list(chemo12), name = "12Chemo")
mPCa$'12-Chemokine' <- mPCa$'12Chemo1'
mPCa$'12Chemo1' <- NULL
# tertiary lymphoid structure (TLS) Imprint (Maylan et al)
tls.mouse <- c(
  # Immunoglobulin genes
  "Igha", "Ighg1", "Ighg2b", "Ighg3", "Ighg2c", "Ighm", "Igkc", "Iglc1", "Iglc2", "Iglc3", "Jchain",
  # B cell markers
  "Cd79a", "Fcrl5", "Mzb1", "Ssr4", "Xbp1",
  # T cell markers
  "Trbc2", "Il7r",
  # Fibroblast markers
  "Cxcl12", "Lum",
  # Complement protein coding genes
  "C1qa", "C7",
  # Other genes
  "Cd52", "Apoe", "Hacd4", "Ptges3", "Pim2", "Derl3" # Ptplad2=hacd4,ptgds=ptgs3 
)
mPCa <- AddModuleScore(mPCa, features=list(tls.mouse), name = "TLS")
mPCa$'TLS-imprint' <- mPCa$TLS1
mPCa$TLS1 <- NULL
# MHC/IG score
mhc <- grep("^H2-", rownames(mPCa), value=T)
mhc1 <- grep("^H2-(K|D|L|T|Q|M)", rownames(mPCa), value=T)
mhc2 <- grep("^H2-(A|E)", rownames(mPCa), value=T)
ig <- grep("^Igl|Igk|Igh", rownames(mPCa), value=T)
ig <- ig[-which(ig=="Iglon5"|ig=="Ighmbp1")]

mPCa <- AddModuleScore(mPCa, list(mhc), name = "MHC-")
mPCa$'MHC' <- mPCa$`MHC-1`
mPCa$`MHC-1` <- NULL
mPCa <- AddModuleScore(mPCa, list(mhc1), name = "MHC-I")
mPCa$'MHC-I' <- mPCa$`MHC-I1`
mPCa$`MHC-I1` <- NULL
mPCa <- AddModuleScore(mPCa, list(mhc2), name = "MHC-II")
mPCa$'MHC-II' <- mPCa$`MHC-II1`
mPCa$`MHC-II1` <- NULL
mPCa <- AddModuleScore(mPCa, list(ig), name = "Ig-")
mPCa$'Ig' <- mPCa$`Ig-1`
mPCa$`Ig-1` <- NULL

##### Save Seurat objects and subsets for plotting #####
# Save main annotated object
saveRDS(mPCa, paste0(output.path,"mPCa_annotated.rds"))

# Save T cell subset
t_cells <- subset(mPCa, celltype_final %in% c("T_CD8", "T_CD4", "T_gd", "T_reg"))
saveRDS(t_cells, paste0(output.path, "mPCa_t_cells.rds"))

# Save B cell subset
b_cells <- subset(mPCa, celltype_final == "B")
saveRDS(b_cells, paste0(output.path, "mPCa_b_cells.rds"))

sessionInfo()
