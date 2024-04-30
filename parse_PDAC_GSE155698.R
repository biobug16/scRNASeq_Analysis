# Parse PDAC_GSE155698
# Multimodal mapping of the tumor and peripheral blood immune landscape in human pancreatic cancer
# https://www.nature.com/articles/s43018-020-00121-4

# Load libraries
library(MAESTRO)
library(Seurat)
library(ggplot2)
library(data.table)
library(Matrix)
library(rhdf5)
library(HDF5Array)
library(GEOquery)

# Set directoies ---------------------------------------------------------------
data_dir <- "/data/Lab/resources/single_cell_datasets/GSE155698_PDAC"

# Set working directory --------------------------------------------------------
setwd(data_dir)

# Read in count matrix
pdac <- getGEO("GSE155698")

# Extract relevant information -------------------------------------------------
# GPL20301
pdac_GPL20301 <- pdac$`GSE155698-GPL20301_series_matrix.txt.gz`@phenoData@data
pdac_GPL20301 <- pdac_GPL20301[c("title", "source_name_ch1", "treatment:ch1", "instrument_model")]
#GPL20301
pdac_GPL24676 <- pdac$`GSE155698-GPL24676_series_matrix.txt.gz`@phenoData@data
pdac_GPL24676 <- pdac_GPL24676[c("title", "source_name_ch1", "treatment:ch1", "instrument_model")]

# merge datasets
pdac_combined <- rbind(pdac_GPL20301, pdac_GPL24676)

# Collect samples --------------------------------------------------------------
seurats <- list()

# Remove problematic samples
remove_samples <- c("PDAC_PBMC_16", "PDAC_TISSUE_14", "PDAC_PBMC_14", "PDAC_TISSUE_16")
pdac_combined <- pdac_combined[!pdac_combined$title %in% remove_samples, ]
number_of_samples <- length(pdac_combined$title)

for (i in 1:number_of_samples){

  dir <- pdac_combined$title[i]
  print(dir)
  print(paste(dir, "/filtered_feature_bc_matrix/", sep=""))
  seurat_data <- Seurat::Read10X(data.dir = paste(dir, "/filtered_feature_bc_matrix/", sep=""))
  seurats[[dir]] <- Seurat::CreateSeuratObject(counts = seurat_data, project = dir, )
  
}

# Merge to single seurat object ------------------------------------------------
seurat_PBMC <- merge(
  seurats$PDAC_PBMC_13, 
  y = list(seurats$PDAC_PBMC_15, seurats$PDAC_PBMC_1, seurats$PDAC_PBMC_2,
           seurats$PDAC_PBMC_3, seurats$PDAC_PBMC_4, seurats$PDAC_PBMC_5,
           seurats$PDAC_PBMC_6, seurats$PDAC_PBMC_7, seurats$PDAC_PBMC_8,
           seurats$PDAC_PBMC_9, seurats$PDAC_PBMC_10A, seurats$PDAC_PBMC_10B,
           seurats$PDAC_PBMC_11, seurats$PDAC_PBMC_12, seurats$Healthy_PBMC_1, 
           seurats$Healthy_PBMC_2, seurats$Healthy_PBMC_3, seurats$Healthy_PBMC_4), 
  add.cell.ids = pdac_combined$title[grep("PBMC", pdac_combined$title)], 
  project = "PDAC_GSE155698_PBMC")

seurat_TISSUE <- merge(
  seurats$PDAC_TISSUE_13, 
  y = list(seurats$PDAC_TISSUE_1, seurats$PDAC_TISSUE_2, seurats$PDAC_TISSUE_3, 
           seurats$PDAC_TISSUE_4, seurats$PDAC_TISSUE_5, seurats$PDAC_TISSUE_6, 
           seurats$PDAC_TISSUE_7, seurats$PDAC_TISSUE_8, seurats$PDAC_TISSUE_9, 
           seurats$PDAC_TISSUE_10, seurats$PDAC_TISSUE_11A, seurats$PDAC_TISSUE_11B, 
           seurats$PDAC_TISSUE_12, seurats$PDAC_TISSUE_15, seurats$AdjNorm_TISSUE_1, 
           seurats$AdjNorm_TISSUE_2, seurats$AdjNorm_TISSUE_3), 
  add.cell.ids = pdac_combined$title[grep("TISSUE", pdac_combined$title)], 
  project = "PDAC_GSE155698_TISSUE")

# PDAC_GSE155698_PBMC ----------------------------------------------------------
# Clustering
RNA_res <- MAESTRO::RNARunSeurat(
  inputMat = seurat_PBMC,
  type = "object",
  project = "PDAC_GSE155698_PBMC", 
  min.c = 10,
  min.g = 500,
  dims.use = 1:30,
  variable.genes = 2000, 
  organism = "GRCh37",
  cluster.res = 1,
  genes.test.use = "presto",
  only.pos = TRUE,
  genes.cutoff = 1e-05)

# Cell-type annotation
RNA_res$RNA <- MAESTRO::RNAAnnotateCelltype(
  RNA = RNA_res$RNA, 
  genes = RNA_res$genes,
  signatures = "human.immune.CIBERSORT",
  min.score = 0.6)

# Check batch effect -----------------------------------------------------------
p <- Seurat::DimPlot(
  object = RNA_res$RNA, 
  group = "orig.ident", 
  label = TRUE, 
  pt.size = 0.1)

ggsave(
  file.path(paste0(RNA_res$RNA@project.name, "_", "orig.ident", ".png")), 
  p, width=6.5, height=5)

p <- Seurat::DimPlot(
  object = RNA_res$RNA, 
  group = "assign.ident", 
  label = TRUE, pt.size = 0.1)

ggsave(
  file.path(paste0(RNA_res$RNA@project.name, "_", "assign_ident", ".png")), 
  p, width=6.5, height=5)

# Reassign pDC to DC -----------------------------------------------------------
RNA_res$RNA@meta.data[RNA_res$RNA@meta.data$assign.ident=="pDC", ]$assign.ident <- "DC"

# Add sample_type information
RNA_res$RNA@meta.data$sample_type <- "PBMC"

# Save expression
saveRDS(RNA_res$RNA@assays$RNA@data, "PBMC_GSE155698_expression.rds")

# Save CellMetainfo_table
umap_df <- as.data.frame(RNA_res$RNA@reductions$umap@cell.embeddings)
meta_col <-  c("assign.ident", "seurat_clusters", "orig.ident", "sample_type")
meta_df <- RNA_res$RNA@meta.data[ ,meta_col]
colnames(meta_df) = c("Celltype (major-lineage)", "Cluster", "Sample", "Tissue")
meta_df$Cluster <- as.character(meta_df$Cluster)

# save to _CellMetainfo_table.tsv
meta_all <- cbind(Cell=row.names(meta_df), umap_df, meta_df)
write.table(meta_all, "PBMC_GSE155698_CellMetainfo_table.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

# PDAC_GSE155698_TISSUE ----------------------------------------------------------
# Clustering
RNA_res <- MAESTRO::RNARunSeurat(
  inputMat = seurat_TISSUE,
  type = "object",
  project = "PDAC_GSE155698_TISSUE", 
  min.c = 10,
  min.g = 500,
  dims.use = 1:30,
  variable.genes = 2000, 
  organism = "GRCh37",
  cluster.res = 1,
  genes.test.use = "presto",
  only.pos = TRUE,
  genes.cutoff = 1e-05)

# Cell-type annotation
RNA_res$RNA <- MAESTRO::RNAAnnotateCelltype(
  RNA = RNA_res$RNA, 
  genes = RNA_res$genes,
  signatures = "human.immune.CIBERSORT",
  min.score = 0.6)

# Check batch effect -----------------------------------------------------------
p <- Seurat::DimPlot(
  object = RNA_res$RNA, 
  group = "orig.ident", 
  label = TRUE, 
  pt.size = 0.1)

ggsave(
  file.path(paste0(RNA_res$RNA@project.name, "_", "orig.ident", ".png")), 
  p, width=6.5, height=5)

p <- Seurat::DimPlot(
  object = RNA_res$RNA, 
  group = "assign.ident", 
  label = TRUE, pt.size = 0.1)

ggsave(
  file.path(paste0(RNA_res$RNA@project.name, "_", "assign_ident", ".png")), 
  p, width=6.5, height=5)

# Reassign pDC to DC -----------------------------------------------------------
RNA_res$RNA@meta.data[RNA_res$RNA@meta.data$assign.ident=="pDC", ]$assign.ident <- "DC"

# Add sample_type information
RNA_res$RNA@meta.data$sample_type <- "TISSUE"

# Save expression
saveRDS(RNA_res$RNA@assays$RNA@data, "PDAC_GSE155698_expression.rds")

# Save CellMetainfo_table
umap_df <- as.data.frame(RNA_res$RNA@reductions$umap@cell.embeddings)
meta_col <-  c("assign.ident", "seurat_clusters", "orig.ident", "sample_type")
meta_df <- RNA_res$RNA@meta.data[ ,meta_col]
colnames(meta_df) = c("Celltype (major-lineage)", "Cluster", "Sample", "Tissue")
meta_df$Cluster <- as.character(meta_df$Cluster)

# save to _CellMetainfo_table.tsv
meta_all <- cbind(Cell=row.names(meta_df), umap_df, meta_df)
write.table(meta_all, "PDAC_GSE155698_CellMetainfo_table.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)



