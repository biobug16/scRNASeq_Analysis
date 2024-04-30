# Parse SC2018 dataset
# Single-cell transcriptomics reveals expansion of cytotoxic CD4 T cells in supercentenarians
# https://www.pnas.org/content/116/48/24242

# Load libraries
library(MAESTRO)
library(Seurat)
library(ggplot2)
library(data.table)
library(Matrix)
library(rhdf5)
library(HDF5Array)
library(SeuratDisk)

# Read in count matrix
data_dir <- "/data/Lab/resources/single_cell_datasets/SC2018/"
sc_raw <- data.table::fread(file=paste(data_dir, "01.UMI.txt", sep=""))

# Read in meta data
meta_data <- read.csv(
  paste(data_dir, "03.Cell.Barcodes.txt", sep=""), 
  sep="\t", header=FALSE)

row.names(meta_data) <- meta_data$V1
meta_data <- meta_data[, 2:3]
colnames(meta_data) <- c("sample_id", "sample_type")

# Get gene name column
features_ensembl <- sc_raw[, 1]$V1

# Remove first column
sc_raw <- sc_raw[, -1]

# Read in full gene names (?)
# For some reason the full matrix has less rows..., some genes 
# where removed.
features_gs <- data.table::fread(
  file=paste(data_dir, "SC1/features.tsv.gz", sep=""),
  header=FALSE)

# Only keep rownames that are present in data matrix
features_tokeep <- features_gs[features_gs$V1 %in% features_ensembl, ]

# Convert to matrix
sc_raw <- as.matrix(sc_raw)

# Add row names
row.names(sc_raw) <- features_tokeep$V2

# Create Seurat object
seuratObj <- Seurat::CreateSeuratObject(
  counts=sc_raw, 
  meta.data=meta_data)

# Free up some space
rm(sc_raw)

# Clustering
RNA_res <- MAESTRO::RNARunSeurat(
  inputMat = seuratObj,
  type = "object",
  project = "PBMC_SC2018", 
  min.c = 10,
  min.g = 500,
  dims.use = 1:30,
  variable.genes = 2000, 
  organism = "GRCh38",
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

# Reassign pDC to DC
RNA_res$RNA@meta.data[RNA_res$RNA@meta.data$assign.ident=="pDC", ]$assign.ident <- "DC"

# Save expression
saveRDS(RNA_res$RNA@assays$RNA@data, "PBMC_SC2018_expression.rds")

# Save CellMetainfo_table
umap_df <- as.data.frame(RNA_res$RNA@reductions$umap@cell.embeddings)
meta_col <-  c("assign.ident", "seurat_clusters", "sample_id", "sample_type")
meta_df <- RNA_res$RNA@meta.data[ ,meta_col]
colnames(meta_df) = c("Celltype (major-lineage)", "Cluster", "Sample", "Tissue")
meta_df$Cluster <- as.character(meta_df$Cluster)

# save to _CellMetainfo_table.tsv
meta_all <- cbind(Cell=row.names(meta_df), umap_df, meta_df)
write.table(meta_all, "PBMC_SC2018_CellMetainfo_table.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

# Check batch effect
p <- DimPlot(
  object = RNA_res$RNA, 
  group = "sample_id", 
  label = TRUE, 
  pt.size = 0.1)

ggsave(
  file.path(paste0(RNA_res$RNA@project.name, "_", "sample_id", ".png")), 
  p, width=6.5, height=5)

p <- DimPlot(
  object = RNA_res$RNA, 
  group = "assign.ident", 
  label = TRUE, pt.size = 0.1)

ggsave(
  file.path(paste0(RNA_res$RNA@project.name, "_", "assign_ident", ".png")), 
  p, width=6.5, height=5)



