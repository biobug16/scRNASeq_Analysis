#!/usr/bin/env Rscript

# If warn is one, warnings are printed as they occur
options(warn=1)

# Set debug variable
debug <- TRUE

# Load libraries
library(rhdf5)
library(Seurat)
library(tidyverse)

# Parse arguments
parser <- ArgumentParser()

parser$add_argument(
  "--seurat-obj", 
  type="character", 
  help="Seurat object.", 
  required=FALSE
)

parser$add_argument(
  "--h5",
  type="character",
  help="h5 file.",
  required=FALSE
)

parser$add_argument(
  "--meta",
  type="character",
  help="Meta cell and sample information.",
  required=FALSE
)

parser$add_argument(
  "--genes", 
  type="character", 
  help="List of genes.", 
  required=TRUE
)

parser$add_argument(
  "--filter-low",
  type="boolean",
  help="Filter low expressed genes (Recommended for Smart-seq2).",
  required=FALSE
)

parser$add_argument(
  "--output",
  type="character",
  help="Output directory",
  required=TRUE
)

args <- parser$parse_args()

if(debug){print(args)}

# Set variables ----------------------------------------------------------------
# Example dataset
#f1 <- "LIHC_GSE125449_aPDL1aCTLA4_expression.h5"
#f2 <- "LIHC_GSE125449_aPDL1aCTLA4_CellMetainfo_table.tsv"

# Set parameters ---------------------------------------------------------------
# Filter for low expressed genes
# Take lower 10% of genes
lower_prob <- 0.05
# Scale factor
scale_factor <- 10000
# Cell types
# Careful, Mono/Macro --> Mono.Macro
custom <- TRUE
custom_cell_types <- c("Mono/Macro", "CD4Tconv", "CD8T", "B", "DC", "NK")

# Print arguments and parameters -----------------------------------------------
# print(args)
#print("Lower prob: ", lower_prob)
#print("Scale factor:", scale_factor)

# Data directory
data_dir <- "/data/Lab/resources/single_cell_datasets/TISCH/raw_downloads/"

metap <- read.table(paste(data_dir, "Parsed_meta_4.csv", sep=""), sep=",", header=TRUE)
filt <- tibble(metap) %>%
  filter(!cancer %in% c("AEL", "ALL", "AML", "MM")) %>% # exclude haematopoietic cancers
  #filter(other %in% c("Smartseq2", "10X", "CellMetainfo")) %>% # only include primary cancers
  filter(complete.cases(Mono.Macro, CD4Tconv, CD8T, B, DC, NK)) %>%  # Tregs not in healthy PBMC !
  filter(Mono.Macro > 5 & CD4Tconv > 5 & CD8T > 5 & B > 5 & DC > 5 & NK > 5)

# Get files
#meta_files <- unique(filt$file_name)
#expr_files <- gsub("CellMetainfo_table.tsv", "expression.h5", meta_files)
#prefix <- gsub("_CellMetainfo_table.tsv", "", meta_files)

# Only process SC2018
#meta_files <- "PBMC_SC2018_CellMetainfo_table.tsv"
#expr_files <- "PBMC_SC2018_expression.rds"
#prefix <- gsub("_CellMetainfo_table.tsv", "", meta_files)

# Only process PBMC_GSE155698
meta_files <- "PBMC_8K_10X_CellMetainfo_table.tsv"
expr_files <- "PBMC_8K_10X_expression.h5"
prefix <- gsub("_CellMetainfo_table.tsv", "", meta_files)

# Read in data -----------------------------------------------------------------
args <- list()
args$filter_low <- FALSE
args$output <- "/data/Lab/mehanib2/RNA-seq-immune-sig/results/"

for (p in c(3:4)){
  
  args$h5 <- paste(data_dir, expr_files[p], sep="")
  args$meta <- paste(data_dir, meta_files[p], sep="")
  args$prefix <- paste(prefix[p], sep="")
  
  print(args$h5)
  print(args$meta)
  print(args$prefix)

  #fnum <-  1
  #args$h5 <- paste(data_dir, expr_files[fnum], sep="")
  #args$meta <- paste(data_dir, meta_files[fnum], sep="")
  
  # Get patients
  #patients <- filter(filt, file_name==meta_files[p])$patient_name
  #patients <- filter(filt, dataset=="SC2018")$patient_name
  #patients <- c("SC1", "SC2", "SC3", "SC4", "SC5", "SC6", "SC7")
  #patients <- filter(filt, dataset=="GSE155698")$patient_name
  patients <- c("S1")
  
  # Load data ------------------------------------------------------------------
  loaded_data <- load_data(args, patients, custom_cell_types, id="Sample")
  
  expression <- loaded_data$expression
  meta_data <- loaded_data$meta_data
  
  # Filter low expressind genes -------------------------------------------------
  if (args$filter_low){
    
    # For each cell phenotype, genes with low average expression (<0.75 transcript per cell) in log2 
    # space were then set to 0 as a quality control filter. Although an expression threshold 
    # of <0.75 was used in this work, downstream results remained comparable when using modestly
    # different thresholds (data not shown).
    expression_filt <- filter_genes(expression, meta_data, patients, custom_cell_types)
  }
  
  # free up some memory
  rm(loaded_data)
  #rm(expression)
  
  # Get aggregate gene expression profiles (GEP) -------------------------------
  print("Start aggregating gene expression profiles")
  geps <- run_geps(expression, meta_data, patients, custom_cell_types, scale_factor=10000)

  # DEG ------------------------------------------------------------------------
  print("Start differential gene expression")
  deg_results <- run_diff_expr(geps, meta_data, patients, custom_cell_types)
  
  # Optimize G -----------------------------------------------------------------
  print("Start optimizing G")
  signature <- run_G(deg_results, meta_data, patients)
  
  # Save results ---------------------------------------------------------------
  saveRDS(signature$sig, paste(args$output, args$prefix, "_signature_matrix.rds", sep=""))
  saveRDS(deg_results, paste(args$output, args$prefix, "_differential_expression_results.rds", sep=""))
  saveRDS(signature$gene_sets, paste(args$output, args$prefix, "_gene_sets_results.rds", sep=""))
  saveRDS(signature$kappa, paste(args$output, args$prefix, "_kappa_results.rds", sep=""))

}


# FUNCTIONS --------------------------------------------------------------------

aggregateGEP <- function(
  matrix,
  meta,
  cell_types,
  scale_factor=scale_factor,
  seed = 1){
  
  # Get aggregated cell type expression.
  # ------------------------------------
  # For each cell type represented by at least three single cells, we selected 50% of all available 
  # single-cell GEPs using random sampling without replacement (fractional sample sizes were rounded 
  # up such that two cells were sampled if only three were available). We then aggregated the profiles 
  # by summation in non-log linear space and normalized each population-level GEP into TPM. This process 
  # was repeated in order to generate five aggregated transcriptome replicates per cell type. 
  # assume availabe cell types are > 3
  
  # set seed
  set.seed(seed)
  
  # Aggregated cell types
  agg_ct <- matrix(nrow=dim(matrix)[1], ncol=length(cell_types) * 5)
  matCol <- 1
  
  for (i in 1:length(cell_types)){
    
    print(cell_types[i])
    
    # Get cell type
    cell_type_specific <- meta[meta$`Celltype (major-lineage)`==cell_types[i], ]
    
    for(j in 1:5){
      
      # Randomly sample 50% 
      sampled_cells <- base::sample(
        cell_type_specific$Cell, 
        size=round(dim(cell_type_specific)[1]*0.5, digits = 0),
        replace=FALSE)
      
      # Aggregate the profiles by summation in non-log linear space
      cell_gep <- Matrix::rowSums(matrix[, sampled_cells])
      
      # Write result to new matrix
      agg_ct[, matCol] <- cell_gep
      matCol <- matCol + 1
      
    }
  }
  
  colnames(agg_ct) <- paste(sort(rep(cell_types, 5)), seq(1:5), sep=".")
  row.names(agg_ct) <- row.names(matrix)
  
  # ..and normalized each population-level GEP into TPM 
  # ?
  # Some how these numbers are off !
  agg_ct_Colsum <- t(scale_factor * (t(agg_ct) / Matrix::colSums(agg_ct)))
  
  return(agg_ct_Colsum)
}

# Function of wilcox test
wilcox_class <- function(
  v, 
  b, 
  alternative="two.sided"){
  
  # Get pvalue
  p <- wilcox.test(v[b],v[!b], alternative = alternative)$p.value
  
  # Get fold change
  # add small pseudocount
  vpseudo <- v + 1
  fc <- log2(mean(na.omit(vpseudo[b])/mean(na.omit(vpseudo[!b]))))
  
  return(c(p, fc))
}

# Remove rows with rowSums == 0
deg <- function(
  gep,
  cell_types){
  
  agg_ct_nz <- gep[rowSums(gep) != 0, ]
  
  diff_expr <- list()
  
  for (i in 1:length(cell_types)){
    
    # Get cell type
    print(cell_types[i])
    
    col_ct <- paste(rep(cell_types[i], 5), seq(1:5), sep=".")
    b <- colnames(gep) %in% col_ct
    
    # Run Wilcoxon tests
    suppressWarnings({
      results <- t(apply(agg_ct_nz, 1, function(x) wilcox_class(x, b, alternative="greater")))
    })
    
    colnames(results) <- c('p_value','log2fc')
    results <- as.data.frame(results)
    
    # Perform test
    results$BH <- p.adjust(results$p_value, "BH")
    diff_expr[[cell_types[i]]] <- results
  }
  return(diff_expr)
}


get_kappa <- function(
  agg_ct,
  diff_expr, 
  cell_types,
  G_min=300,
  G_max=500){
  
  kappa_cond_num <- c()
  
  for (G in G_min:G_max){
    
    print(G)
    
    sub_diff <- list()
    
    for(i in 1:length(cell_types)){
      
      diff_ct <- diff_expr[[cell_types[i]]]
      
      # Genes with a q value <0.01 (false discovery rate) were considered significant.
      diff_sig <- diff_ct[diff_ct$BH < 0.01, ]
      
      # Sort by logfc
      diff_sorted <- diff_sig[order(diff_sig$log2fc, decreasing=TRUE), ]
      
      # Select top G genes and calculate kappa
      diff_top <- head(diff_sorted, G)
      sub_diff[[cell_types[i]]] <- diff_top
    }
    
    # Combined dataframes to single matrix
    # Replace missing genes with 1
    
    # Get genes
    genes <- c()
    for(i in 1:length(cell_types)){
      print(cell_types[i])
      #print(length(rownames(sub_diff[[types]])))
      genes <- c(genes, rownames(sub_diff[[cell_types[i]]]))
    }
    print(length(genes))
    genes_uniq <- unique(sort(genes))
    print(length(genes_uniq))
    
    # Summarize matrix to single average expression
    av_expr <- data.frame(row.names = row.names(agg_ct))
    
    for(i in 1:length(cell_types)){
      col_ct <- paste(rep(cell_types[i], 5), seq(1:5), sep=".")
      av_expr[[cell_types[i]]] <- rowMeans(agg_ct[, col_ct])
    }
    av_expr_sel <- av_expr[genes, ]
    
    # For each column, set non-sig genes to 1
    #for(i in 1:length(cell_types)){
    #  types <- names(cell_types)[i]
    #  av_expr_sel[! row.names(av_expr_sel) %in% row.names(sub_diff[[types]]), types ] <- 1
    #}
    
    #print(dim(av_expr_sel))
    mat_kappa <- as.matrix(av_expr_sel)
    kappa_cond_num <- c(kappa_cond_num, kappa(mat_kappa))
  }
  
  kappa_df <- data.frame(
    gene_num = G_min:G_max,
    kappa = kappa_cond_num
  )
  
  kappa_df <- kappa_df[order(kappa_df$kappa), ]
  return(kappa_df)
}


get_signature_matrix <- function(
  agg_ct,
  diff_expr, 
  cell_types,
  G){
  
  sub_diff <- list()
  
  for(i in 1:length(cell_types)){
    
    diff_ct <- diff_expr[[cell_types[i]]]
    
    # Genes with a q value <0.01 (false discovery rate) were considered significant.
    diff_sig <- diff_ct[diff_ct$BH < 0.01, ]
    
    # Sort by logfc
    diff_sorted <- diff_sig[order(diff_sig$log2fc, decreasing=TRUE), ]
    
    # Select top G genes and calculate kappa
    diff_top <- head(diff_sorted, G)
    sub_diff[[cell_types[i]]] <- diff_top
  }
  
  # Combined dataframes to single matrix
  # Replace missing genes with 1
  
  # Get genes
  genes <- c()
  for(i in 1:length(cell_types)){
    print(cell_types[i])
    #print(length(rownames(sub_diff[[types]])))
    genes <- c(genes, rownames(sub_diff[[cell_types[i]]]))
  }
  print(length(genes))
  genes_uniq <- unique(sort(genes))
  print(length(genes_uniq))
  
  # Summarize matrix to single average expression
  av_expr <- data.frame(row.names = row.names(agg_ct))
  
  for(i in 1:length(cell_types)){
    col_ct <- paste(rep(cell_types[i], 5), seq(1:5), sep=".")
    av_expr[[cell_types[i]]] <- rowMeans(agg_ct[, col_ct])
  }
  av_expr_sel <- av_expr[genes, ]
  return(av_expr_sel)
}


get_gene_set <- function(
  diff_expr,
  cell_types,
  G){
  
  sub_diff <- list()
  
  for(i in 1:length(cell_types)){
    
    diff_ct <- diff_expr[[cell_types[i]]]
    
    # Genes with a q value <0.01 (false discovery rate) were considered significant.
    diff_sig <- diff_ct[diff_ct$BH < 0.01, ]
    
    # Sort by logfc
    diff_sorted <- diff_sig[order(diff_sig$log2fc, decreasing=TRUE), ]
    
    # Select top G genes and calculate kappa
    diff_top <- head(diff_sorted, G)
    sub_diff[[cell_types[i]]] <- diff_top
  }
  return(sub_diff)
}


load_data <- function(
  args, 
  patients,
  custom_cell_types,
  id="Sample"){
  
  # Load data --------------------------------------------------------------------
  h5 <- Seurat::Read10X_h5(args$h5)
  #h5 <- readRDS(args$h5)
  meta <- read_tsv(args$meta)
  
  # Select specific patients
  #pn <- filter(filt, cancer %in% "HNSC" & !patient_name %in% c("HD1", "HD2", "HD4", "HD5", "HD6"))$patient_name
  #half <- pn[12:25]
  cells <- filter(meta, Sample %in% patients)$Cell
  h5_selected <- h5[, cells]
  unlog <- exp(h5_selected) - 1
  
  # unlog normalized data
  # Result will be a matrix with equal colSums --> 10.000
  #unlog <- exp(h5) - 1
  
  # free up some memory
  rm(h5)
  
  # Create list of matrices and meta data, each entry one patient ----------------
  expression <- list()
  meta_data <- list()
  #patients <- half
  
  if(id=="Patient"){
    
    for (j in 1:length(patients)){
      
      print(patients[j])
      
      cells <- filter(
        meta, Patient %in% patients[j] & 
          `Celltype (major-lineage)` %in% custom_cell_types)$Cell
      expression[[patients[j]]] <- unlog[, cells]
      meta_data[[patients[j]]] <- as.data.frame(meta[meta$Patient %in% patients[j], ])
    }
  }
  
  if(id=="Sample"){
    
    for (i in 1:length(patients)){
      # Only select cells with the right 
      # cell type and patient name
      
      print(patients[i])
      
      cells <- filter(
        meta, Sample %in% patients[i] & 
          `Celltype (major-lineage)` %in% custom_cell_types)$Cell
      expression[[patients[i]]] <- unlog[, cells]
      meta_data[[patients[i]]] <- as.data.frame(meta[meta$Sample %in% patients[i], ])
    }
  }
  
  # free up some memory
  rm(unlog)
  
  return(list(expression=expression, meta_data=meta_data))
}  


filter_genes <- function(
  expression, 
  meta_data, 
  patients,
  custom_cell_types,
  custom=TRUE,
  lower_prob=0.05){
  
  expression_filt <- list()
  
  for (i in 1:length(patients)){
    
    if(custom){
      cell_types <- custom_cell_types
    }else{
      cell_types <- table(meta_data[[patients[i]]]$`Celltype (major-lineage)`)
    }
    
    # Only keep cell types with at least 3 cells
    #cell_types <- cell_types[cell_types > 2]
    
    meta_sub <- meta_data[[patients[i]]]
    
    # Collect low expresed genes for each cell phenotype
    genes <- integer()
    
    # verbose
    print(cell_types)
    
    for (j in 1:length(cell_types)){
      
      # For each cell phenotype ...
      cells <- meta_sub[meta_sub$`Celltype (major-lineage)` %in% cell_types[j], ]$Cell
      
      # Get row mean
      cells_pheno_mean <- Matrix::rowMeans(log2(expression[[patients[i]]][, cells] + 1))
      
      # Calculate the lowest 10%
      lower_percent <- quantile(cells_pheno_mean, probs=c(lower_prob))
      
      # What are the genes ?
      #print(head(genes))
      genes <- c(genes, which(cells_pheno_mean <= lower_percent))
      
      # Set genes to 0 with low average expression
      uniq_genes <- unique(genes)
      expression_filt[[patients[i]]][uniq_genes, cells] <- 0 
    }
  }
  
  return(expression_filt)
}


run_geps <- function(
  expression, 
  meta_data, 
  patients, 
  custom_cell_types, 
  scale_factor=10000){
  
  geps <- list()
  
  for (k in 1:length(patients)){
    
    print(patients[k])
    
    expr_matrix <- expression[[k]]
    meta <- meta_data[[k]]
    
    if(custom){
      cell_types <- custom_cell_types
    }else{
      cell_types <- table(meta_data[[patients[k]]]$`Celltype (major-lineage)`)
    }
    
    # Only keep cell types with at least 3 cells
    # cell_types <- cell_types[cell_types > 2]
    
    geps[[patients[k]]] <- aggregateGEP(expr_matrix, meta, cell_types = cell_types, scale_factor=10000)
  }
  
  return(geps)  
}


run_diff_expr <- function(
  geps, 
  meta_data, 
  patients, 
  custom_cell_types){
  
  deg_results <- list()
  
  for (i in 1:length(patients)){
    
    print(patients[i])
    
    matrix_gep <- geps[[i]]
    meta <- meta_data[[i]]
    
    if(custom){
      cell_types <- custom_cell_types
    }else{
      cell_types <- table(meta$`Celltype (major-lineage)`)
    }
    
    
    # Only keep cell types with at least 3 cells
    # cell_types <- cell_types[cell_types > 2]
    
    res <- deg(gep=as.matrix(matrix_gep), cell_types)
    deg_results[[patients[i]]] <- res
  }
  
  return(deg_results)
}


run_G <- function(
  deg_results,
  meta_data,
  patients,
  single_cell=TRUE,
  k_max=999,
  q_value=0.01,
  G_min=300,
  G_max=500,
  filter=FALSE){
  
  sig <- list()
  gene_sets <- list()
  
  for (n in 1:length(patients)){
    
    gep <- geps[[n]]
    deg_res <- deg_results[[n]]
    meta <- meta_data[[n]]
    
    if(custom){
      cell_types <- custom_cell_types
    }else{
      cell_types <- table(meta$`Celltype (major-lineage)`)
    }
    
    # Only keep cell types with at least 3 cells
    #cell_types <- cell_types[cell_types > 2]
    
    kappa_df <- get_kappa(
      agg_ct=gep, 
      diff_expr = deg_res, 
      cell_types = cell_types,
      G_min=300,
      G_max=500)
    
    min_g <- kappa_df[kappa_df$kappa == min(kappa_df$kappa), ]$gene_num
    
    # Create matrix with optimal number of genes
    sigmat <- get_signature_matrix(agg_ct=gep, diff_expr = deg_res, cell_types = cell_types, G=min_g)
    
    # Get gene sets (with logfc)
    gene_sets[[patients[n]]] <- get_gene_set(
      diff_expr=deg_res, 
      cell_types=cell_types,
      G=min_g)
    
    sig[[patients[n]]] <- sigmat
  }
  return(list(gene_sets=gene_sets, sig=sig, kappa=kappa_df))
}
