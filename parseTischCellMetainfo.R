# Parse files from TISCH/biowulf -------------

# Settings
verbose <- TRUE

# Load libraries
library(tidyverse)

# Read in meta data files
files <- read.csv("CellMetainfo_table_4.tsv", header=FALSE)

meta_data <- list()
for (i in 1:length(files$V1)){
  print(files$V1[i])
  meta_data[[files$V1[i]]] <- read.csv(files$V1[i], sep="\t")
  print(head(meta_data[[i]]))
}

patient_meta <- data.frame()

for (i in 1:length(meta_data)){
  
  # Filename
  fn <- files$V1[i]
  cancer <- unlist(strsplit(files$V1[i], split="_"))[1]
  dataset <- unlist(strsplit(files$V1[i], split="_"))[2]
  other <- unlist(strsplit(files$V1[i], split="_"))[3]
  
  if(verbose){
    print(fn)
    print(cancer)
    print(dataset)
    print(other)}
  
  # Meta
  meta <- meta_data[[i]]
  
  # Patients
  if("Patient" %in% colnames(meta)){
    patients <- table(meta$Patient)
  }else{
    patients <- table(meta$Sample)
  }
  
  if(verbose){print(patients)}
  
  for (patient in names(patients)){
    
    if(verbose){print(patient)}
    
    if("Patient" %in% colnames(meta)){
      patient_df <- meta[meta$Patient==patient, ]
    }else{
      patient_df <- meta[meta$Sample==patient, ]
    }
    
    parsed_data <- as.data.frame(table(patient_df$Celltype..major.lineage.))
    
    parsed_data$patient_name <- patient
    parsed_data$cancer <- cancer
    parsed_data$dataset <- dataset
    parsed_data$other <- other
    parsed_data$file_name <- fn
    
    if ("Age" %in% colnames(meta)){
      parsed_data$Age <- unique(patient_df$Age)[1]
    }else{
      parsed_data$Age <- NA
    }
    
    if ("Stage" %in% colnames(meta)){
      # Sometimes the following appears
      # > unique(patient_df$Stage)
      # [1] "Metastatic" NA
      parsed_data$Stage <- unique(patient_df$Stage)[1]
    }else{
      parsed_data$Stage <- NA
    }
    
    if ("Gender" %in% colnames(meta)){
      parsed_data$Gender <- unique(patient_df$Gender)[1]
    }else{
      parsed_data$Gender <- NA
    }
    
    if(verbose){print(parsed_data)}
    
    patient_meta <- rbind(patient_meta, parsed_data)
  }
}

pivot_meta <- tibble(patient_meta) %>%
  tidyr::pivot_wider(names_from = Var1, values_from = Freq)

head(as.data.frame(pivot_meta), 100)

write_csv(pivot_meta, "Parsed_meta_4.csv")

# Back to Desktop -------------- 
metap <- read.table("~/src/RNA-seq-immune-sig/data/TISCH/Parsed_meta_3.csv", sep=",", header=TRUE)

# Get numbers
# How many datasets?
files <- paste(metap$cancer, metap$dataset, metap$other, sep="_")
length(unique(files))
#[1] 66

# How many patients?
dim(metap)
#[1] [1] 684  55

# Most frequent cell types?
cell_type_matrix <- as.matrix(metap[8:length(colnames(metap))])
# Num of patient with at least X cells 
numpat <- apply(cell_type_matrix, 2, function(x) {sum(na.omit(x) > 5)})
sort(numpat, decreasing=TRUE)

filt <- tibble(metap) %>%
  filter(!cancer %in% c("AEL", "ALL", "AML", "MM")) %>% # exclude haematopoietic cancers
  #filter(other %in% c("Smartseq2", "10X", "CellMetainfo")) %>% # only include primary cancers
  filter(complete.cases(Mono.Macro, CD4Tconv, CD8T, B, DC, NK)) %>%  # Tregs not in healthy PBMC !
  filter(Mono.Macro > 5 & CD4Tconv > 5 & CD8T > 5 & B > 5 & DC > 5 & NK > 5)

# How many cancer types remain ?
dim(filt)
unique(filt$cancer)
table(filt$Gender)

filt <- tibble(metap) %>%
  filter(!cancer %in% c("AEL", "ALL", "AML", "MM")) %>% # exclude haematopoietic cancers
  #filter(other %in% c("Smartseq2", "10X", "CellMetainfo")) %>% # only include primary cancers
  filter(complete.cases(Mono.Macro, CD4Tconv, CD8T, B)) %>%  # Tregs not in healthy PBMC !
  filter(Mono.Macro > 5 & CD4Tconv > 5 & CD8T > 5 & B > 5)

# How many cancer types remain ?
unique(filt$cancer)
table(filt$Gender)

