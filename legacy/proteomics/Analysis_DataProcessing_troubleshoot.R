# Data Processing, including NA imputation and log transformation
# Written by Luc Pegram, with assistance from Claude
# 2026/05/12


rm(list=ls())

#---------------------------
## Shoot Tissue
Shoot = T; Root = F 
#---------------------------
## Root Tissue
Root = T; Shoot = F
#---------------------------


if (Shoot) {tissue <- "Shoot"} else {tissue <- "Root"}


  
  
  input_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                       tissue,
                       "Tissue_pg_filtered.csv")
  output_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                        tissue,
                        "Tissue_pg_imputed&log.csv")
  
  df <- read.csv(input_file, header=T)
  
  treatment_groups <- list(
    treatment1 = c("No.0_0mM_1", "No.0_0mM_2", "No.0_0mM_3"),
    treatment2 = c("No.0_75mM_1", "No.0_75mM_2", "No.0_75mM_3"),
    treatment3 = c("NIT2_0mM_1", "NIT2_0mM_2", "NIT2_0mM_3"),
    treatment4 = c("NIT2_75mM_1", "NIT2_75mM_2", "NIT2_75mM_3"))
  
  
  # For each treatment, replace NA with the smallest value in the data set only if all three replicates are NA.
  val <- min(df[, 8:19], na.rm = T)
  for (cols in 1:length(treatment_groups)){
    t=treatment_groups[[cols]]
    all_na_rows <- rowSums(is.na(df[, t])) == length(t)
    df[all_na_rows, t] <- val
  }
  
  path_genes = c(
    "APC2",
    "PILS5",
    "At3g47420",
    "NCS1",
    "DTX36",
    "NHX5")
  
  
  
  
  val_df <- as.matrix(df[,-c(1:7)])
  sel=df[,"Genes"] %in% path_genes
  #wsel=which(sel)
  #nn=length(wsel)
  #wsel=c(wsel,1:1450)
  
  #knn <- impute::impute.knn(t(val_df[wsel,]), k=3, rowmax=0.7)  #performs imputation
  
  #data <- t(knn$data)
  #data[1:nn,]
  
  #val_df[1400:1450,]
  
  val_df_log <- log2(val_df)
  
  knn <- impute::impute.knn(t(val_df_log), k=3, rowmax=0.7)  #performs imputation
  
  data <- t(knn$data)
  
  #  #log transformation
  
  df_imputed = df
  df_imputed[,-c(1:7)] <- data  #matrix
  
  df_imputed[df_imputed$Genes=="APC2",]
  
  write.csv(df_imputed, output_file, row.names = F)
