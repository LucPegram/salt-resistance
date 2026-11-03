# Adds gene descriptions
# Written by Luc Pegram.
# 2026/06/11
#===========================
library(dplyr)
#===========================

rm(list=ls())

#---------------------------
## Shoot Tissue
Shoot = T; Root = F 
#---------------------------
## Root Tissue
Root = T; Shoot = F
#---------------------------


if (Shoot) {tissue <- "Shoot"} else {tissue <- "Root"}


for (tissue in c("Shoot", "Root")) {
  
  
  input_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                        tissue,
                        "Tissue_summary.csv")
  
  output_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                        tissue,
                        "Tissue_summary_NotImputed.csv")
  
  df <- read.csv(input_file, header = T, check.names = F)
  
  val_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                         tissue,
                         "Tissue_pg_filtered.csv")
  
  val_df <- read.csv(val_file, header=T, check.names = F, stringsAsFactors = F)
  val_df <- val_df[, -c(4:7)]
  
  val_df <- df[, -c(5:ncol(df))] |>
    left_join(val_df, by = c("Protein.Names" = "Protein.Names"))
  
  df[, 14:25] <- val_df[, 7:18]
  
  write.csv(df, output_file, row.names = F, quote = TRUE) #Write filtered file.
  
}
