# Adds AGI names as the first column and replaces old file.
# Written by Luc Pegram, assistance from Claude AI.
# 2026/05/09
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

# Data import and pre-processing
#---------------------------
if (Shoot) {tissue <- "Shoot"} else {tissue <- "Root"}

file_name <- paste0("D:/Masters_PhD/Data/Proteomics/2026/ForAnalysis/",
                    tissue,
                    "Tissue_pg_filtered.csv")

file_name <- paste0("D:/Masters_PhD/Data/Proteomics/2026/ForAnalysis/R_Exports/",
                    tissue,
                    "Tissue_VolcanoPlot_Table.csv")

df <- read.csv(file_name, header = T, check.names = F)


# Maps Protein Names to AGI Names -------------------------------

mapping_file <- read.table(
  "D:/Masters_PhD/Data/Proteomics/2026/ForAnalysis/TAIR_files/TAIR2UniprotMapping.txt",
  col.names = c("uniprot", "locus", "AGI"),
  header = T,
  sep = '\t',
  stringsAsFactors = F)


mapping_file <- mapping_file[, -2]  #removes second column

mapping_file$uniprot <- sub("-.*", "", mapping_file$uniprot)  #removes gene variants

mapping_file <- distinct(mapping_file)  #removed duplicated rows.

df <- df |>
  left_join(mapping_file, by = c("Protein.Group" = "uniprot"))

df <- df |> relocate("AGI")

write.csv(df, file_name, row.names = F) #Write filtered file.
