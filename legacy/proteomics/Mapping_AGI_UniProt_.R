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

if (Shoot) {tissue <- "Shoot"} else {tissue <- "Root"}


#=====================================================

for (tissue in c("Shoot", "Root")) {


file_name <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                    tissue,
                    "Tissue_pg_filtered.csv")


df <- read.csv(file_name, header = T, check.names = F)


# Maps Protein Names to AGI Names -------------------------------

mapping_file <- read.delim(
  "D:/Masters_PhD/Data/Proteomics/2026/UniprotMapping/ARATH_3702_idmapping.dat",
  col.names = c("uniprot", "type_name", "AGI"),
  header = T,
  sep = '\t',
  stringsAsFactors = F)


mapping_file <- mapping_file[mapping_file$type_name == "TAIR", ]
mapping_file <- mapping_file[, -2]  #removes second column
mapping_file <- distinct(mapping_file, uniprot, .keep_all = T)

df <- df |>
  left_join(mapping_file, by = c("Protein.Group" = "uniprot"))

df <- df |> relocate("AGI")

write.csv(df, file_name, row.names = F) #Write filtered file.

}
