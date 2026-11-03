# Adds AGI names as the first column and replaces old file.
# Written by Luc Pegram, assistance from Claude AI.
# 2026/05/09
#===========================
library(dplyr)
#===========================


rm(list=ls())


# Data import and pre-processing
#---------------------------

file_name <- "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/phospho_normalized_final_140526.csv"


df <- read.csv(file_name, header = T, check.names = F, stringsAsFactors = F)
df <- df[, -1]

#df$Protein.Group <- sub(";.*", "", df$Protein.Group)  #removes gene variants


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
