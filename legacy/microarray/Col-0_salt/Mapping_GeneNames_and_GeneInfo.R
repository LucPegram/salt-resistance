# Adds AGI names as the first column and replaces old file.
# Written by Luc Pegram, assistance from Claude AI.
# 2026/05/09
#===========================
library(dplyr)
#===========================

rm(list=ls())

mapping_file <- read.delim(
  "D:/Masters_PhD/Data/Proteomics/2026/UniprotMapping/ARATH_3702_idmapping.dat",
  col.names = c("uniprot", "type_name", "name"),
  header = T,
  sep = '\t',
  stringsAsFactors = F)

tair_mapping_file <- mapping_file[mapping_file$type_name == "TAIR", ]
tair_mapping_file <- tair_mapping_file[, -2]  #removes second column
tair_mapping_file <- distinct(tair_mapping_file, name, .keep_all = T)

gene_mapping_file <- mapping_file[mapping_file$type_name == "Gene_Name", ]
gene_mapping_file <- gene_mapping_file[, -2]  #removes second column
gene_mapping_file <- distinct(gene_mapping_file, uniprot, .keep_all = T)

rm(mapping_file)

des_file <- read.csv("D:/Masters_PhD/Data/Microarray Data_From Lara/Stu_MappingDescription.csv",
                     header = T)

for (tissue in c("Shoot", "Root")) {
  
  
  input_file <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/Col-0_SaltHydroponics/",
                       tissue,
                       "/R_output/",
                       "linear_trend_results_significant.csv")
  
#  input_file <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/Col-0_SaltHydroponics/",
#                       tissue,
#                       "/R_output/poly2/",
#                       "poly2_results_significant.csv")
  
  output_file = input_file
  
  
  df <- read.csv(input_file, header = T, check.names = F)
  
  # Maps Description to AGI Names -------------------------------
  
  df <- df |>
    left_join(des_file, by = c("AGI MAPPING" = "Gene.ID"))
  
  # Maps Protein Names to AGI Names -------------------------------
  
  df$"AGIs" <- sub("\\..*", "", df$`AGI MAPPING`)  #removes gene variants
  
  df <- df |>
    left_join(tair_mapping_file, by = c("AGIs" = "name"))
  
  df <- df |>
    left_join(gene_mapping_file, by = c("uniprot" = "uniprot"))
  
  names(df)[names(df) == "name"] <- "Gene_Name"
  
  df$AGIs <- NULL
  
  df <- df |> relocate("Gene_Name", "uniprot")
  
  write.csv(df, output_file, row.names = F) #Write filtered file.
  
}