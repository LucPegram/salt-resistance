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

for (test in c("ordinary_t-test", "limma_eBayes_t-test")) {
  
  input_file <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/NIT2vsNo-0_Treatments/R_analysis/",
                       test,
                       "/genotype_comparison_significant.csv")
  
  output_file = input_file
  
  
  df <- read.csv(input_file, header = T, check.names = F)
  
  
  # Maps Protein Names to AGI Names -------------------------------
  
  tair_mapping_file <- mapping_file[mapping_file$type_name == "TAIR", ]
  tair_mapping_file <- tair_mapping_file[, -2]  #removes second column
  tair_mapping_file <- distinct(tair_mapping_file, name, .keep_all = T)
  
  df$"AGIs" <- sub("\\..*", "", df$`AGI MAPPING`)  #removes gene variants
  
  df <- df |>
    left_join(tair_mapping_file, by = c("AGIs" = "name"))
  
  rm(tair_mapping_file)
  
  gene_mapping_file <- mapping_file[mapping_file$type_name == "Gene_Name", ]
  gene_mapping_file <- gene_mapping_file[, -2]  #removes second column
  gene_mapping_file <- distinct(gene_mapping_file, uniprot, .keep_all = T)
  
  df <- df |>
    left_join(gene_mapping_file, by = c("uniprot" = "uniprot"))
  
  names(df)[names(df) == "name"] <- "Gene_Name"
  
  df$AGIs <- NULL
  
  df <- df |> relocate("Gene_Name", "uniprot")
  
  write.csv(df, output_file, row.names = F) #Write filtered file.
  
}