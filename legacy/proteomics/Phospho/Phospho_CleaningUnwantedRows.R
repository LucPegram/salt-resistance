# Filter rows that have at least two non-NA values in at least one treatment group.
# Written by Claude AI. Adjusted by Luc Pegram
# 2026/05/09


rm(list=ls())

input_file <- "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/phospho_normalized_extracted.csv"

output_file <- "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/Phospho_filtered.csv"

df <- read.csv(input_file, header=T, stringsAsFactors = F, check.names = F)

na_threshold = 3

# Define the 4 treatment groups by column names.
treatments <- df[, 12:23]

keep <- rowSums(is.na(df)) <= na_threshold

df_filtered <- df[keep, ]

cat("Rows kept:   ", nrow(df_filtered), "\n")
cat("Rows removed:", nrow(df)-nrow(df_filtered), "\n")

write.csv(df_filtered, output_file, row.names = F) #Write filtered file.
