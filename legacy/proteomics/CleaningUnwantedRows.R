# Filter rows that have at least two non-NA values in at least one treatment group.
# Written by Claude AI. Adjusted by Luc Pegram
# 2026/05/09


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


for (tissue in c("Shoot", "Root")) {


input_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/",
                     tissue,
                     "Tissue_pg.tsv")

output_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                       tissue,
                       "Tissue_pg_filtered.csv")

df <- read.delim(input_file, header=T, sep="\t", check.names = F)
df$Protein.Group <- sub(";.*", "", df$Protein.Group)  #removes additional protein names.


# Define the 4 treatment groups by column names.
treatment_groups <- list(
  treatment1 = c("No-0_0mM_1", "No-0_0mM_2", "No-0_0mM_3"),
  treatment2 = c("No-0_75mM_1", "No-0_75mM_2", "No-0_75mM_3"),
  treatment3 = c("NIT2_0mM_1", "NIT2_0mM_2", "NIT2_0mM_3"),
  treatment4 = c("NIT2_75mM_1", "NIT2_75mM_2", "NIT2_75mM_3"))

# Filter. Returns TRUE if at least one treatment group has >= min_val non-NA values
row_passes <- function(row, groups, min_val) {
  any(sapply(groups, function(cols) {
    sum(!is.na(row[cols])) >= min_val
  }))
}

keep <- apply(df, 1, row_passes, groups = treatment_groups, min_val = 2)

filtered_df <- df[keep, ]

sgl_pep <- filtered_df$N.Sequences > 1

filtered_df <- filtered_df[sgl_pep, ]

cat("Rows kept:   ", nrow(filtered_df), "\n")
cat("Rows removed:", nrow(df)-nrow(filtered_df), "\n")


write.csv(filtered_df, output_file, row.names = F) #Write filtered file.

}
