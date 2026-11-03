# Performing two-way anova analysis
# Written by Luc Pegram, with assistance from Claude
# 2026/05/13

library(car)
library(tidyr)
library(dplyr)


rm(list=ls())


#===========================
## Shoot Tissue
Shoot = T; Root = F 
#===========================
## Root Tissue
Root = T;  Shoot = F
#===========================

if (Shoot) {tissue <- "Shoot"} else {tissue <- "Root"}


for (tissue in c("Shoot", "Root")) {


input_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                     tissue,
                     "Tissue_pg_logFC(min).csv")

output_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                      tissue,
                      "Tissue_AnovaPvalues.csv")

og_df <- read.csv(input_file, header = T, check.names = F)

df <- og_df
df[, c(14:25)] <- log(og_df[ , 14:25])  #log transformation

# Protein & Gene column name
protein_col <- "Protein.Group"


# 2. Prepare data set----------------------------------------

#sample_cols <- colnames(df)[colnames(df) != protein_col]
sample_cols <- colnames(df)[14:25]

long_df <- pivot_longer(
  df,
  cols = all_of(sample_cols),
  names_to = "sample",
  values_to = "values")

long_df <- separate_wider_delim(
  long_df,
  cols = "sample",
  names = c("genotype", "treatment", "sample"),
  delim = "_")

#head(long_df, 6)
write.csv(long_df, paste0(
  "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/", tissue, "Tissue_long.csv"),
  row.names = FALSE)  #saves long format for use later if needed.

#long_df <- left_join(long_df, metadata, by = "sample")
#long_df <- arrange(long_df, .data[[protein_col]])

proteins   <- unique(long_df[[protein_col]])  #extracts unique proteins
n_proteins <- length(proteins)  #total number of proteins
protein_list <- split(long_df, long_df[[protein_col]])  #data frame per protein


long_df$values <- as.numeric(long_df$values)
long_df$genotype  <- as.factor(long_df$genotype)
long_df$treatment <- as.factor(long_df$treatment)


# Anova of full matrix
mdl <- lm(values ~ genotype*treatment, data = long_df)
A <- Anova(mdl, type = "II"); A


# Individual Anova's per protein-----------------------------
results <- data.frame(
  protein       = proteins,
  p_genotype    = NA_real_,
  p_treatment   = NA_real_,
  p_interaction = NA_real_,
  stringsAsFactors = F)


for (i in seq_len(n_proteins)) {
  
  prot_df <- protein_list[[proteins[i]]]  #data frame for one protein
  
  mdl <- lm(values ~ genotype*treatment, data = prot_df)
  
  anova_table <- Anova(mdl, type = "II")
  
  results$p_genotype[i]    <- anova_table["genotype", "Pr(>F)"]
  results$p_treatment[i]   <- anova_table["treatment", "Pr(>F)"]
  results$p_interaction[i] <- anova_table["genotype:treatment", "Pr(>F)"]
}

# 5. FDR (Benjamini-Hochberg, applied per term) correction
results$p_FDR_genotype    <- p.adjust(results$p_genotype,    method = "BH")
results$p_FDR_treatment   <- p.adjust(results$p_treatment,   method = "BH")
results$p_FDR_interaction <- p.adjust(results$p_interaction, method = "BH")

og_df <- og_df |> left_join(results, by = c("Protein.Group" = "protein"))


# Export file ---------------------------------------------

write.csv(og_df, output_file, row.names = FALSE)

}
