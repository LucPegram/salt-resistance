# Export the P-val significant proteins according the Anova already run.
# Luc Pegram
# 17/05/2026


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


input_file  <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                      tissue,
                      "Tissue_summary_NotImputed.csv")

output_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                       tissue,
                       "Tissue_Sig(0.1)Summary.csv")

df <- read.csv(input_file, header = T, stringsAsFactors = F, check.names = F)

# Significance threshold for ORA
q_threshold <- 0.1


sig_rows <- df$p_FDR_genotype    < q_threshold
sig_df_g <- df[sig_rows,] 

sig_rows <- df$p_FDR_treatment   < q_threshold
sig_df_t <- df[sig_rows,] 

sig_rows <- df$p_FDR_interaction < q_threshold
sig_df_i <- df[sig_rows,]


sig_df <- rbind(sig_df_g, sig_df_t, sig_df_i)
sig_df <- dplyr::distinct(sig_df)


write.csv(sig_df, output_file, row.names = FALSE)

}