# Calculating Log Change
# Written by Luc Pegram
# 2026/05/13



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


input_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                     tissue,
                     "Tissue_pg_imputed&log.csv")

output_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                       tissue,
                       "Tissue_pg_logFC(min).csv")

df_imputed <- read.csv(input_file, header=T, check.names = T)

#df$Protein.Group <- sub(";.*", "", df$Protein.Group)  #removes the additional protein names.
#df$Genes <- sub(";.*", "", df$Genes)  #removes the additional gene names.
#df <- df[, -c(2,4:6)]  #removes unnecessary columns.

treatment_groups <- list(
  treatment1 = c("No.0_0mM_1", "No.0_0mM_2", "No.0_0mM_3"),
  treatment2 = c("No.0_75mM_1", "No.0_75mM_2", "No.0_75mM_3"),
  treatment3 = c("NIT2_0mM_1", "NIT2_0mM_2", "NIT2_0mM_3"),
  treatment4 = c("NIT2_75mM_1", "NIT2_75mM_2", "NIT2_75mM_3"))


#Fold Change Calculation
#-------------------------------

df_notImputed <- 2^(df_imputed[, 8:19])

attach(df_notImputed)
  
fc_Salt_vs_Untreated <- rowMeans(cbind(No.0_75mM_1 / No.0_0mM_1,
                                       No.0_75mM_2 / No.0_0mM_2,
                                       No.0_75mM_3 / No.0_0mM_3,
                                       NIT2_75mM_1 / NIT2_0mM_1,
                                       NIT2_75mM_2 / NIT2_0mM_2,
                                       NIT2_75mM_3 / NIT2_0mM_3), na.rm = F)

fc_NIT2_vs_No0       <- rowMeans(cbind(NIT2_0mM_1  / No.0_0mM_1,
                                       NIT2_0mM_2  / No.0_0mM_2,
                                       NIT2_0mM_3  / No.0_0mM_3,
                                       NIT2_75mM_1 / No.0_75mM_1,
                                       NIT2_75mM_2 / No.0_75mM_2,
                                       NIT2_75mM_3 / No.0_75mM_3), na.rm = F)

fc_no0vno0_saltvUnt   <- rowMeans(cbind(No.0_75mM_1 / No.0_0mM_1,
                                        No.0_75mM_2 / No.0_0mM_2,
                                        No.0_75mM_3 / No.0_0mM_3), na.rm = F)

fc_nit2vnit2_saltvUnt <- rowMeans(cbind(NIT2_75mM_1 / NIT2_0mM_1,
                                        NIT2_75mM_2 / NIT2_0mM_2,
                                        NIT2_75mM_3 / NIT2_0mM_3), na.rm = F)

fc_nit2vno0_UntvUnt   <- rowMeans(cbind(NIT2_0mM_1 / No.0_0mM_1,
                                        NIT2_0mM_2 / No.0_0mM_2,
                                        NIT2_0mM_3 / No.0_0mM_3),  na.rm = F)

fc_nit2vno0_saltvsalt <- rowMeans(cbind(NIT2_75mM_1 / No.0_75mM_1,
                                        NIT2_75mM_2 / No.0_75mM_2,
                                        NIT2_75mM_3 / No.0_75mM_3), na.rm = F)

detach(df_notImputed)



fc_Salt_vs_Untreated  <- log2(fc_Salt_vs_Untreated)

fc_NIT2_vs_No0        <- log2(fc_NIT2_vs_No0)

fc_no0vno0_saltvUnt   <- log2(fc_no0vno0_saltvUnt)

fc_nit2vnit2_saltvUnt <- log2(fc_nit2vnit2_saltvUnt)

fc_nit2vno0_UntvUnt   <- log2(fc_nit2vno0_UntvUnt)

fc_nit2vno0_saltvsalt <- log2(fc_nit2vno0_saltvsalt)


fc_df <- dplyr::bind_cols(df_imputed[, 1:5],
                          "log.Genotype"  = fc_NIT2_vs_No0,
                          "log.Treatment" = fc_Salt_vs_Untreated,
                          "log.No.0_salt/untreated" = fc_no0vno0_saltvUnt,
                          "log.NIT2_salt/untreated" = fc_nit2vnit2_saltvUnt,
                          "log.NIT2/No.0_untreated/untreated" = fc_nit2vno0_UntvUnt,
                          "log.NIT2/No.0_salt/salt" = fc_nit2vno0_saltvsalt,
                          df_imputed[, 6:7],
                          df_notImputed)


#fc_df[fc_df == "Inf"] <- NA  #sets Inf to NA
#for (i in 1:ncol(fc_df)) {  #sets NaN to NA
#  fc_df[,i][is.nan(fc_df[,i])] <- NA
#}

names(fc_df)[14:25] = c("No-0_0mM_1",  "No-0_0mM_2",  "No-0_0mM_3",
                       "No-0_75mM_1", "No-0_75mM_2", "No-0_75mM_3",
                       "NIT2_0mM_1",  "NIT2_0mM_2",  "NIT2_0mM_3",
                       "NIT2_75mM_1", "NIT2_75mM_2", "NIT2_75mM_3")

write.csv(fc_df, output_file, row.names = F)

}
