# Adds gene descriptions
# Written by Luc Pegram, assistance from Claude AI.
# 2026/05/09
#===========================
library(dplyr)
#===========================

rm(list=ls())

# Data import and pre-processing
# ── Helper: truncate long text columns for Excel compatibility ─────────────────
# Excel's hard cell limit is 32,767 characters; 32,000 gives a safe margin.

EXCEL_CHAR_LIMIT <- 32000

truncate_col <- function(df, col_name, limit) {
  if (col_name %in% names(df)) {
    long_vals <- nchar(df[[col_name]]) > limit & !is.na(df[[col_name]])
    df[[col_name]][long_vals] <- paste0(
      substr(df[[col_name]][long_vals], 1, limit),
      "... [truncated]"
    )
    cat(sprintf("  '%s': %d value(s) truncated\n", col_name, sum(long_vals)))
  } else {
    cat(sprintf("  '%s': column not found, skipping\n", col_name))
  }
  return(df)
}


for (test in c("ordinary_t-test", "limma_eBayes_t-test")) {
  
  input_file <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/NIT2vsNo-0_Treatments/R_analysis/",
                       test,
                       "/genotype_comparison_significant.csv")
  
  output_file = input_file
  
  df <- read.csv(input_file, header = T, check.names = F)
  
  map_df <- read.delim("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/TAIR_files/TAIR10_functional_descriptions.tsv",
                       header = T,
                       sep = '\t',
                       stringsAsFactors = F)
  
  map_df <- distinct(map_df, Model_name, .keep_all = TRUE)  #removed duplicated rows.
  
  df <- df |>
    left_join(map_df, by = c("AGI MAPPING" = "Model_name"))
  
  
  # ── Truncate long text columns ───────────────────────────────────────────────
  cat(sprintf("\n[%s] Checking columns for Excel character limit (%d chars):\n",
              test, EXCEL_CHAR_LIMIT))
  
  df <- truncate_col(df, "Curator_summary",           EXCEL_CHAR_LIMIT)
  df <- truncate_col(df, "Computational_description", EXCEL_CHAR_LIMIT)
  
  write.csv(df, output_file, row.names = F)
}