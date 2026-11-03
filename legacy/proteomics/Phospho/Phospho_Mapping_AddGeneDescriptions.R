# Add TAIR gene function descriptions
# Luc Pegram
# 17/05/2026

library(dplyr)


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


input_file <- "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/phospho_normalized_final_140526.csv"

output_file <- "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/phospho_normalized.csv"


df <- read.csv(input_file, header = T, check.names = F, stringsAsFactors = F)

map_df <- read.delim("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/TAIR_files/TAIR10_functional_descriptions.tsv",
                     header = T,
                     sep = '\t',
                     stringsAsFactors = F)


map_df$Model_name <- sub("\\..*", "", map_df$Model_name)  #removes gene variants

map_df <- distinct(map_df, Model_name, .keep_all = TRUE)  #removed duplicated rows.

df <- df |>
  left_join(map_df, by = c("AGI" = "Model_name"))


# ── Truncate long text columns ───────────────────────────────────────────────
cat(sprintf("\nChecking columns for Excel character limit (%d chars):\n",
            EXCEL_CHAR_LIMIT))

df <- truncate_col(df, "Curator_summary",           EXCEL_CHAR_LIMIT)
df <- truncate_col(df, "Computational_description", EXCEL_CHAR_LIMIT)


write.csv(df, output_file, row.names = F, quote = TRUE) #Write filtered file.
