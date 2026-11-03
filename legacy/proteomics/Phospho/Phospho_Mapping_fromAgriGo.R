# Adds AGI names as the first column and replaces old file.
# Written by Luc Pegram, assistance from Claude AI.
# 2026/07/30
#===========================
library(dplyr)
#===========================


rm(list=ls())


#=====================================================


file_name <- "D:/Masters_PhD/Data/Proteomics/2026/AgriGO/Root/Phospho/PhosphoPeptides_vs_DataSet/NA_4_or_less/"


df <- read.delim(
  paste0(file_name, "results.txt"),
  sep = '\t', header = T, check.names = F, stringsAsFactors = F)


# Maps Protein Names to AGI Names -------------------------------

mapping_file <- read.csv("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/phospho_normalized_extracted.csv",
                         header = T, check.names = F, stringsAsFactors = F)

text <- "Transporter Activity"
go <- 'GO:0005215'

text <- "Response to Lipid"
go <- 'GO:0033993'

text <- "Transport"
go <- 'GO:0006810'

row <- which(df$GO_acc == go)

agi <- df$entries[row]

agi_ids <- trimws(strsplit(agi, "//")[[1]])
agi_ids <- agi_ids[agi_ids != ""]

# Match each agi_id to its row in mapping_file, keeping all columns
matched_rows <- mapping_file[match(agi_ids, mapping_file$AGI), ]

agi_df <- data.frame(
  AGI = agi_ids,
  matched_rows,
  row.names = NULL,
  stringsAsFactors = FALSE
)


write.csv(agi_df, paste0(file_name, text, ".csv"), row.names = F) #Write filtered file.
          