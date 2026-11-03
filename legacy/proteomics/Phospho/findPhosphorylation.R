#check phosphorylation percent in data

input_file <- "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/phospho_normalized_final_140526.csv"

df <- read.csv(input_file, header = T, check.names = F, stringsAsFactors = F)

phos <- grep("UniMod:21", df$Modified.Sequence)
total <- nrow(df)

print(paste0(100*phos/total, "%"))

df_phos <- df[phos,]

output_file <- "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/phospho_normalized_extracted.csv"

write.csv(df_phos, output_file, row.names = F, quote = TRUE) #Write filtered file.
