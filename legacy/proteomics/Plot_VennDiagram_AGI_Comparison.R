# Export Venn Diagram of shared GOs
# Luc Pegram
# 12/06/2026

library(VennDiagram)


rm(list=ls())

#===========================
## Shoot Tissue
Shoot = T; Root = F 
#===========================
## Root Tissue
Root = T;  Shoot = F
#===========================


if (Shoot) {tissue <- "Shoot"} else {tissue <- "Root"}


input_file1  <- paste0("D:/Masters_PhD/Data/Proteomics/2026/AgriGO/",
                       tissue,
                       "/All/results.txt")

input_file2  <- paste0("D:/Masters_PhD/Data/Proteomics/2026/AgriGO/",
                       tissue,
                       "/sig_pVal(0.1)/FC_Genotype/1.5fold_UP&DOWN/vsWholeGenome/results.txt")

output_dir <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/VennDiagrams/",
                     tissue,
                     "/AGI_Comparison/")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df_full   <- read.delim(input_file1, sep = "\t", header = T, stringsAsFactors = F, check.names = F)
df_subset <- read.delim(input_file2, sep = "\t", header = T, stringsAsFactors = F, check.names = F)

keep_full   <- df_full$FDR   < 0.05
keep_subset <- df_subset$FDR < 0.05

go_full  <- df_full$GO_acc  [keep_full]
go_subset <- df_subset$GO_acc [keep_subset]

d   <- 0  #number of enriched GOs unique to data Subset
g_d <- 0  #number of enriched GOs shared

df_subset_uniq <- data.frame()

for (i in 1:length(go_subset)) {
  if (go_subset[i] %in% go_full) { g_d = g_d + 1 }
  else {
    d = d + 1
    r = which(df_subset$GO_acc == go_subset[i])
    df_subset_uniq <- rbind(df_subset_uniq, df_subset[r, ])
    }
}

g = length(go_full) - g_d  #number of enriched GOs unique to full Data Set

write.csv(df_subset_uniq, paste0(output_dir, "Unique_to_Subset.csv"), row.names = F)

# Suppress log file generation by VennDiagram
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn_file <- paste0(output_dir, tissue, "_fullDataset_vs_DataSubset_Venn.png")

png(filename = venn_file, width = 800, height = 700, res = 120)

venn.plot <- draw.pairwise.venn(
  area1          = g + g_d,          # Total full Dataset
  area2          = d + g_d,          # Total Sybset
  cross.area     = g_d,              # Shared
  category       = c("Full Dataset", "Subset"),
  fill           = c("#4E84C4", "#D95F02"),
  alpha          = 0.5,
  lwd            = 4,
  lty            = "solid",
  col            = c("#2C4F8A", "#8B3A00"),
  cex            = 2.5,              # Numbers font size
  cat.cex        = 2,              # Category label font size
  cat.col        = c("#2C4F8A", "#8B3A00"),
  cat.fontface   = "bold",
  cat.dist       = 0.06,
  cat.pos        = c(-20, 20),
  margin         = 0.08,
  ind            = FALSE
)

grid.draw(venn.plot)

dev.off()

cat("Venn Diagram saved to:", venn_file, "\n")

