# Export Venn Diagram of shared GOs
# Luc Pegram
# 12/06/2026

library(VennDiagram)


rm(list=ls())


input_file1 <- paste0("D:/Masters_PhD/Data/Proteomics/2026/AgriGO/",
                      "Root",
                      "/All/results.txt")

input_file2 <- "D:/Masters_PhD/Data/Proteomics/2026/AgriGO/Root/Phospho/PhosphoPeptides_vs_Genome/results.txt"

output_dir <- "D:/Masters_PhD/Data/Proteomics/2026/AgriGO/Root/Phospho/PhosphoPeptides_vs_Genome/"

#dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df_full <- read.delim(input_file1, sep = "\t", header = T, stringsAsFactors = F, check.names = F)
df_phos <- read.delim(input_file2, sep = "\t", header = T, stringsAsFactors = F, check.names = F)

keep_full <- df_full$FDR < 0.05
keep_phos <- df_phos$FDR < 0.05

go_full <- df_full$GO_acc [keep_full]
go_phos <- df_phos$GO_acc [keep_phos]

d   <- 0  #number of enriched GOs unique to data phos
g_d <- 0  #number of enriched GOs shared

df_phos_uniq <- data.frame()

for (i in 1:length(go_phos)) {
  if (go_phos[i] %in% go_full) { g_d = g_d + 1 }
  else {
    d = d + 1
    r = which(df_phos$GO_acc == go_phos[i])
    df_phos_uniq <- rbind(df_phos_uniq, df_phos[r, ])
  }
}

g = length(go_full) - g_d  #number of enriched GOs unique to full Data Set

write.csv(df_phos_uniq, paste0(output_dir, "Unique_to_phos.csv"), row.names = F)

# Suppress log file generation by VennDiagram
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn_file <- paste0(output_dir, "Unique_to_phos_Venn.png")

png(filename = venn_file, width = 800, height = 700, res = 120)

venn.plot <- draw.pairwise.venn(
  area1          = g + g_d,          # Total full Dataset
  area2          = d + g_d,          # Total Sybset
  cross.area     = g_d,              # Shared
  category       = c("Full Dataset", "phos"),
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

