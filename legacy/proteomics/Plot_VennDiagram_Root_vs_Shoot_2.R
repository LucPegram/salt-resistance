# Export Venn Diagram of shared proteins between GenotypeSignificant(0.05), |FC| > 1.5,
# between Root and Shoot tissue.
# Luc Pegram
# 10/08/2026

library(VennDiagram)
library(grid)

rm(list=ls())

root_input_file   <- "D:/Masters_PhD/Data/Proteomics/2026/AgriGO/Root/Genotype_sig_pVal(0.05)/GenotypeFC_log2(1.5)/vs_Genome/AGIs.txt"
shoot_input_file  <- "D:/Masters_PhD/Data/Proteomics/2026/AgriGO/Shoot/Genotype_sig_pVal(0.05)/GenotypeFC_log2(1.5)/vs_Genome/AGIs.txt"


output_dir <- "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/VennDiagrams/"

root_df  <- read.delim(root_input_file, header = F, sep = '\t')
shoot_df <- read.delim(shoot_input_file, header = F, sep = '\t')

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

head(root_df)
head(shoot_df)

# Pull the AGI vectors out of each df, deduplicated
root_genes  <- unique(root_df$V1)
shoot_genes <- unique(shoot_df$V1)

root_only  <- length(setdiff(root_genes, shoot_genes))
shoot_only <- length(setdiff(shoot_genes, root_genes))
shared     <- length(intersect(root_genes, shoot_genes))

venn_out <- paste0(output_dir, "root_vs_shoot_GenSig(0.05)_FC(1.5).png")

png(venn_out, width = 800, height = 700, res = 120)

venn.plot <- draw.pairwise.venn(
  area1          = root_only + shared,   # Total Root
  area2          = shoot_only + shared,  # Total Shoot
  cross.area     = shared,               # Shared
  category       = c("Root", "Shoot"),
  fill           = c("#4E84C4", "#4CAF50"),
  alpha          = 0.5,
  lwd            = 4,
  lty            = "solid",
  col            = c("#2C4F8A", "#2E7D32"),
  cex            = 2.5,              # Numbers font size
  cat.cex        = 2,                # Category label font size
  cat.col        = c("#2C4F8A", "#2E7D32"),
  cat.fontface   = "bold",
  cat.dist       = 0.06,
  cat.pos        = c(-20, 20),
  margin         = 0.08,
  ind            = FALSE
)

grid.draw(venn.plot)
dev.off()

cat("Venn diagram saved to:", venn_out, "\n")
cat("Root n =", root_only + shared,
    "| Shoot n =", shoot_only + shared,
    "| Overlap n =", shared, "\n")
