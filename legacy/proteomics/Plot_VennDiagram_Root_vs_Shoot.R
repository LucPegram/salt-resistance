# Export Venn Diagram of shared diff reg proteins between root and shoot tissue.
# Luc Pegram
# 16/07/2026

library(VennDiagram)


rm(list=ls())


root_input_file  <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                           "Root",
                           "Tissue_summary.csv")

shoot_input_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                           "Shoot",
                           "Tissue_summary.csv")

output_dir <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/VennDiagrams/")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

root_df   <- read.csv(root_input_file,  header = T, stringsAsFactors = F, check.names = F)
shoot_df  <- read.csv(shoot_input_file, header = T, stringsAsFactors = F, check.names = F)

root_keep  <- (root_df$p_FDR_genotype  < 0.05) & (abs(root_df$log.Genotype) > log2(1.5))
shoot_keep <- (shoot_df$p_FDR_genotype  < 0.05) & (abs(shoot_df$log.Genotype) > log2(1.5))

root_df  <- root_df[root_keep, ]
shoot_df <- shoot_df[shoot_keep, ]

r   <- 0  #number of root not in shoot
r_s <- 0  #number of root and shoot

for (i in 1:nrow(root_df)) {
  if (root_df$AGI[i] %in% shoot_df$AGI) { r_s = r_s + 1 }
  else { r = r + 1  }
}

s <- length(shoot_df$AGI) - r_s  #number of shoot not in root


# Suppress log file generation by VennDiagram
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn_file <- paste0(output_dir, "Venn_Root_vs_Shoot.png")

png(filename = venn_file, width = 800, height = 700, res = 120)

venn.plot <- draw.pairwise.venn(
  area1          = r + r_s,          # Total root
  area2          = s + r_s,          # Total shoot
  cross.area     = r_s,              # Shared
  category       = c("Root", "Shoot"),
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
