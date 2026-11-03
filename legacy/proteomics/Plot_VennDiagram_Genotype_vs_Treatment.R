# Export Venn Diagram of shared proteins between genotype and treatment effects
# Luc Pegram
# 04/06/2026

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



for (tissue in c("Shoot", "Root")) {
  
  
  input_file  <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                        tissue,
                        "Tissue_Sig(0.1)Summary.csv")
  
  output_dir <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/VennDiagrams/",
                       tissue,
                       "/")
  
  df <- read.csv(input_file, header = T, stringsAsFactors = F, check.names = F)
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Significance threshold for ORA
  q_threshold <- 0.05

  
  g   <- 0  #sig prot unique to genotype
  t   <- 0  #sig prot unique to treatment
  g_t <- 0  #sig prot shared
  
  attach(df)
  
  for (i in 1:nrow(df)) {
    
    if ((p_FDR_genotype[i] <= q_threshold) && (p_FDR_treatment[i] > q_threshold)) {
      g = g+1
    }
    else if ((p_FDR_genotype[i] > q_threshold) && (p_FDR_treatment[i] <= q_threshold)) {
      t = t+1
    }
    else if ((p_FDR_genotype[i] <= q_threshold) && (p_FDR_treatment[i] <= q_threshold)) {
      g_t = g_t+1
    }
  }  
  
  detach(df)

  # Suppress log file generation by VennDiagram
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  venn_file <- paste0(output_dir, tissue, "_Genotype_vs_Treatment_Venn(0.05).png")
  
  png(filename = venn_file, width = 800, height = 700, res = 120)
  
  venn.plot <- draw.pairwise.venn(
    area1          = g + g_t,          # Total Genotype significant
    area2          = t + g_t,          # Total Treatment significant
    cross.area     = g_t,              # Shared
    category       = c("Genotype", "Treatment"),
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
}
