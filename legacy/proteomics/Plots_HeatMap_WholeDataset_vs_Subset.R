# Exports Heat-Map
# Written by Luc Pegram & Claude AI.
# 2026/06/09
#===========================
library(pheatmap)
library(RColorBrewer)
library(dplyr)
#===========================


rm(list=ls())

#===========================
## Shoot Tissue
Shoot = T; Root = F 
#===========================
## Root Tissue
Root = T;  Shoot = F
#===========================

if (Shoot) {tissue <- "Shoot"} else {tissue <- "Root"}

input_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/VennDiagrams/",
                     tissue,
                     "/AGI_Comparison/Unique_to_Subset.csv")

output_dir <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/VennDiagrams/",
                     tissue,
                     "/AGI_Comparison/Heatmap/")


# ==============================================================================
# LOAD AND PREPARE DATA
# ==============================================================================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)


for (i in 1:nrow(df)) {
  
  go   = df$GO_acc[i]
  file = df$Term[i]
  #file <- gsub(",", "", file)
  if (nchar(file) > 153) { file = strsplit(file, " ")[[1]][1]}
  
  r <- df[df$GO_acc == go,]
  
  entries <- trimws(unlist(strsplit(r$entries, "//")))
  entries <- entries[nchar(entries) > 0]
  
  path_df <- data.frame(AGI = entries)
  
  mapping_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                         tissue,
                         "Tissue_summary.csv")
  
  NA_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                         tissue,
                         "Tissue_summary_NotImputed.csv")
  
  map_df <- read.csv(mapping_file, header = T, check.names = F, stringsAsFactors = F)
  df_NA <- read.csv(NA_file, header = T, check.names = F, stringsAsFactors = F)
  
  path_df <- path_df |>
    left_join(map_df, by = c("AGI" = "AGI"))
  
  path_df_NA <- data.frame(AGI = entries)
  path_df_NA <- path_df_NA |>
    left_join(df_NA, by = c("AGI" = "AGI"))
  
  write.csv(path_df_NA, paste0(output_dir, file, ".csv"), row.names = F, quote = TRUE)
  
  rm(path_df_NA)
  
  # ── Heatmap of relative protein abundance ────────────────────────────────────
  
  # 1. Extract and scale the intensity matrix (columns 8–19)
  heat_mat <- as.matrix(path_df[, 14:25])
  heat_mat = log2(heat_mat)
  
  rownames(heat_mat) <- path_df$Genes
  
  # 2. Z-score scale across rows (per-protein relative abundance)
  heat_mat_scaled <- t(scale(t(heat_mat)))
  
  # 3. Build a clean column annotation from sample names
  sample_meta <- data.frame(
    Genotype  = sub("_[0-9]+mM_[0-9]+$", "", colnames(heat_mat)),
    Treatment = paste0(sub("^[^_]+_([0-9]+mM)_[0-9]+$", "\\1", colnames(heat_mat))),
    row.names = colnames(heat_mat)
  )
  
  # 4. Colour palette: blue (low) → white (mid) → red (high)
  colour_scale <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  
  
  # 5. Draw and export
  pheatmap(
    heat_mat_scaled,
    color             = colour_scale,
    annotation_col    = sample_meta,        # genotype + treatment sidebar
    cluster_rows      = TRUE,               # hierarchical clustering on proteins
    cluster_cols      = FALSE,              # preserve original column order
    show_rownames     = TRUE,               # hierarchical clustering on samples
    show_colnames     = TRUE,
    fontsize_row      = 6,                  # shrink row labels for large datasets
    fontsize_col      = 8,
    border_color      = NA,
    main              = file,
    filename          = file.path(output_dir, paste0("heatmap_", file, ".pdf")),
    width             = 10,
    height            = max(6, nrow(heat_mat_scaled) * 0.15)  # scales with n proteins
  )
  
  message("Heatmap written to: ", file.path(output_dir))
  
}
