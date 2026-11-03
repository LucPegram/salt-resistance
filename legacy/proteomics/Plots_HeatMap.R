# Exports Heat-Map
# Written by Luc Pegram & Claude AI.
# 2026/06/09
#===========================
library(pheatmap)
library(RColorBrewer)
library(dplyr)
#===========================

rm(list=ls())

tissue = "Shoot"

#for (tissue in c("Root", "Shoot")) {
#
#input_file <- paste0(
#  "D:/Masters_PhD/Data/Proteomics/2026/AgriGO/",
#  tissue,
#  "/sig_pVal(0.1)/FC_Genotype/1.5fold_UP&DOWN/vsWholeGenome/",
#  "results.txt")

#output_dir <- paste0(
#  "D:/Masters_PhD/Data/Proteomics/2026/AgriGO/",
#  tissue,
#  "/sig_pVal(0.1)/FC_Genotype/1.5fold_UP&DOWN/vsWholeGenome/",
#  "HeatMap/")



input_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/AgriGO/",
                     tissue,
                     "/Genotype_sig_pVal(0.05)/GenotypeFC_log2(1.5)/vs_Genome/results.txt")

output_dir <- paste0("D:/Masters_PhD/Data/Proteomics/2026/AgriGO/",
                     tissue,
                     "/Genotype_sig_pVal(0.05)/GenotypeFC_log2(1.5)/vs_Genome/Heatmap/")


#input_file <- "D:/Masters_PhD/Data/Proteomics/2026/AgriGO/Root/Phospho/PhosphoPeptides_vs_DataSet/NA_4_or_less/results.txt"

#output_dir <-  "D:/Masters_PhD/Data/Proteomics/2026/AgriGO/Root/Phospho/PhosphoPeptides_vs_DataSet/NA_4_or_less/Heatmap/"


# ==============================================================================
# LOAD AND PREPARE DATA
# ==============================================================================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


df <- read.delim(input_file, sep = "\t",
                 header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

#GO = "GO:0033993"; file = "Response to Lipid"
#GO = "GO:0006952"; file = "Defense Response"
#GO = "GO:0006970"; file = "Response to Osmotic Stress"
#GO = "GO:0009651"; file = "Response to Salt Stress"
#GO = "GO:0016491"; file = "Oxidoreductase activity"
#GO = "GO:0004601"; file = "Peroxidase activity"
#GO = "GO:0030003"; file = "Cellular Cation Homeostasis"
#GO = "GO:0006813"; file = "Potassium Ion Transport"
#GO = "GO:0006810"; file = "Transport"
#GO = "GO:0005524"; file = "ATP_Binding"
#GO = "GO:0006796"; file = "Phosphate-Containing_Compound_Metabolic_Process"
#GO = "GO:0006468"; file = "Protein_Phosphorylation"
#GO = "GO:0022891"; file = "Substrate_Specific_Transmembrane_Transporter_Activity"
#GO = "GO:0006811"; file = "Ion_Transport"
#GO = "GO:0055085"; file = "Transmembrane_Transport"
#GO = "GO:0022857"; file = "Transmembrane Transporter Activity"
GO = "GO:0005215"; file = "Transporter activity"

r <- df[df$GO_acc == GO,]

entries <- trimws(unlist(strsplit(r$entries, "//")))
entries <- entries[nchar(entries) > 0]

path_df <- data.frame(AGI = entries)

mapping_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                       tissue,
                       "Tissue_summary.csv")

map_df <- read.csv(mapping_file, header = T, check.names = F, stringsAsFactors = F)

path_df <- path_df |>
  left_join(map_df, by = c("AGI" = "AGI"))

write.csv(path_df, paste0(output_dir, file, ".csv"), row.names = F, quote = TRUE)


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
  Treatment = sub("^[^_]+_([0-9]+mM)_[0-9]+$", "\\1", colnames(heat_mat)),
  row.names = colnames(heat_mat)
)

# 4. Colour palette: blue (low) → white (mid) → red (high) for the heatmap body
colour_scale <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# 5. Explicit annotation colours — Genotype and NaCl Treatment sidebars
ann_colors <- list(
  Genotype  = c("No-0" = "#4CAF50", "NIT2" = "#8E44AD"),
  Treatment = c("0mM" = "#00BFA5", "75mM" = "#F4511E")
)

# 6. Draw and export
pheatmap(
  heat_mat_scaled,
  color             = colour_scale,
  annotation_col    = sample_meta,        # genotype + treatment sidebar
  annotation_colors = ann_colors,
  cluster_rows      = TRUE,               # hierarchical clustering on proteins
  cluster_cols      = FALSE,              # preserve original column order
  show_rownames     = TRUE,               # hierarchical clustering on samples
  show_colnames     = TRUE,
  fontsize_row      = 6,                  # shrink row labels for large datasets
  fontsize_col      = 8,
  border_color      = NA,
  main              = file,
  filename          = file.path(output_dir, paste0("heatmap_", file, ".png")),
  width             = 10,
  height            = max(6, nrow(heat_mat_scaled) * 0.15), # scales with n proteins
  res               = 300
)

message("Heatmap written to: ", file.path(output_dir))

#}