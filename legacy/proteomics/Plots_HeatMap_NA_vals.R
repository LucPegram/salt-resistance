# Exports Heat-Map
# Written by Luc Pegram & Claude AI.
# 2026/06/09
#===========================
library(pheatmap)
library(RColorBrewer)
library(dplyr)
#===========================

rm(list=ls())

tissue = "Root"

input_file <- paste0(
  "D:/Masters_PhD/Data/Proteomics/2026/AgriGO/",
  tissue,
  "/sig_pVal(0.1)/FC_Genotype/1.5fold_UP&DOWN/vsWholeGenome/",
  "results.txt")

output_dir <- paste0(
  "D:/Masters_PhD/Data/Proteomics/2026/AgriGO/",
  tissue,
  "/sig_pVal(0.1)/FC_Genotype/1.5fold_UP&DOWN/vsWholeGenome/",
  "HeatMap")


# ==============================================================================
# LOAD AND PREPARE DATA
# ==============================================================================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


df <- read.delim(input_file, sep = "\t",
                 header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)


#agi = "GO:0022891"; file = "Substrate_Specific_Transmembrane_Transporter_Activity"
#agi = "GO:0006811"; file = "Ion_Transport"
agi = "GO:0055085"; file = "Transmembrane_Transport"

r <- df[df$GO_acc == agi,]

entries <- trimws(unlist(strsplit(r$entries, "//")))
entries <- entries[nchar(entries) > 0]

path_df <- data.frame(AGI = entries)

mapping_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                       tissue,
                       "Tissue_pg_filtered.csv")

map_df <- read.csv(mapping_file, header=T, check.names = F, stringsAsFactors = F)
map_df$Genes <- sub(";.*", "", map_df$Genes)  #removes additional gene names.


path_df <- path_df |>
  left_join(map_df, by = c("AGI" = "AGI"))

# ── Heatmap of relative protein abundance ────────────────────────────────────

# 1. Extract and scale the intensity matrix (columns 8–19)
heat_mat <- as.matrix(path_df[, 8:19])
heat_mat = log2(heat_mat)

rownames(heat_mat) <- path_df$Genes

# Keep rows with NA in no more than X samples (adjust min_present to suit)
min_present <- 3                          # e.g. present in at least 3 of 12 samples
keep        <- rowSums(!is.na(heat_mat)) >= min_present
heat_mat    <- heat_mat[keep, ]

# Impute remaining NAs with row mean before scaling
heat_mat <- t(apply(heat_mat, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

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
  main              = "Relative Protein Abundance",
  filename          = file.path(output_dir, paste0("heatmap_", file, "_NA.pdf")),
  width             = 10,
  height            = max(6, nrow(heat_mat_scaled) * 0.15)  # scales with n proteins
)


message("Heatmap written to: ", file.path(output_dir))
