# Data Processing, including NA imputation and log transformation
# Written by Luc Pegram, with assistance from Claude
# 2026/05/12


rm(list=ls())

input_file  <- "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/Phospho_filtered.csv"

output_file <- "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/Phospho_anova.csv"

df_og <- read.csv(input_file, header = T, check.names = F)

df <- df_og

treatment_groups <- list(
  treatment1 = c("No.0_0mM_1", "No.0_0mM_2", "No.0_0mM_3"),
  treatment2 = c("No.0_75mM_1", "No.0_75mM_2", "No.0_75mM_3"),
  treatment3 = c("NIT2_0mM_1", "NIT2_0mM_2", "NIT2_0mM_3"),
  treatment4 = c("NIT2_75mM_1", "NIT2_75mM_2", "NIT2_75mM_3"))

df_imputed <- df

for (grp in treatment_groups) {
  grp_cols <- intersect(grp, colnames(df_imputed))          # guard against missing cols
  grp_means <- rowMeans(df_imputed[, grp_cols], na.rm = TRUE)
  
  na_mask <- is.na(df_imputed[, grp_cols])                  # logical matrix: rows x cols
  df_imputed[, grp_cols][na_mask] <- grp_means[row(na_mask)[na_mask]]
}

#write.csv(df_imputed, output_file, row.names = F)

df[, c(12:23)] <- df_imputed[, c(12:23)]; rm(df_imputed)

#================================================================

library(car)
library(tidyr)
library(dplyr)


# Protein & Gene column name
protein_col <- "Protein.Group"


# 2. Prepare data set----------------------------------------

#sample_cols <- colnames(df)[colnames(df) != protein_col]
sample_cols <- colnames(df)[12:23]

long_df <- pivot_longer(
  df,
  cols = all_of(sample_cols),
  names_to = "sample",
  values_to = "values")

long_df <- separate_wider_delim(
  long_df,
  cols = "sample",
  names = c("genotype", "treatment", "sample"),
  delim = "_")

#head(long_df, 6)

#long_df <- left_join(long_df, metadata, by = "sample")
#long_df <- arrange(long_df, .data[[protein_col]])

proteins   <- unique(long_df[[protein_col]])  #extracts unique proteins
n_proteins <- length(proteins)  #total number of proteins
protein_list <- split(long_df, long_df[[protein_col]])  #data frame per protein


long_df$values <- as.numeric(long_df$values)
long_df$genotype  <- as.factor(long_df$genotype)
long_df$treatment <- as.factor(long_df$treatment)


# Anova of full matrix
mdl <- lm(values ~ genotype*treatment, data = long_df)
A <- Anova(mdl, type = "II"); A

r <- residuals(mdl)
car::qqPlot(r) #qqnorm(r); qqline(r)

hist(r, freq = FALSE, breaks = 20, xlim=c(min(r),max(r)))
s = seq(-60,60,1/2)
lines(dnorm(s, sd=sqrt(var(r)))~s, lty=3, lwd=2.5, col='dark red')

# Individual Anova's per protein-----------------------------
results <- data.frame(
  protein       = proteins,
  p_genotype    = NA_real_,
  p_treatment   = NA_real_,
  p_interaction = NA_real_,
  stringsAsFactors = F)


for (i in seq_len(n_proteins)) {
  
  prot_df <- protein_list[[proteins[i]]]  #data frame for one protein
  
  mdl <- lm(values ~ genotype*treatment, data = prot_df)
  
  anova_table <- Anova(mdl, type = "II")
  
  results$p_genotype[i]    <- anova_table["genotype", "Pr(>F)"]
  results$p_treatment[i]   <- anova_table["treatment", "Pr(>F)"]
  results$p_interaction[i] <- anova_table["genotype:treatment", "Pr(>F)"]
}

# 5. FDR (Benjamini-Hochberg, applied per term) correction
results$p_FDR_genotype    <- p.adjust(results$p_genotype,    method = "BH")
results$p_FDR_treatment   <- p.adjust(results$p_treatment,   method = "BH")
results$p_FDR_interaction <- p.adjust(results$p_interaction, method = "BH")

df_og <- df_og |> left_join(results, by = c("Protein.Group" = "protein"))


# Export file ---------------------------------------------

#write.csv(df_og, output_file, row.names = FALSE)

#========================================================================================
#===========================
library(pheatmap)
library(RColorBrewer)
library(dplyr)
#===========================

# ── Heatmap of relative protein abundance ────────────────────────────────────

# 1. Extract and scale the intensity matrix (columns 8–19)
heat_mat <- as.matrix(df[, 12:23])

rownames(heat_mat) <- df$Genes

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

heat_dir <- "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/"

# 5. Draw and export
pheatmap(
  heat_mat_scaled,
  color             = colour_scale,
  annotation_col    = sample_meta,        # genotype + treatment sidebar
  cluster_rows      = TRUE,               # hierarchical clustering on proteins
  cluster_cols      = T,              # preserve original column order
  show_rownames     = TRUE,               # hierarchical clustering on samples
  show_colnames     = TRUE,
  fontsize_row      = 6,                  # shrink row labels for large datasets
  fontsize_col      = 8,
  border_color      = NA,
  main              = "Heatmap of Phospho-peptides",
  filename          = file.path(heat_dir, "Phospho_Heatmap.pdf"),
  width             = 10,
  height            = max(6, nrow(heat_mat_scaled) * 0.15)  # scales with n proteins
)

message("Heatmap written to: ", file.path(heat_dir))
