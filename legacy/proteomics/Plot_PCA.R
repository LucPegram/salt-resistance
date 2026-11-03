# Proteomics PCA Analysis
# Input: protein name column + 12 sample columns
#        (2 genotypes x 2 treatments x 3 replicates)
#        Data assumed to be imputed and log-transformed
# Written by Claude, refined by Luc Pegram
# 2-026/05/02
# -----------------------------------------------------------------------------

library(ggplot2)
library(ggrepel)

rm(list=ls())

#---------------------------
## Shoot Tissue
Shoot = T; Root = F
#---------------------------
## Root Tissue
Root = T; Shoot = F
#---------------------------

if (Shoot) {tissue <- "Shoot"} else {tissue <- "Root"}

# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------
# Replace with your actual file path and adjust sep/header as needed
df <- read.csv(paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                      tissue,
                      "Tissue_pg_imputed&log.csv"),
               header = T, check.names = F)

# Adjust this to match your actual protein column name
protein_col <- "Protein.Group"
gene_col <- "Genes"
# -----------------------------------------------------------------------------
# 2. DEFINE SAMPLE METADATA
# -----------------------------------------------------------------------------
sample_cols <- colnames(df)[8:19]  #extracts sample column names

# Build a metadata table from your column names
# Edit genotype/treatment labels to match your experiment
metadata <- data.frame(
  sample    = sample_cols,
  genotype  = rep(c("No-0", "NIT2"), each = 6),
  treatment = rep(rep(c("untreated", "salt"), each = 3), times = 2),
  stringsAsFactors = FALSE)

metadata$group <- paste(metadata$genotype, metadata$treatment, sep = "_")

# -----------------------------------------------------------------------------
# 3. PREPARE MATRIX FOR PCA
# Rows = samples, columns = proteins (prcomp expects this orientation)
# Remove any proteins with zero variance (can't be scaled)
# -----------------------------------------------------------------------------
expr_matrix <- t(df[, sample_cols])          # samples x proteins
colnames(expr_matrix) <- df[[protein_col]]

# -----------------------------------------------------------------------------
# 4. RUN PCA
# center = TRUE, scale. = TRUE: recommended for proteomics to account for
# proteins measured on different abundance scales
# -----------------------------------------------------------------------------
pca_result <- prcomp(expr_matrix, center = TRUE, scale. = TRUE)

# Extract variance explained
pca_var        <- pca_result$sdev^2
pct_var        <- round(100 * pca_var / sum(pca_var), 1)  #<-- check and understand**
pc1_label      <- paste0("PC1 (", pct_var[1], "%)")
pc2_label      <- paste0("PC2 (", pct_var[2], "%)")

# Combine scores with metadata
scores <- as.data.frame(pca_result$x[, 1:3])
scores <- cbind(metadata, scores)  #<-- check order

# -----------------------------------------------------------------------------
# 5. DEFINE AESTHETICS
# Edit colours to suit your preferences
# -----------------------------------------------------------------------------
group_colours <- c("No-0_untreated" = "#2166AC",
                   "No-0_salt" = "#92C5DE",
                   "NIT2_untreated" = "#D6604D",
                   "NIT2_salt" = "#F4A582")

group_shapes <- c("No-0_untreated" = 16,
                  "No-0_salt" = 17,
                  "NIT2_untreated" = 16,
                  "NIT2_salt" = 17)

# -----------------------------------------------------------------------------
# 6. PLOT: PC1 vs PC2
# -----------------------------------------------------------------------------
p_pca <- ggplot(scores, aes(x = PC1, y = PC2, colour = group, shape = group)) +
  geom_polygon(
    data = do.call(rbind, lapply(split(scores, scores$group), function(g) {
      g[chull(g$PC1, g$PC2), ]})),
    aes(x = PC1, y = PC2, fill = group, group = group),
    alpha = 0.12, linetype = "dashed", linewidth = 0.5,
    show.legend = FALSE) +
  
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE,
                  box.padding = 0.4, max.overlaps = 20) +
  scale_colour_manual(values = group_colours, name = "Group") +
  scale_fill_manual(values   = group_colours) +
  scale_shape_manual(values  = group_shapes,  name = "Group") +
  labs(
    title    = "PCA of Proteomics Data",
    subtitle = paste0(tissue, " Tissue"),
    x        = pc1_label,
    y        = pc2_label) +
  
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"))


print(p_pca)

# -----------------------------------------------------------------------------
# 7. SCREE PLOT — variance explained per PC
# -----------------------------------------------------------------------------
scree_df <- data.frame(
  PC             = factor(paste0("PC", seq_along(pct_var)),levels = paste0("PC", seq_along(pct_var))),
  pct_variance   = pct_var,
  cumulative_var = cumsum(pct_var))


# Show only first 10 PCs
scree_df_sub <- head(scree_df, 10)

p_scree <- ggplot(scree_df_sub, aes(x = PC)) +
  geom_col(aes(y = pct_variance), fill = "#4393C3", alpha = 0.85) +
  geom_line(aes(y = cumulative_var, group = 1),
            colour = "#D6604D", linewidth = 1) +
  geom_point(aes(y = cumulative_var),
             colour = "#D6604D", size = 3) +
  scale_y_continuous(
    name     = "Variance Explained (%)",
    sec.axis = sec_axis(~., name = "Cumulative Variance (%)")) +
  
  labs(
    title = "Scree Plot",
    x     = "Principal Component") +
  
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"))


print(p_scree)

# -----------------------------------------------------------------------------
# 8. OPTIONAL: PC1 vs PC3
# -----------------------------------------------------------------------------
pc3_label <- paste0("PC3 (", pct_var[3], "%)")


p_pc13 <- ggplot(scores, aes(x = PC1, y = PC3,
                             colour = group, shape = group)) +
  geom_polygon(
    data = do.call(rbind, lapply(split(scores, scores$group), function(g) {
      g[chull(g$PC1, g$PC3), ]
    })),
    aes(x = PC1, y = PC3, fill = group, group = group),
    alpha = 0.12, linetype = "dashed", linewidth = 0.5,
    show.legend = FALSE
  ) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE,
                  box.padding = 0.4, max.overlaps = 20) +
  scale_colour_manual(values = group_colours, name = "Group") +
  scale_fill_manual(values   = group_colours) +
  scale_shape_manual(values  = group_shapes,  name = "Group") +
  labs(
    title    = "PCA of Proteomics Data (PC1 vs PC3)",
    subtitle = paste0(tissue, " Tissue"),
    x        = pc1_label,
    y        = pc3_label
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

print(p_pc13)

# -----------------------------------------------------------------------------
# 9. SAVE PLOTS
# -----------------------------------------------------------------------------
ggsave(paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/PCA/",
              tissue,
              "Tissue_pc1-pc2.png"),
       plot = p_pca,   width = 7, height = 5.5, dpi = 300)
ggsave(paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/PCA/",
              tissue,
              "Tissue_pca_scree.png"),
       plot = p_scree, width = 6, height = 4, dpi = 300)
ggsave(paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/PCA/",
              tissue,
              "Tissue_pc1-pc3.png"),
       plot = p_pc13,  width = 7, height = 5.5, dpi = 300)

# -----------------------------------------------------------------------------
# 10. PRINT SUMMARY
# -----------------------------------------------------------------------------
cat("\n--- Top 5 PCs: Variance Explained ---\n")
print(head(scree_df, 5))

cat("\n--- PC Loadings: Top 10 proteins driving PC1 ---\n")
loadings_pc1 <- sort(abs(pca_result$rotation[, "PC1"]), decreasing = TRUE)
print(head(loadings_pc1, 10))
