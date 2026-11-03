## ============================================================
## ORA Visualisations — Proteomic Data
## Five plot types:
##   1. Dot Plot (Bubble Plot)
##   2. Bar Chart
##   3. Volcano-style Scatter Plot
##   4. Heatmap
##   5. Network Enrichment Map
## ============================================================

# --- Required packages --------------------------------------------------------
# Install any missing packages before running:

library(ggplot2)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(scales)
library(stringr)
library(tidyr)

# ==============================================================================
# USER INPUT — Set these before running
# ==============================================================================

#===========================
## Shoot Tissue
Shoot = T; Root = F 
#===========================
## Root Tissue
Root = T;  Shoot = F
#===========================


if (Shoot) {tissue <- "Shoot"} else {tissue <- "Root"}

## File Paths
############################

input_file <- paste0(
  "D:/Masters_PhD/Data/Proteomics/2026/ForAnalysis/R_Exports/ORA/",
  tissue, "/Genotype/", tissue, "_ORA_Map.csv"
)

output_dir <- paste0(
  "D:/Masters_PhD/Data/Proteomics/2026/ForAnalysis/R_Exports/ORA/",
  tissue, "/Plots2/"
)

# Path to the ORA output .csv file you want to visualise.
# This can be a single-category file (e.g. Shoot_ORA_BP.csv) or
# the combined AllCategories file (Shoot_ORA_AllCategories.csv).
input_file <- "D:/Masters_PhD/Data/Proteomics/2026/ForAnalysis/R_Exports/ORA/Shoot/Genotype/Shoot_ORA_AllCategories.csv"

# Directory where all plots will be saved
output_dir <- "D:/Masters_PhD/Data/Proteomics/2026/ForAnalysis/R_Exports/ORA/Shoot/Genotype/Plots2/"

# Maximum number of gene sets to show per plot / per category facet
top_n <- 15

# Significance threshold for colouring / filtering
fdr_threshold <- 0.1

# Plot dimensions (inches)
plot_width  <- 12
plot_height <- 8

# ==============================================================================
# LOAD AND PREPARE DATA
# ==============================================================================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ora <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE,
                check.names = FALSE)

# Detect whether a Category column is present (combined file) or not
has_category <- "Category" %in% names(ora)

# Coerce numeric columns (in case they were read as character)
ora$p_value    <- as.numeric(ora$p_value)
ora$p_FDR      <- as.numeric(ora$p_FDR)
ora$Odds_Ratio <- as.numeric(ora$Odds_Ratio)
ora$N_Overlap  <- as.numeric(ora$N_Overlap)
ora$N_GeneSet  <- as.numeric(ora$N_GeneSet)
ora$N_Sig      <- as.numeric(ora$N_Sig)

# Remove rows with NA p-values (gene sets absent from background)
ora <- ora[!is.na(ora$p_value), ]

# Add derived columns used across multiple plots
ora$neg_log10_FDR <- -log10(ora$p_FDR)
ora$Gene_Ratio     <- ora$N_Overlap / ora$N_Sig        # fraction of sig proteins in set
ora$Significant    <- ora$p_FDR <= fdr_threshold

# Wrap long gene set names for axis legibility (max 50 chars per line)
ora$Gene_Set_Wrapped <- str_wrap(ora$Gene_Set, width = 50)

# --- Helper: select top N rows per category by p_FDR -------------------------
top_per_cat <- function(df, n = top_n) {
  if (has_category) {
    df %>%
      group_by(Category) %>%
      slice_min(order_by = p_FDR, n = n, with_ties = FALSE) %>%
      ungroup()
  } else {
    df %>% slice_min(order_by = p_FDR, n = n, with_ties = FALSE)
  }
}

# Shared colour scale for FDR
fdr_colours <- scale_colour_gradientn(
  colours = c("#d73027", "#fc8d59", "#fee090", "#e0f3f8", "#91bfdb", "#4575b4"),
  name    = expression(-log[10](FDR)),
  trans   = "identity"
)

cat("Data loaded:", nrow(ora), "gene sets\n")
if (has_category) cat("Categories:", paste(unique(ora$Category), collapse = ", "), "\n")
cat("\n")

# ==============================================================================
# SECTION 1: DOT PLOT (BUBBLE PLOT)
# ==============================================================================
# x-axis  = Odds Ratio (enrichment strength)
# y-axis  = Gene Set name (top N per category)
# size    = N_Overlap (number of significant proteins in the set)
# colour  = -log10(p_FDR)
# ==============================================================================

cat("Plotting Section 1: Dot Plot...\n")

dot_data <- top_per_cat(ora)

# Reorder gene sets within each category by Odds Ratio
dot_data <- dot_data %>%
  arrange(if (has_category) Category else NULL, Odds_Ratio) %>%
  mutate(Gene_Set_Wrapped = factor(Gene_Set_Wrapped,
                                   levels = unique(Gene_Set_Wrapped)))

dot_plot <- ggplot(dot_data,
                   aes(x     = Odds_Ratio,
                       y     = Gene_Set_Wrapped,
                       size  = N_Overlap,
                       colour = neg_log10_FDR)) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(name   = "Overlap\n(proteins)",
                        range  = c(3, 12),
                        breaks = pretty_breaks(n = 4)) +
  scale_colour_gradientn(
    colours = c("#4575b4", "#91bfdb", "#fee090", "#fc8d59", "#d73027"),
    name    = expression(-log[10](FDR))
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  labs(
    title    = "ORA Dot Plot — Top Gene Sets by Category",
    subtitle = paste0("Top ", top_n, " gene sets per category  |  FDR threshold: ", fdr_threshold),
    x        = "Odds Ratio",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(colour = "grey40", size = 10),
    axis.text.y      = element_text(size = 8),
    legend.position  = "right",
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "grey92", colour = NA),
    strip.text         = element_text(face = "bold", size = 10)
  )

if (has_category) {
  dot_plot <- dot_plot + facet_wrap(~ Category, scales = "free_y", ncol = 2)
  h <- max(plot_height, length(unique(dot_data$Category)) * 3)
} else {
  h <- plot_height
}

ggsave(file.path(output_dir, "01_DotPlot.png"),
       plot = dot_plot, width = plot_width, height = h, dpi = 300, limitsize = FALSE)
cat("  Saved: 01_DotPlot.png\n\n")

# ==============================================================================
# SECTION 2: BAR CHART
# ==============================================================================
# Bar length = Odds Ratio
# Fill       = -log10(p_FDR)
# Faceted by Category (if present)
# Top N gene sets per category, ordered by Odds Ratio
# ==============================================================================

cat("Plotting Section 2: Bar Chart...\n")

bar_data <- top_per_cat(ora) %>%
  arrange(if (has_category) Category else NULL, Odds_Ratio) %>%
  mutate(Gene_Set_Wrapped = factor(Gene_Set_Wrapped,
                                   levels = unique(Gene_Set_Wrapped)))

bar_plot <- ggplot(bar_data,
                   aes(x    = Odds_Ratio,
                       y    = Gene_Set_Wrapped,
                       fill = neg_log10_FDR)) +
  geom_col(width = 0.75) +
  scale_fill_gradientn(
    colours = c("#4575b4", "#91bfdb", "#fee090", "#fc8d59", "#d73027"),
    name    = expression(-log[10](FDR))
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  labs(
    title    = "ORA Bar Chart — Top Gene Sets by Category",
    subtitle = paste0("Top ", top_n, " gene sets per category  |  Bar length = Odds Ratio"),
    x        = "Odds Ratio",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(colour = "grey40", size = 10),
    axis.text.y        = element_text(size = 8),
    legend.position    = "right",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "grey92", colour = NA),
    strip.text         = element_text(face = "bold", size = 10)
  )

if (has_category) {
  bar_plot <- bar_plot + facet_wrap(~ Category, scales = "free_y", ncol = 2)
  h <- max(plot_height, length(unique(bar_data$Category)) * 3)
} else {
  h <- plot_height
}

ggsave(file.path(output_dir, "02_BarChart.png"),
       plot = bar_plot, width = plot_width, height = h, dpi = 300, limitsize = FALSE)
cat("  Saved: 02_BarChart.png\n\n")

# ==============================================================================
# SECTION 3: VOLCANO-STYLE SCATTER PLOT
# ==============================================================================
# x-axis  = log2(Odds Ratio) — centred on 0 (no enrichment)
# y-axis  = -log10(p_FDR)
# colour  = Significant (TRUE/FALSE)
# size    = N_Overlap
# Labels on top hits (ggrepel, non-overlapping)
# ==============================================================================

cat("Plotting Section 3: Volcano-style Scatter Plot...\n")

vol_data <- ora %>%
  filter(!is.na(Odds_Ratio), Odds_Ratio > 0) %>%
  mutate(
    log2_OR   = log2(Odds_Ratio),
    label_col = case_when(
      Significant & log2_OR > 0 ~ "Enriched (FDR \u2264 0.1)",
      Significant & log2_OR < 0 ~ "Depleted (FDR \u2264 0.1)",
      TRUE                       ~ "Not significant"
    )
  )

# Label only the top N by FDR (to avoid overcrowding)
label_data <- vol_data %>%
  filter(Significant) %>%
  slice_min(order_by = p_FDR, n = top_n, with_ties = FALSE)

point_colours <- c(
  "Enriched (FDR \u2264 0.1)" = "#d73027",
  "Depleted (FDR \u2264 0.1)" = "#4575b4",
  "Not significant"             = "grey70"
)

vol_plot <- ggplot(vol_data,
                   aes(x      = log2_OR,
                       y      = neg_log10_FDR,
                       colour = label_col,
                       size   = N_Overlap)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(fdr_threshold),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_text_repel(
    data         = label_data,
    aes(label    = str_trunc(Gene_Set, 40)),
    size         = 2.8,
    colour       = "black",
    segment.colour = "grey60",
    segment.size = 0.3,
    max.overlaps = 20,
    box.padding  = 0.4
  ) +
  scale_colour_manual(values = point_colours, name = NULL) +
  scale_size_continuous(name  = "Overlap\n(proteins)",
                        range = c(1.5, 8),
                        breaks = pretty_breaks(n = 4)) +
  labs(
    title    = "ORA Volcano Plot — Enrichment vs Significance",
    subtitle = paste0("Dashed lines: FDR = ", fdr_threshold, "  |  log\u2082(OR) = 0 (no enrichment)"),
    x        = expression(log[2](Odds~Ratio)),
    y        = expression(-log[10](FDR))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(colour = "grey40", size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

if (has_category) {
  vol_plot <- vol_plot +
    facet_wrap(~ Category, scales = "free", ncol = 2) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text        = element_text(face = "bold", size = 10))
  h <- max(plot_height, length(unique(vol_data$Category)) * 3.5)
} else {
  h <- plot_height
}

ggsave(file.path(output_dir, "03_VolcanoPlot.png"),
       plot = vol_plot, width = plot_width, height = h, dpi = 300, limitsize = FALSE)
cat("  Saved: 03_VolcanoPlot.png\n\n")

# ==============================================================================
# SECTION 4: HEATMAP
# ==============================================================================
# Rows    = Top significant gene sets (pooled across categories)
# Columns = Categories (if multiple) OR a single value column
# Values  = -log10(p_FDR); grey = not significant / not tested
#
# If only one category is present, columns become the top gene sets and
# the single column shows -log10(p_FDR) as a ranked bar heatmap.
# ==============================================================================

cat("Plotting Section 4: Heatmap...\n")

if (has_category && length(unique(ora$Category)) > 1) {
  
  # Multi-category: pivot so rows = gene sets, columns = categories
  sig_sets <- ora %>%
    filter(Significant) %>%
    group_by(Category) %>%
    slice_min(order_by = p_FDR, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    pull(Gene_Set) %>%
    unique()
  
  if (length(sig_sets) == 0) {
    warning("No significant gene sets found at FDR <= ", fdr_threshold,
            ". Heatmap will use top ", top_n, " by p_FDR regardless of significance.")
    sig_sets <- ora %>%
      group_by(Category) %>%
      slice_min(order_by = p_FDR, n = ceiling(top_n / 2), with_ties = FALSE) %>%
      ungroup() %>%
      pull(Gene_Set) %>%
      unique()
  }
  
  heat_wide <- ora %>%
    filter(Gene_Set %in% sig_sets) %>%
    select(Gene_Set, Category, neg_log10_FDR) %>%
    pivot_wider(names_from  = Category,
                values_from = neg_log10_FDR,
                values_fill = 0)
  
  heat_mat <- as.matrix(heat_wide[, -1])
  rownames(heat_mat) <- str_trunc(heat_wide$Gene_Set, 55)
  
  sig_line <- -log10(fdr_threshold)
  
  # Replace 0 (not tested in that category) with NA for white display
  heat_mat[heat_mat == 0] <- NA
  
  colour_breaks <- seq(0, max(heat_mat, na.rm = TRUE), length.out = 101)
  heat_colours  <- colorRampPalette(c("white", "#fee090", "#fc8d59",
                                      "#d73027", "#800026"))(100)
  
  png(file.path(output_dir, "04_Heatmap.png"),
      width = plot_width * 100, height = max(plot_height, nrow(heat_mat) * 0.3 + 3) * 100,
      res = 100)
  
  pheatmap(
    heat_mat,
    color            = heat_colours,
    breaks           = colour_breaks,
    na_col           = "grey92",
    cluster_rows     = TRUE,
    cluster_cols     = TRUE,
    border_color     = "white",
    fontsize_row     = 8,
    fontsize_col     = 10,
    angle_col        = 45,
    main             = paste0("ORA Heatmap  |  -log\u2081\u2080(FDR)  |  Grey = not tested"),
    legend_breaks    = round(seq(0, max(heat_mat, na.rm = TRUE), length.out = 6), 1),
    legend_labels    = round(seq(0, max(heat_mat, na.rm = TRUE), length.out = 6), 1),
    annotation_legend = TRUE
  )
  dev.off()
  cat("  Saved: 04_Heatmap.png\n\n")
  
} else {
  
  # Single category: simple ranked heatmap (1 column)
  heat_data <- ora %>%
    slice_min(order_by = p_FDR, n = top_n * 2, with_ties = FALSE) %>%
    arrange(p_FDR) %>%
    mutate(Gene_Set_Short = str_trunc(Gene_Set, 55),
           Gene_Set_Short = factor(Gene_Set_Short,
                                   levels = rev(unique(Gene_Set_Short))))
  
  heat_single <- ggplot(heat_data,
                        aes(x = 1, y = Gene_Set_Short, fill = neg_log10_FDR)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_fill_gradientn(
      colours = c("white", "#fee090", "#fc8d59", "#d73027", "#800026"),
      name    = expression(-log[10](FDR))
    ) +
    geom_text(aes(label = round(neg_log10_FDR, 2)), size = 3, colour = "black") +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title    = "ORA Heatmap — Significance Ranked",
         subtitle = paste0("Top ", top_n * 2, " gene sets by FDR"),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(colour = "grey40", size = 10),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      legend.position = "right"
    )
  
  ggsave(file.path(output_dir, "04_Heatmap.png"),
         plot = heat_single, width = 9, height = plot_height, dpi = 300)
  cat("  Saved: 04_Heatmap.png\n\n")
}

#==============================================================================
  # SECTION 5: NETWORK ENRICHMENT MAP
  # ==============================================================================
# Nodes   = significant gene sets (top N per category)
# Edges   = Jaccard similarity of overlapping proteins between pairs of sets
#           (edge drawn only when Jaccard >= jaccard_threshold)
# Node size    = N_Overlap
# Node colour  = -log10(p_FDR)
# Node shape   = Category (if present)
# Edge width   = Jaccard similarity
# ==============================================================================

cat("Plotting Section 5: Network Enrichment Map...\n")

# Minimum Jaccard similarity to draw an edge between two gene sets
jaccard_threshold <- 0.10

# Top N significant gene sets to include as nodes
network_top_n <- 40

net_data <- ora %>%
  filter(Significant) %>%
  arrange(p_FDR) %>%
  slice_head(n = network_top_n)

if (nrow(net_data) < 2) {
  warning("Fewer than 2 significant gene sets — cannot build network. ",
          "Relaxing to top ", network_top_n, " regardless of significance.")
  net_data <- ora %>%
    arrange(p_FDR) %>%
    slice_head(n = network_top_n)
}

# Parse protein membership lists for each node
parse_proteins <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  # Strip [truncated] tag if present
  x <- sub("; \\[truncated\\]", "", x)
  trimws(strsplit(x, "; ", fixed = TRUE)[[1]])
}

net_data$protein_list <- lapply(net_data$Mapped_Proteins, parse_proteins)

# Compute pairwise Jaccard similarity
n_nodes <- nrow(net_data)
edge_list <- data.frame(from = integer(), to = integer(),
                        jaccard = numeric(), stringsAsFactors = FALSE)

for (i in seq_len(n_nodes - 1)) {
  for (j in (i + 1):n_nodes) {
    a <- net_data$protein_list[[i]]
    b <- net_data$protein_list[[j]]
    inter <- length(intersect(a, b))
    if (inter == 0) next
    jacc  <- inter / length(union(a, b))
    if (jacc >= jaccard_threshold) {
      edge_list <- rbind(edge_list,
                         data.frame(from = i, to = j, jaccard = jacc,
                                    stringsAsFactors = FALSE))
    }
  }
}

# Build igraph object using make_empty_graph + add_edges.
# This avoids the V<-/graph_from_data_frame indexing mismatch when integer
# edge indices are passed alongside a separate named-vertices dataframe.

g <- make_empty_graph(n = n_nodes, directed = FALSE)

# Assign vertex attributes directly by position (safe with make_empty_graph)
V(g)$label         <- str_trunc(net_data$Gene_Set, 35)
V(g)$N_Overlap     <- net_data$N_Overlap
V(g)$neg_log10_FDR <- net_data$neg_log10_FDR
if (has_category && "Category" %in% names(net_data)) {
  V(g)$Category <- net_data$Category
}

# Add edges and their Jaccard weights (if any)
has_edges <- nrow(edge_list) > 0
if (has_edges) {
  # add_edges expects a flat interleaved vector: c(from1, to1, from2, to2, ...)
  g <- add_edges(g, as.vector(rbind(edge_list$from, edge_list$to)))
  E(g)$jaccard <- edge_list$jaccard
}

# Compute layout outside of ggraph so it is always a plain numeric matrix,
# preventing the "unknown column x" warning on malformed graph objects.
set.seed(42)
if (has_edges) {
  lay <- create_layout(g, layout = "fr")
} else {
  lay <- create_layout(g, layout = "grid")
}

# Build plot layers conditionally to avoid scale conflicts when there are no edges
edge_layers <- if (has_edges) {
  list(
    geom_edge_link(aes(width = jaccard, alpha = jaccard),
                   colour = "grey60", show.legend = FALSE),
    scale_edge_width_continuous(range = c(0.3, 2.5), guide = "none"),
    scale_edge_alpha_continuous(range = c(0.2, 0.8),  guide = "none")
  )
} else {
  list()
}

net_plot <- ggraph(lay) +
  edge_layers +
  geom_node_point(aes(size   = N_Overlap,
                      colour = neg_log10_FDR),
                  alpha = 0.85) +
  geom_node_text(aes(label = label),
                 size           = 2.5,
                 repel          = TRUE,
                 colour         = "black",
                 segment.colour = "grey70",
                 segment.size   = 0.2,
                 max.overlaps   = 30) +
  scale_colour_gradientn(
    colours = c("#4575b4", "#91bfdb", "#fee090", "#fc8d59", "#d73027"),
    name    = expression(-log[10](FDR))
  ) +
  scale_size_continuous(name   = "Overlap\n(proteins)",
                        range  = c(3, 14),
                        breaks = pretty_breaks(n = 4)) +
  labs(
    title    = "ORA Network Enrichment Map",
    subtitle = paste0("Top ", nrow(net_data), " gene sets  |  ",
                      "Edges: Jaccard similarity \u2265 ", jaccard_threshold,
                      if (!has_edges) "  |  No edges at this threshold — try lowering jaccard_threshold" else "")
  ) +
  theme_graph(base_size = 11, background = "white") +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(colour = "grey40", size = 10),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "05_NetworkMap.png"),
       plot   = net_plot,
       width  = plot_width + 2,
       height = plot_height + 2,
       dpi    = 300)
cat("  Saved: 05_NetworkMap.png\n\n")

# ==============================================================================
# DONE
# ==============================================================================

cat("All plots saved to:", output_dir, "\n")
cat("Files generated:\n")
cat("  01_DotPlot.png\n")
cat("  02_BarChart.png\n")
cat("  03_VolcanoPlot.png\n")
cat("  04_Heatmap.png\n")
cat("  05_NetworkMap.png\n")
