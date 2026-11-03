## ============================================================
## ORA Volcano Plots — One plot per annotation Category
## Input : Shoot_ORA_AllCategories.csv
## Output: One .png per Category saved to output_dir_volcano
## ============================================================

library(ggplot2)
library(dplyr)
library(ggrepel)
library(scales)
library(stringr)

# ==============================================================================
# USER INPUT
# ==============================================================================

rm(list=ls())

tissue = "Root"
type = "Combined"


for (tissue in c("Shoot", "Root")) {
  
#  for (type in c("Genotype", "Treatment", "Interaction")) {
  
  
  input_file <- paste0(
    "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/ORA/",
    tissue,
    "/",
    type,
    "/WholeGenome/",
    tissue,
    "_ORA_AllCategories.csv")
  
  output_dir_volcano <- paste0(
    "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/ORA/",
    tissue,
    "/",
    type,
    "/WholeGenome/VolcanoPlots/")
  
  output_dir_bar <- paste0(
    "D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/ORA/",
    tissue,
    "/",
    type,
    "/WholeGenome/BarCharts/")
  
  # FDR threshold — used for significance colouring and the horizontal dashed line
  fdr_threshold <- 0.1
  
  # Maximum number of gene set labels to draw per plot (top hits by FDR)
  top_n_labels <- 20
  
  # Plot dimensions (inches)
  plot_width  <- 12
  plot_height <- 8
  
  # ==============================================================================
  # LOAD AND PREPARE DATA
  # ==============================================================================
  
  dir.create(output_dir_volcano, showWarnings = FALSE, recursive = TRUE)
  
  ora <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE,
                  check.names = FALSE)
  
  if (!"Category" %in% names(ora)) {
    stop("No 'Category' column found. Please supply the AllCategories combined file.")
  }
  
  # Coerce numeric columns
  ora$p_value    <- as.numeric(ora$p_value)
  ora$p_FDR      <- as.numeric(ora$p_FDR)
  ora$Odds_Ratio <- as.numeric(ora$Odds_Ratio)
  ora$N_Overlap  <- as.numeric(ora$N_Overlap)
  
  # Remove rows with NA p-values or non-positive Odds Ratios
  ora <- ora[!is.na(ora$p_value) & !is.na(ora$Odds_Ratio) & ora$Odds_Ratio > 0, ]
  
  # Derived columns
  ora$log2_OR       <- log2(ora$Odds_Ratio)
  ora$neg_log10_FDR <- -log10(ora$p_FDR)
  ora$Significant   <- !is.na(ora$p_FDR) & ora$p_FDR <= fdr_threshold
  
  # --- Cap Inf values -----------------------------------------------------------
  # Fisher's exact test returns Odds_Ratio = Inf for perfect enrichment
  # (log2(Inf) = Inf), and p_FDR = 0 gives -log10(0) = Inf.
  # ggplot2's scale internals cannot handle Inf limits, so we cap at the
  # largest finite value + 10% buffer and flag capped points for transparency.
  
  finite_log2_OR <- ora$log2_OR[is.finite(ora$log2_OR)]
  finite_nlFDR   <- ora$neg_log10_FDR[is.finite(ora$neg_log10_FDR)]
  
  log2_OR_cap <- if (length(finite_log2_OR) > 0) max(abs(finite_log2_OR)) * 1.15 else 10
  nlFDR_cap   <- if (length(finite_nlFDR)   > 0) max(finite_nlFDR)        * 1.15 else 10
  
  ora$log2_OR_capped       <- FALSE
  ora$neg_log10_FDR_capped <- FALSE
  
  ora$log2_OR_capped[!is.finite(ora$log2_OR)] <- TRUE
  ora$log2_OR[ora$log2_OR >  log2_OR_cap]     <-  log2_OR_cap
  ora$log2_OR[ora$log2_OR < -log2_OR_cap]     <- -log2_OR_cap
  
  ora$neg_log10_FDR_capped[!is.finite(ora$neg_log10_FDR)] <- TRUE
  ora$neg_log10_FDR[ora$neg_log10_FDR > nlFDR_cap]        <- nlFDR_cap
  
  n_capped <- sum(ora$log2_OR_capped | ora$neg_log10_FDR_capped)
  if (n_capped > 0) {
    cat(sprintf("Note: %d row(s) had Inf values in log2(OR) or -log10(FDR) and were capped.\n\n",
                n_capped))
  }
  
  ora$label_col <- case_when(
    ora$Significant ~ paste0("Enriched (FDR \u2264 ", fdr_threshold, ")"),
    TRUE            ~ "Not significant"
  )
  
  #point_colours <- c(
  #  paste0("Enriched (FDR \u2264 ", fdr_threshold, ")") = "#d73027",
  #  "Not significant"                                    = "grey70"
  #)
  
  point_colours <- c(
    "Enriched (FDR \u2264 0.1)" = "#d73027",
    "Not significant"             = "grey70"
  )
  
  categories <- unique(ora$Category)
  cat("Categories found:", paste(categories, collapse = ", "), "\n")
  cat("Generating", length(categories), "volcano plot(s)...\n\n")
  
  # ==============================================================================
  # GENERATE ONE VOLCANO PLOT PER CATEGORY
  # ==============================================================================
  
  for (cat_name in categories) {
    
    cat_data <- ora[ora$Category == cat_name, ]
    
    # Gene sets to label: top hits by FDR among significant points
    label_data <- cat_data %>%
      filter(Significant) %>%
      slice_min(order_by = p_FDR, n = top_n_labels, with_ties = FALSE)
    
    # Axis limits — use only finite values with a small margin.
    # y_max must be at least as large as sig_line, otherwise the significance
    # threshold line is clipped outside the scale range and triggers a warning.
    x_abs_max <- max(abs(cat_data$log2_OR[is.finite(cat_data$log2_OR)]),
                     na.rm = TRUE) * 1.1
    sig_line  <- -log10(fdr_threshold)
    y_max     <- max(
      max(cat_data$neg_log10_FDR[is.finite(cat_data$neg_log10_FDR)], na.rm = TRUE) * 1.1,
      sig_line * 1.2
    )
    
    # Fallback if all values are somehow non-finite
    if (!is.finite(x_abs_max) || x_abs_max == 0) x_abs_max <- 5
    if (!is.finite(y_max)     || y_max      == 0) y_max     <- 5
    
    n_enriched <- sum(cat_data$Significant, na.rm = TRUE)
    n_ns       <- sum(!cat_data$Significant, na.rm = TRUE)
    
    p <- ggplot(cat_data,
                aes(x      = log2_OR,
                    y      = neg_log10_FDR,
                    colour = label_col,
                    size   = N_Overlap)) +
      
      # Reference lines
      geom_hline(yintercept = sig_line,
                 linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      geom_vline(xintercept = 0,
                 linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      
      # Points
      geom_point(alpha = 0.75) +
      
      # Non-overlapping labels for top significant gene sets
      geom_text_repel(
        data           = label_data,
        aes(label      = str_trunc(Gene_Set, 45)),
        size           = 2.8,
        colour         = "black",
        segment.colour = "grey55",
        segment.size   = 0.3,
        segment.alpha  = 0.7,
        box.padding    = 0.45,
        point.padding  = 0.3,
        max.overlaps   = 25,
        min.segment.length = 0.2
      ) +
      
      # Scales
      scale_colour_manual(values = point_colours, name = NULL) +
      scale_size_continuous(
        name   = "Overlap\n(proteins)",
        range  = c(1.5, 8),
        breaks = pretty_breaks(n = 4)
      ) +
      scale_x_continuous(limits = c(-x_abs_max, x_abs_max),
                         breaks = pretty_breaks(n = 8)) +
      scale_y_continuous(limits = c(0, y_max),
                         breaks = pretty_breaks(n = 6)) +
      
      # Annotation: enriched count in top-right corner
      annotate("text",
               x = x_abs_max * 0.95, y = y_max * 0.97,
               label = paste0("Enriched: ", n_enriched),
               colour = "#d73027", hjust = 1, size = 3.2, fontface = "bold") +
      
      # Labels
      labs(
        title    = paste0("ORA Volcano Plot — ", cat_name),
        subtitle = paste0(
          "FDR threshold: ", fdr_threshold,
          "  |  log\u2082(OR) = 0: no enrichment",
          "  |  Labelled: top ", min(top_n_labels, nrow(label_data)),
          " significant gene sets"
        ),
        x = expression(log[2](Odds~Ratio)),
        y = expression(-log[10](FDR))
      ) +
      
      # Theme
      theme_bw(base_size = 11) +
      theme(
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(colour = "grey40", size = 9),
        legend.position  = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92")
      )
    
    # Save — sanitise category name for use in filename
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", cat_name)
    out_file  <- file.path(output_dir_volcano, paste0("VolcanoPlot_", safe_name, ".png"))
    
    tryCatch({
      ggsave(out_file, plot = p, width = plot_width, height = plot_height, dpi = 300)
      cat(sprintf("  [%s]  Saved: %s\n",
                  cat_name,
                  basename(out_file)))
      cat(sprintf("         Gene sets: %d total  |  %d enriched  |  %d NS\n\n",
                  nrow(cat_data), n_enriched, n_ns))
    }, error = function(e) {
      # Ensure any open graphics device is closed before continuing
      try(dev.off(), silent = TRUE)
      message(sprintf("  [%s]  ERROR — plot not saved: %s\n", cat_name, conditionMessage(e)))
    })
  }
  
  cat("Done. All volcano plots saved to:\n", output_dir_volcano, "\n")
  
  # ==============================================================================
  # SECTION 2: BAR CHARTS — One per Category
  # ==============================================================================
  # Bar length = Odds Ratio
  # Fill       = -log10(p_FDR)
  # Top N gene sets per category, ordered by Odds Ratio
  # Colour anchors (fixed across all plots):
  #   #4575b4 (blue)       = 0
  #   #91bfdb (light blue) = -log10(0.05) = 1.301
  #   #fee090 (yellow)     = -log10(0.01) = 2
  #   #fc8d59 (orange)     = 3.5
  #   #d73027 (red)        = 5  [values > 5 squished to red]
  # ==============================================================================
  dir.create(output_dir_bar, showWarnings = FALSE, recursive = TRUE)
  
  # Number of top gene sets to show per bar chart
  top_n_bar <- 20
  
  cat("\nGenerating", length(categories), "bar chart(s)...\n\n")
  
  for (cat_name in categories) {
    
    cat_data <- ora[ora$Category == cat_name, ]
    
    # Select top N by FDR, then order bars by Odds Ratio for clean presentation
    bar_data <- cat_data %>%
      slice_min(order_by = p_FDR, n = top_n_bar, with_ties = FALSE) %>%
      arrange(Odds_Ratio) %>%
      mutate(Gene_Set_Wrapped = str_wrap(Gene_Set, width = 50),
             Gene_Set_Wrapped = factor(Gene_Set_Wrapped,
                                       levels = unique(Gene_Set_Wrapped)))
    
    b <- ggplot(bar_data,
                aes(x    = Odds_Ratio,
                    y    = Gene_Set_Wrapped,
                    fill = neg_log10_FDR)) +
      
      geom_col(width = 0.75) +
      
      geom_vline(xintercept = 1,
                 linetype = "dashed", colour = "grey30", linewidth = 0.4) +
      
      scale_fill_gradientn(
        colours = c("#4575b4", "#91bfdb", "#fee090", "#fc8d59", "#d73027"),
        values  = rescale(c(0, -log10(0.1), -log10(0.05), 3.5, 5)),
        limits  = c(0, 5),
        oob     = squish,
        name    = expression(-log[10](FDR))
      ) +
      
      scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                         breaks = pretty_breaks(n = 6)) +
      
      labs(
        title    = paste0("ORA Bar Chart — ", cat_name),
        subtitle = paste0("Top ", nrow(bar_data), " gene sets by FDR  |  ",
                          "Bar length = Odds Ratio  |  ",
                          "Fill = -log\u2081\u2080(FDR)"),
        x        = "Odds Ratio",
        y        = NULL
      ) +
      
      theme_bw(base_size = 11) +
      theme(
        plot.title         = element_text(face = "bold", size = 13),
        plot.subtitle      = element_text(colour = "grey40", size = 9),
        axis.text.y        = element_text(size = 8),
        legend.position    = "right",
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank()
      )
    
    # Scale plot height to the number of bars so labels are never cramped
    bar_height <- max(plot_height, nrow(bar_data) * 0.4 + 2)
    
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", cat_name)
    out_file  <- file.path(output_dir_bar, paste0("BarChart_", safe_name, ".png"))
    
    tryCatch({
      ggsave(out_file, plot = b, width = plot_width, height = bar_height, dpi = 300)
      cat(sprintf("  [%s]  Saved: %s\n", cat_name, basename(out_file)))
      cat(sprintf("         Gene sets shown: %d\n\n", nrow(bar_data)))
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      message(sprintf("  [%s]  ERROR — plot not saved: %s\n", cat_name, conditionMessage(e)))
    })
  }
  
  cat("Done. All bar charts saved to:\n", output_dir_bar, "\n")
  
#  }
}
