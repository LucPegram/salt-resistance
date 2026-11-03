# Volcano Plots — Log Fold Change vs P-value
# Output: PNG of volcano plots and table of proteins
# Written by Claude, modified by Luc Pegram
# 2026/05/22
#===========================
library(ggplot2)
library(ggrepel)
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



#===========================
## Genotype
Genotype  = T
#===========================
## Treatment
#Genotype  = F; Treatment = T
#===========================



for (tissue in c("Shoot", "Root")) {
  
  
  anova <- read.csv(
    paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",tissue,"Tissue_AnovaPvalues.csv"),
    header = T,
    check.names = F,
    stringsAsFactors = F
  )
  
  #anova = distinct(anova[,-1]) #removes duplicated rows due to multiple AGI matching
  
  
  long_df <- read.csv(
    paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",tissue,"Tissue_long.csv"),
    header = T,
    check.names = F,
    stringsAsFactors = F
  )  #for plotting function legends
  
  
  # Adjust to match your actual protein column name
  protein_col <- "Protein.Names"
  log_FC      <- 'log.Genotype'                                        ##<-- Change
  
  
  if (Genotype) {p_col = "p_genotype"} else {p_col = "p_treatment"}
  
  
  output_pdf <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/Volcano/",
                       tissue,
                       "Tissue_VolcanoPlots_Top30Pval_",
                       "Genotype",                                    ##<-- Change
                       ".png")          
  
  output_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/Volcano/",
                        tissue,
                        "Tissue_VolcanoPlot_",
                        "Genotype",                                   ##<-- Change
                        "_Table.csv")    
  
  
  # Significance thresholds
  p_threshold  <- 0.05     # raw p-value cutoff for significance line
  fc_threshold <- log2(1.5)      # log2 fold change cutoff for significance line
  
  # Number of top proteins to label per direction (up/down)
  n_labels <- 30
  
  
  classify <- function(df, p_col, log_col) {
    df$significance <- "Not significant"
    df$significance[df[[p_col]]  < p_threshold &
                      df[[log_col]] >  fc_threshold] <- "Up"
    df$significance[df[[p_col]]  < p_threshold &
                      df[[log_col]] < -fc_threshold] <- "Down"
    df$significance <- factor(df$significance,
                              levels = c("Up", "Down", "Not significant"))
    df
  }
  
  anova  <- classify(anova,  p_col = p_col, log_col = log_FC)
  
  
  # Print counts
  cat("--- Genotype significant proteins ---\n")
  print(table(anova$significance))
  
  
  # SELECT PROTEINS TO LABEL --------------------------------------------------
  # Label the top n_labels Up and Down FC proteins by p-value in each direction
  
  label_proteins <- function(df, sort, Genes) {
    
    df$label <- ""
    
    # Get row indices for each direction
    up_idx   <- which(df$significance == "Up")
    down_idx <- which(df$significance == "Down")
    
    # Sort each set by fold change and keep the top n_labels
    up_sorted   <- up_idx[order(df[[sort]][up_idx],     decreasing = F)]
    down_sorted <- down_idx[order(df[[sort]][down_idx], decreasing = F)]
    
    up_top   <- head(up_sorted,   n_labels)
    down_top <- head(down_sorted, n_labels)
    
    label_idx <- c(up_top, down_top)
    
    lbl <- df$Genes[label_idx]
    df$label[label_idx] <- sub(";.*", "", lbl)
    
    df
  }
  
  anova  <- label_proteins(anova, sort = p_col, Genes = "Genes")
  
  # PLOT FUNCTION --------------------------------------------------
  
  sig_colours <- c(
    "Up"              = "#D6604D",
    "Down"            = "#2166AC",
    "Not significant" = "grey70"
  )
  
  
  make_volcano <- function(df, p_col, title, x_label, log2FC) {
    
    # Dynamic axis limits with a small margin
    x_lim <- max(abs(df[[log2FC]]), na.rm = TRUE) * 1.05
    y_lim <- max(df$neg_log10_p,  na.rm = TRUE) * 1.08
    
    n_up   <- sum(df$significance == "Up")
    n_down <- sum(df$significance == "Down")
    
    ggplot(df, aes(x = .data[[log2FC]], y = neg_log10_p,
                   colour = significance, label = label)) +
      
      # Threshold lines
      geom_hline(yintercept = -log10(p_threshold),
                 linetype = "dashed", colour = "grey40", linewidth = 0.5) +
      geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                 linetype = "dashed", colour = "grey40", linewidth = 0.5) +
      
      # Points — non-significant first so significant sit on top
      geom_point(data = df[df$significance == "Not significant", ],
                 size = 1.5, alpha = 0.5) +  
      geom_point(data = df[df$significance != "Not significant", ],
                 size = 2.0, alpha = 0.8) +  #alpha is transparency
      
      # Labels
      geom_text_repel(
        data          = df[df$label != "", ],
        size          = 3,
        box.padding   = 0.5,
        max.overlaps  = 20,
        segment.size  = 0.3,
        segment.colour = "grey50",
        show.legend   = FALSE
      ) +
      # Significance counts in corners
      annotate("text", x = -x_lim * 0.97, y = y_lim * 0.97,
               label = sprintf("Down: %d", n_down),
               colour = sig_colours["Down"], hjust = 0, size = 4,
               fontface = "bold") +
      annotate("text", x =  x_lim * 0.97, y = y_lim * 0.97,
               label = sprintf("Up: %d", n_up),
               colour = sig_colours["Up"],   hjust = 1, size = 4,
               fontface = "bold") +
      scale_colour_manual(values = sig_colours[1:2], name = NULL) +
      scale_x_continuous(limits = c(-x_lim, x_lim), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0,      y_lim), expand = c(0, 0)) +
      labs(
        title    = title,
        subtitle = sprintf("p < %.2f  |  |FC| > log2(1.5)  |  Labels: top %d per direction",
                           p_threshold, n_labels),
        x        = x_label,
        y        = expression(-log[10](p))
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title       = element_text(face = "bold", size = 22),
        plot.subtitle    = element_text(colour = "grey40", size = 10),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom"
      )
  }
  
  
  # GENERATE AND SAVE PLOTS -------------------------------------------
  
  genotype_levels  <- unique(long_df$genotype)
  treatment_levels <- unique(long_df$treatment)
  
  if (Genotype) {
    genotype_ref  <- genotype_levels[1]
    genotype_comp <- genotype_levels[2]
  } else {
    genotype_ref  <- treatment_levels[1]
    genotype_comp <- treatment_levels[2]
  }
  anova$neg_log10_p <- -log10(anova[[p_col]])  #Transforms p-values for display in volcano plot.
  
  # ------------------------------------------------------------------
  
  p_geno <- make_volcano(
    df      = anova,
    p_col   = p_col,
    title   = bquote(.(tissue) * ": " * italic("NIT2-OE") * " vs " * .(genotype_ref)),
    x_label = bquote(log[2] ~ "fold change (" * .(genotype_comp) * " / " * .(genotype_ref) * ")"),
    log2FC  = log_FC
  )
  
  #print(p_geno)
  
  
  png(output_pdf, width = 12, height = 10, units = "in", res = 300)
  print(p_geno)
  dev.off()
  
  cat("\nVolcano plots saved to ", output_pdf, "\n", sep = "")
  
  # EXPORT SIGNIFICANT PROTEIN TABLES ----------------------------
  
  sig_geno <- anova[anova$significance != "Not significant",
                    !colnames(anova) %in% c("neg_log10_p", "label")]
  
  sig_geno <- sig_geno[order(sig_geno[[p_col]]), ]
  
  write.csv(sig_geno, output_file, row.names = FALSE)
  
}
