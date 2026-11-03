# Creates and exports lists of Gene Set Annotations
# Written by Claude, modified by Luc Pegram
# 2026/05/26
#===========================



rm(list=ls())


# Read in Gene Set Annotation file
gsa_df <- read.delim(
  "D:/Masters_PhD/Data/Proteomics/2026/Ath_GeneSetAnnotations_All.txt",
  header = T, 
  sep = '\t',
  check.names = F,
  stringsAsFactors = F
)


# Split by Category, then for each Category create a named list
# where each element is a Gene set name containing its associated Genes

for (cat in unique(gsa_df$Category)) {
  
  # Subset to this Category
  cat_df <- gsa_df[gsa_df$Category == cat, ]
  
  # For each Gene, get all associated Gene set names
  gsa_list <- tapply(cat_df$Gene, cat_df$`Gene set name`, FUN = c)
  
  # Assign to an object named after the Category value
  assign(cat, gsa_list)
}



#===========================
## Shoot Tissue
Shoot = T; Root = F 
#===========================
## Root Tissue
Root = T;  Shoot = F
#===========================

if (Shoot) {tissue <- "Shoot"} else {tissue <- "Root"}



for (tissue in c("Shoot", "Root")) {

  
# --- Input directory ---
input_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/ForAnalysis/R_Exports/",
                      tissue,
                      "Tissue_AnovaPvalues.csv")

# --- Output directory ---
output_dir <- paste0(
  "D:/Masters_PhD/Data/Proteomics/2026/ForAnalysis/R_Exports/ORA/",
  tissue, "/Genotype/"
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


df <- read.csv(
  input_file,
  header = T,
  check.names = F,
  stringsAsFactors = F
)

df <- dplyr::distinct(df, AGI, .keep_all = TRUE)  # removes duplicated AGI.

## ============================================================
## Over Enrichment Analysis (ORA) — Proteomic Data
## Fisher's Exact Test (one-sided: greater) + BH-FDR correction
## One .csv exported per annotation Category + one combined file
## ============================================================

# --- Define significant proteins (query) and background ---
sig_proteins        <- df$AGI[df$p_FDR_genotype <= 0.05 & 
                                (abs(df$log.Genotype) > 0)]           ## <--Set
background_proteins <- df$AGI

cat("Background proteins :", length(background_proteins), "\n")
cat("Significant proteins:", length(sig_proteins), "\n\n")

if (length(sig_proteins) == 0) {
  stop("No proteins passed the ORA filter. 
       Check column names and thresholds.")
}

# --- Initialise combined results collector ---
combined_results <- list()

# --- Run ORA for each Category ---
for (cat in unique(gsa_df$Category)) {
  
  message("Running ORA for category: ", cat)
  
  # Retrieve the annotation mapping list created in the previous loop
  annot_list <- get(cat)
  
  n_background <- length(background_proteins)
  n_sig        <- length(sig_proteins)
  
  # Pre-allocate results list for efficiency
  results_list <- vector("list", length(annot_list))
  gs_names     <- names(annot_list)
  
  for (i in seq_along(annot_list)) {
    
    gs_genes <- annot_list[[i]]
    
    # Proteins from the full dataset that are members of this gene set
    bg_in_gs   <- background_proteins[background_proteins %in% gs_genes]
    n_gs_in_bg <- length(bg_in_gs)
    
    # Skip gene sets with no representation in the background
    if (n_gs_in_bg == 0) {
      results_list[[i]] <- data.frame(
        Gene_Set            = gs_names[i],
        N_Background        = n_background,
        N_GeneSet           = 0L,
        N_Sig               = n_sig,
        N_Overlap           = 0L,
        Expected            = NA_real_,
        Odds_Ratio          = NA_real_,
        p_value             = NA_real_,
        Mapped_Proteins     = NA_character_,
        Sig_Proteins        = NA_character_,
        stringsAsFactors    = FALSE
      )
      next
    }
    
    # --- 2x2 Contingency Table ---
    #                  | In Gene Set | Not in Gene Set |
    # Significant      |      a      |        b        |
    # Not Significant  |      c      |        d        |
    
    sig_in_gs  <- sig_proteins[sig_proteins %in% gs_genes]                 # sig & in GS
    a <- length(sig_in_gs)
    b <- n_sig - a                                                           # sig & not in GS
    c <- sum(background_proteins[!background_proteins %in% sig_proteins]    # not-sig & in GS
             %in% gs_genes)
    d <- n_background - n_sig - c                                           # not-sig & not in GS
    
    contingency_mat <- matrix(c(a, b, c, d), nrow = 2,
                              dimnames = list(c("Sig", "Not_Sig"),
                                              c("In_GS", "Not_In_GS")))
    
    fisher_res <- fisher.test(contingency_mat, alternative = "greater")
    
    expected <- (n_sig * n_gs_in_bg) / n_background
    
    results_list[[i]] <- data.frame(
      Gene_Set            = gs_names[i],
      N_Background        = n_background,
      N_GeneSet           = n_gs_in_bg,
      N_Sig               = n_sig,
      N_Overlap           = a,
      Expected            = round(expected, 4),
      Odds_Ratio          = round(as.numeric(fisher_res$estimate), 4),
      p_value             = fisher_res$p.value,
      Mapped_Proteins     = paste(bg_in_gs,  collapse = "; "),
      Sig_Proteins        = if (a > 0) paste(sig_in_gs, collapse = "; ") else NA_character_,
      stringsAsFactors    = FALSE
    )
  }
  
  # Combine all gene set results for this category
  results_df <- do.call(rbind, results_list)
  
  # BH-FDR correction across all gene sets in this category
  # NA p-values (gene sets absent from background) are excluded from correction
  results_df$p_FDR <- NA_real_
  valid_idx <- !is.na(results_df$p_value)
  results_df$p_FDR[valid_idx] <- p.adjust(results_df$p_value[valid_idx], method = "BH")
  
  # Sort by raw p-value
  results_df <- results_df[order(results_df$p_value, na.last = TRUE), ]
  
  # Round p-values for readability
  results_df$p_value <- signif(results_df$p_value, 4)
  results_df$p_FDR   <- signif(results_df$p_FDR,   4)
  
  # --- Export per-category file ---
  out_file <- paste0(output_dir, tissue, "_ORA_", cat, ".csv")
  write.csv(results_df, out_file, row.names = FALSE, quote = TRUE)
  
  message("  Exported: ", out_file)
  message("  Gene sets tested : ", nrow(results_df))
  message("  Significant (FDR <= 0.05): ",
          sum(results_df$p_FDR <= 0.05, na.rm = TRUE), "\n")
  
  # --- Store for combined export (add Category column) ---
  results_df$Category     <- cat
  combined_results[[cat]] <- results_df
}

# --- Export combined file ---
combined_df <- do.call(rbind, combined_results)

# Reorder columns so Category appears first
combined_df <- combined_df[, c("Category", setdiff(names(combined_df), "Category"))]

# Sort by p_value across all categories
combined_df <- combined_df[order(combined_df$p_value, na.last = TRUE), ]

combined_file <- paste0(output_dir, tissue, "_ORA_AllCategories.csv")
write.csv(combined_df, combined_file, row.names = FALSE, quote = TRUE)

message("Combined file exported : ", combined_file)
message("Total gene sets across all categories: ", nrow(combined_df))
message("ORA complete. Results saved to: ", output_dir)

}
