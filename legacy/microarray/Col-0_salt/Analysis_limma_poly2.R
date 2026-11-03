# Per-gene test for a polynomial (degree 2) trend of expression with salt (NaCl) treatment
#
# Input  : long-format file with columns
#          `AGI MAPPING` | treatment (PN, 50, 75, 100) | batch (B1-B3) | values
# Model  : expression ~ linear + quadratic (+ batch), with PN coded as dose 0
#          (linear and quadratic are orthogonal polynomial terms from poly(dose, 2))
# Method : limma (empirical-Bayes moderated tests) with Benjamini-Hochberg FDR across all genes
#          Main test : joint F-test of linear + quadratic terms (is there ANY dose response?)
#          Follow-up : quadratic term alone (is the response curved?)
# Output : full results table, significant-genes table, and diagnostic plots,
#          all saved in the directory set by `out_dir`

rm(list=ls())

# ---- Packages ---------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(ggplot2)

# ---- Settings ---------------------------------------------------------------
#tissue = "Root"
for (tissue in c("Shoot", "Root")) {
  
  infile         <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/Col-0_SaltHydroponics/",
                           tissue,
                           "/Col-0_",
                           tissue,
                           "_long.csv")
  out_dir        <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/Col-0_SaltHydroponics/",
                           tissue,
                           "/R_output/poly2")
  
  fdr_cutoff <- 0.05    # adjusted P-value threshold
  use_batch  <- FALSE   # TRUE: batch as a blocking covariate; FALSE: ignore batch
  min_obs    <- 12      # minimum non-missing values per gene to be tested (max = 12)
  
  # ---- Output directory -------------------------------------------------------
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, winslash = "/")
  
  # ---- Read data; code treatment as a numeric dose (PN = 0) -------------------
  long <- read_csv(
    infile,
    col_types = cols(
      `AGI MAPPING` = col_character(),
      treatment     = col_character(),
      batch         = col_character(),
      values        = col_double()
    )
  ) %>%
    mutate(dose   = as.numeric(sub("^PN$", "0", treatment)),
           sample = paste(treatment, batch, sep = "_"))
  
  # ---- Long -> expression matrix (genes x samples) + sample table -------------
  meta <- long %>%
    distinct(sample, batch, dose) %>%
    arrange(dose, batch) %>%
    mutate(batch = factor(batch))
  
  expr <- long %>%
    select(`AGI MAPPING`, sample, values) %>%
    pivot_wider(names_from = sample, values_from = values) %>%
    column_to_rownames("AGI MAPPING") %>%
    as.matrix()
  expr <- expr[, meta$sample]                         # match column order to `meta`
  
  # Keep genes with enough observations
  expr <- expr[rowSums(!is.na(expr)) >= min_obs, ]
  message(tissue, ": ", nrow(expr), " genes tested")
  
  # ---- Design matrix: degree-2 orthogonal polynomial of dose ------------------
  dose_poly      <- poly(meta$dose, degree = 2)
  meta$linear    <- dose_poly[, 1]
  meta$quadratic <- dose_poly[, 2]
  
  design <- if (use_batch) {
    model.matrix(~ linear + quadratic + batch, data = meta)
  } else {
    model.matrix(~ linear + quadratic, data = meta)
  }
  
  # ---- Fit: moderated tests of the polynomial terms ---------------------------
  fit <- lmFit(expr, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  # Joint F-test of linear + quadratic (any dose response); rows stay in the order of `fit`
  joint <- topTable(fit, coef = c("linear", "quadratic"), number = Inf,
                    sort.by = "none", adjust.method = "BH")
  
  # Measured change in mean expression from dose 0 to dose 100 (raw data, not model-fitted)
  change_0_to_100 <- rowMeans(expr[, meta$dose == 100, drop = FALSE], na.rm = TRUE) -
    rowMeans(expr[, meta$dose == 0,   drop = FALSE], na.rm = TRUE)
  
  res <- tibble(
    `AGI MAPPING`   = rownames(fit),
    mean_expression = joint$AveExpr,
    coef_linear     = fit$coefficients[, "linear"],
    coef_quadratic  = fit$coefficients[, "quadratic"],
    change_0_to_100 = change_0_to_100,           # measured change in mean values, dose 0 to 100
    F_joint         = joint$F,
    P_joint         = joint$P.Value,
    adj.P_joint     = joint$adj.P.Val,
    P_linear        = fit$p.value[, "linear"],
    P_quadratic     = fit$p.value[, "quadratic"],
    adj.P_quadratic = p.adjust(fit$p.value[, "quadratic"], method = "BH")
  ) %>%
    mutate(
      direction             = if_else(change_0_to_100 > 0, "increasing", "decreasing"),
      curvature             = if_else(coef_quadratic > 0, "convex", "concave"),
      significant           = adj.P_joint < fdr_cutoff,
      significant_curvature = adj.P_quadratic < fdr_cutoff
    ) %>%
    arrange(P_joint)
  
  sig <- filter(res, significant)
  message(tissue, ": ", nrow(sig), " genes with a significant dose response (BH FDR < ", fdr_cutoff, "): ",
          sum(sig$direction == "increasing"), " increasing, ",
          sum(sig$direction == "decreasing"), " decreasing; ",
          sum(sig$significant_curvature), " of these with significant curvature")
  
  write_csv(res, file.path(out_dir, "poly2_results_all.csv"))
  write_csv(sig, file.path(out_dir, "poly2_results_significant.csv"))
  
  # ---- Diagnostics and plots --------------------------------------------------
  png(file.path(out_dir, "pvalue_histogram.png"), width = 1200, height = 900, res = 200)
  hist(res$P_joint, breaks = 50, main = "P-value distribution (joint test, linear + quadratic)",
       xlab = "P-value", col = "grey80", border = "white")
  dev.off()
  
  png(file.path(out_dir, "mean_variance_trend.png"), width = 1200, height = 900, res = 200)
  plotSA(fit, main = "Residual SD vs mean expression")
  dev.off()
  
  # Volcano plot
  p_volcano <- ggplot(res, aes(change_0_to_100, -log10(P_joint), colour = significant)) +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "firebrick")) +
    labs(x = "Measured change in mean values, dose 0 to 100", y = expression(-log[10](P[joint])),
         colour = paste0("FDR < ", fdr_cutoff)) +
    theme_bw()
  ggsave(file.path(out_dir, "volcano_plot.png"), p_volcano, width = 7, height = 5, dpi = 300)
  
  # Top 12 genes: data with fitted degree-2 curves
  top_ids <- head(res$`AGI MAPPING`, 12)
  p_top <- long %>%
    filter(`AGI MAPPING` %in% top_ids) %>%
    mutate(`AGI MAPPING` = factor(`AGI MAPPING`, levels = top_ids)) %>%
    ggplot(aes(dose, values)) +
    geom_point(aes(colour = batch), size = 1.8) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, colour = "black", linewidth = 0.5) +
    facet_wrap(~ `AGI MAPPING`, scales = "free_y") +
    labs(x = "Treatment (PN = 0)", y = "Expression (log-transformed)") +
    theme_bw()
  ggsave(file.path(out_dir, "top12_trends.png"), p_top, width = 10, height = 7, dpi = 300)
  
}