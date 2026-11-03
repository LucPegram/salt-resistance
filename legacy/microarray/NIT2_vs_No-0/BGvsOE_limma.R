# Per-gene comparison of expression between genotypes BG and OE
#
# Input  : output of convert_to_long_genotype.R (log2-transformed, one row per
#          AGI MAPPING x treatment x genotype x replicate), with columns
#          `AGI MAPPING` | treatment (PN, 100) | genotype (BG, OE) | replicate | values
# Design : the four treatment x replicate samples within each genotype (PN 1, PN 2,
#          100 NaCl 1, 100 NaCl 2) are treated as 4 biological replicates per genotype
# Model  : values ~ genotype  (reference = BG, so positive log2FC = higher in OE)
# Method : limma linear model, Benjamini-Hochberg FDR across all genes
#          (empirical-Bayes moderated t-test, or ordinary t-test: see `moderated`)
# Output : results tables and diagnostic plots, all saved in `out_dir`

# ---- Packages ---------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(ggplot2)

rm(list=ls())

# ---- Settings ---------------------------------------------------------------
infile               <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/NIT2vsNo-0_Treatments/",
                               "NIT2vsNo-0_long.csv")
out_dir              <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/NIT2vsNo-0_Treatments/",
                               "R_analysis/limma_eBayes_t-test")

fdr_cutoff           <- 0.05    # adjusted P-value threshold
min_abs_log2fc       <- 0       # additional fold-change threshold (0 = P-value only; 1 = 2-fold)
moderated            <- FALSE   # TRUE: empirical-Bayes moderated t-test (limma eBayes)
# FALSE: ordinary per-gene t-test, no variance shrinkage
adjust_for_treatment <- TRUE    # FALSE: treatments pooled as replicates (as specified)
# TRUE: model values ~ genotype + treatment
min_obs_per_genotype <- 2       # minimum non-missing values per genotype to test a gene (max = 4)

# ---- Output directory -------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/")

# ---- Read data --------------------------------------------------------------
long <- read_csv(
  infile,
  col_types = cols(
    `AGI MAPPING` = col_character(),
    treatment     = col_character(),
    genotype      = col_character(),
    replicate     = col_character(),
    values        = col_double()
  )
) %>%
  mutate(sample = paste(genotype, treatment, replicate, sep = "_"))

# ---- Long -> expression matrix (genes x samples) + sample table -------------
meta <- long %>%
  distinct(sample, genotype, treatment) %>%
  mutate(genotype  = factor(genotype,  levels = c("BG", "OE")),
         treatment = factor(treatment, levels = c("PN", "100"))) %>%
  arrange(genotype, treatment, sample)

expr <- long %>%
  select(`AGI MAPPING`, sample, values) %>%
  pivot_wider(names_from = sample, values_from = values) %>%
  column_to_rownames("AGI MAPPING") %>%
  as.matrix()
expr <- expr[, meta$sample]                         # match column order to `meta`

is_BG <- meta$genotype == "BG"
is_OE <- meta$genotype == "OE"

# Keep genes with enough observations in both genotypes
n_BG <- rowSums(!is.na(expr[, is_BG]))
n_OE <- rowSums(!is.na(expr[, is_OE]))
expr <- expr[n_BG >= min_obs_per_genotype & n_OE >= min_obs_per_genotype, ]
message(nrow(expr), " genes tested")

# ---- Design matrix ----------------------------------------------------------
design <- if (adjust_for_treatment) {
  model.matrix(~ genotype + treatment, data = meta)
} else {
  model.matrix(~ genotype, data = meta)
}

# ---- Fit linear model -------------------------------------------------------
fit <- lmFit(expr, design)
if (moderated) {
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
}

# ---- Test the genotype coefficient (OE - BG) --------------------------------
log2FC <- fit$coefficients[, "genotypeOE"]

if (moderated) {
  std_error <- sqrt(fit$s2.post) * fit$stdev.unscaled[, "genotypeOE"]
  df        <- fit$df.total
} else {
  std_error <- fit$sigma * fit$stdev.unscaled[, "genotypeOE"]
  df        <- fit$df.residual
}

t_stat  <- log2FC / std_error
p_value <- 2 * pt(-abs(t_stat), df = df)
t_crit  <- qt(0.975, df = df)

res <- tibble(
  `AGI MAPPING` = rownames(expr),
  mean_BG       = rowMeans(expr[, is_BG], na.rm = TRUE),
  mean_OE       = rowMeans(expr[, is_OE], na.rm = TRUE),
  log2FC        = log2FC,                         # OE minus BG (values are log2)
  CI_low        = log2FC - t_crit * std_error,
  CI_high       = log2FC + t_crit * std_error,
  t_stat        = t_stat,
  P.Value       = p_value,
  adj.P.Val     = p.adjust(p_value, method = "BH")
) %>%
  mutate(
    direction   = if_else(log2FC > 0, "higher in OE", "higher in BG"),
    significant = adj.P.Val < fdr_cutoff & abs(log2FC) >= min_abs_log2fc
  ) %>%
  arrange(P.Value)

sig <- filter(res, significant)
message(nrow(sig), " genes differ between genotypes (BH FDR < ", fdr_cutoff, "): ",
        sum(sig$direction == "higher in OE"), " higher in OE, ",
        sum(sig$direction == "higher in BG"), " higher in BG")

write_csv(res, file.path(out_dir, "genotype_comparison_all.csv"))
write_csv(sig, file.path(out_dir, "genotype_comparison_significant.csv"))

# ---- Diagnostics and plots --------------------------------------------------
png(file.path(out_dir, "pvalue_histogram.png"), width = 1200, height = 900, res = 200)
hist(res$P.Value, breaks = 50, main = "P-value distribution (OE vs BG)",
     xlab = "P-value", col = "grey80", border = "white")
dev.off()

# Volcano plot
p_volcano <- ggplot(res, aes(log2FC, -log10(P.Value), colour = significant)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "firebrick")) +
  labs(x = "log2 fold change (OE - BG)", y = expression(-log[10](P)),
       colour = paste0("FDR < ", fdr_cutoff)) +
  theme_bw()
ggsave(file.path(out_dir, "volcano_plot.png"), p_volcano, width = 7, height = 5, dpi = 300)

# Top 12 genes: individual samples by genotype, coloured by treatment, with genotype means
top_ids <- head(res$`AGI MAPPING`, 12)
p_top <- long %>%
  filter(`AGI MAPPING` %in% top_ids) %>%
  mutate(`AGI MAPPING` = factor(`AGI MAPPING`, levels = top_ids),
         genotype      = factor(genotype,  levels = c("BG", "OE")),
         treatment     = factor(treatment, levels = c("PN", "100"))) %>%
  ggplot(aes(genotype, values)) +
  geom_jitter(aes(colour = treatment), width = 0.1, height = 0, size = 2) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymin = after_stat(y), ymax = after_stat(y)), width = 0.4) +
  facet_wrap(~ `AGI MAPPING`, scales = "free_y") +
  labs(x = "Genotype", y = "Expression (log-transformed)") +
  theme_bw()
ggsave(file.path(out_dir, "top12_genes.png"), p_top, width = 10, height = 7, dpi = 300)