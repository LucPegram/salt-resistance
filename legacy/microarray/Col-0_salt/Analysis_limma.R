# Per-gene test for a linear trend of expression with salt (NaCl) treatment
#
# Input  : long-format file with columns
#          `AGI MAPPING` | treatment (PN, 50, 75, 100) | batch (B1-B3) | values
# Model  : expression ~ dose (+ batch), with PN coded as dose 0
# Method : limma (fast, vectorised; empirical-Bayes moderated t-test on the dose slope)
#          with Benjamini-Hochberg FDR across all genes
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

tissue = 'Root'

for (tissue in c("Shoot", "Root")) {

infile         <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/Col-0_SaltHydroponics/",
                         tissue,
                         "/Col-0_",
                         tissue,
                         "_long.csv")
out_dir        <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/Col-0_SaltHydroponics/",
                         tissue,
                         "/R_output")

fdr_cutoff <- 0.05    # adjusted P-value threshold
use_batch  <- FALSE    # TRUE: batch as a blocking covariate; FALSE: ignore batch
min_obs    <- 12       # minimum non-missing values per gene to be tested (max = 12)

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
message(nrow(expr), " genes tested")

# ---- Design matrix ----------------------------------------------------------
design <- if (use_batch) {model.matrix(~ dose + batch, data = meta)} else {
  model.matrix(~ dose, data = meta)}

# ---- Fit: moderated test of the dose slope ----------------------------------
fit <- lmFit(expr, design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)

res <- topTable(fit, coef = "dose", number = Inf, adjust.method = "BH",
                sort.by = "P", confint = TRUE) %>%
  rownames_to_column("AGI MAPPING") %>%
  as_tibble() %>%
  transmute(
    `AGI MAPPING`,
    slope_per_unit   = logFC,            # change in values per 1 unit of treatment
    slope_per_100    = logFC * 100,      # change in values across 0 -> 100
    CI_low_per_unit  = CI.L,
    CI_high_per_unit = CI.R,
    mean_expression  = AveExpr,
    t_moderated      = t,
    P.Value,
    adj.P.Val,
    direction        = if_else(logFC > 0, "increasing", "decreasing"),
    significant      = adj.P.Val < fdr_cutoff
  )

sig <- filter(res, significant)
message(nrow(sig), " genes with a significant linear trend (BH FDR < ", fdr_cutoff, "): ",
        sum(sig$direction == "increasing"), " increasing, ",
        sum(sig$direction == "decreasing"), " decreasing")

write_csv(res, file.path(out_dir, "linear_trend_results_all.csv"))
write_csv(sig, file.path(out_dir, "linear_trend_results_significant.csv"))

# ---- Diagnostics and plots --------------------------------------------------
png(file.path(out_dir, "pvalue_histogram.png"), width = 1200, height = 900, res = 200)
hist(res$P.Value, breaks = 50, main = "P-value distribution (dose slope)",
     xlab = "P-value", col = "grey80", border = "white")
dev.off()

png(file.path(out_dir, "mean_variance_trend.png"), width = 1200, height = 900, res = 200)
plotSA(fit, main = "Residual SD vs mean expression")
dev.off()

# Volcano plot
p_volcano <- ggplot(res, aes(slope_per_unit, -log10(P.Value), colour = significant)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "firebrick")) +
  labs(x = "Slope (change in values per unit treatment)", y = expression(-log[10](P)),
       colour = paste0("FDR < ", fdr_cutoff)) +
  theme_bw()
ggsave(file.path(out_dir, "volcano_plot.png"), p_volcano, width = 7, height = 5, dpi = 300)

# Top 12 genes: data with fitted lines
top_ids <- head(res$`AGI MAPPING`, 12)
p_top <- long %>%
  filter(`AGI MAPPING` %in% top_ids) %>%
  mutate(`AGI MAPPING` = factor(`AGI MAPPING`, levels = top_ids)) %>%
  ggplot(aes(dose, values)) +
  geom_point(aes(colour = batch), size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.5) +
  facet_wrap(~ `AGI MAPPING`, scales = "free_y") +
  labs(x = "Treatment (PN = 0)", y = "Expression (log-transformed)") +
  theme_bw()
ggsave(file.path(out_dir, "top12_trends.png"), p_top, width = 10, height = 7, dpi = 300)

}