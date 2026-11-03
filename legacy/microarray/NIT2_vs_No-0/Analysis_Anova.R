# Two-way ANOVA per gene: genotype (BG, OE) x treatment (PN, 100)
# Input : output of convert_to_long_genotype.R
# Terms : genotype, treatment, genotype:treatment interaction
# FDR   : Benjamini-Hochberg across all genes, separately for each term

# ---- Packages ---------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(tibble)

# ---- Settings ---------------------------------------------------------------
infile <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/NIT2vsNo-0_Treatments/",
                 "NIT2vsNo-0_long.csv")
out_dir <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/NIT2vsNo-0_Treatments/",
                  "R_analysis/Two-wayAnova")
fdr_cutoff <- 0.05

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Read data and build expression matrix (genes x samples) ----------------
long <- read_csv(infile, show_col_types = FALSE) %>%
  mutate(sample = paste(genotype, treatment, replicate, sep = "_"))

meta <- long %>%
  distinct(sample, genotype, treatment) %>%
  mutate(genotype  = factor(genotype,  levels = c("BG", "OE")),
         treatment = factor(treatment, levels = c("PN", "100")))

expr <- long %>%
  select(`AGI MAPPING`, sample, values) %>%
  pivot_wider(names_from = sample, values_from = values) %>%
  column_to_rownames("AGI MAPPING") %>%
  as.matrix()
expr <- expr[, meta$sample]

# ---- Two-way ANOVA for one gene ---------------------------------------------
anova_for_one_gene <- function(y) {
  tab <- anova(lm(y ~ genotype * treatment, data = meta))
  c(F_genotype    = tab["genotype",           "F value"],
    P_genotype    = tab["genotype",           "Pr(>F)"],
    F_treatment   = tab["treatment",          "F value"],
    P_treatment   = tab["treatment",          "Pr(>F)"],
    F_interaction = tab["genotype:treatment", "F value"],
    P_interaction = tab["genotype:treatment", "Pr(>F)"])
}

# ---- Run on every gene ------------------------------------------------------
anova_results <- t(apply(expr, 1, anova_for_one_gene))

res <- as_tibble(anova_results, rownames = "AGI MAPPING") %>%
  mutate(
    adj.P_genotype    = p.adjust(P_genotype,    method = "BH"),
    adj.P_treatment   = p.adjust(P_treatment,   method = "BH"),
    adj.P_interaction = p.adjust(P_interaction, method = "BH"),
    mean_BG = rowMeans(expr[, meta$genotype == "BG"], na.rm = TRUE),
    mean_OE = rowMeans(expr[, meta$genotype == "OE"], na.rm = TRUE),
    log2FC_OE_vs_BG = mean_OE - mean_BG
  ) %>%
  arrange(P_genotype)

# ---- Save and summarise -----------------------------------------------------
write_csv(res, file.path(out_dir, "two_way_anova_results.csv"))

message(sum(res$adj.P_genotype    < fdr_cutoff, na.rm = TRUE), " genes with a significant genotype effect")
message(sum(res$adj.P_treatment   < fdr_cutoff, na.rm = TRUE), " genes with a significant treatment effect")
message(sum(res$adj.P_interaction < fdr_cutoff, na.rm = TRUE), " genes with a significant genotype x treatment interaction")