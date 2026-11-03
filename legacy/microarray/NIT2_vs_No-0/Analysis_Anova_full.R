# Two-way ANOVA across the entire dataset (all genes pooled) with a gene term
# Model: values ~ gene + genotype * treatment
library(readr)
library(dplyr)
library(car)

infile  <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/NIT2vsNo-0_Treatments/",
                  "NIT2vsNo-0_long.csv")

out_dir <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/NIT2vsNo-0_Treatments/",
                  "R_analysis/Anova_of_full_microarray")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

long_df <- read_csv(infile, show_col_types = FALSE)
long_df$genotype  <- factor(long_df$genotype,  levels = c("BG", "OE"))
long_df$treatment <- factor(long_df$treatment, levels = c("PN", "100"))

# Gene term: subtract each gene's own mean (same as adding gene as a fixed factor,
# without building a 60,000-column model matrix)
long_df <- long_df %>%
  group_by(`AGI MAPPING`) %>%
  mutate(values_centred = values - mean(values)) %>%
  ungroup()

# Anova of full matrix
mdl <- lm(values_centred ~ genotype * treatment, data = long_df)
A <- Anova(mdl, type = "II")

# Correct the residual df: the gene means use up one df per gene
n_genes <- n_distinct(long_df$`AGI MAPPING`)
df_res  <- nrow(long_df) - n_genes - 3     # 3 = genotype, treatment, interaction
ss_res  <- A["Residuals", "Sum Sq"]
terms   <- rownames(A) != "Residuals"

A[terms, "F value"] <- (A[terms, "Sum Sq"] / A[terms, "Df"]) / (ss_res / df_res)
A[terms, "Pr(>F)"]  <- pf(A[terms, "F value"], A[terms, "Df"], df_res, lower.tail = FALSE)
A["Residuals", "Df"] <- df_res
A

write.csv(A, file.path(out_dir, "full_matrix_anova.csv"))