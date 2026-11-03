# Check normality of residuals for anova test.
# 2026/04/15
# Luc Pegram
#------------------------ 

rm(list=ls())


library(car)
library(tidyr)


#===========================
## Shoot Tissue
Shoot = T; Root = F 
#===========================
## Root Tissue
Root = T;  Shoot = F
#===========================


if (Shoot) {tissue <- "Shoot"} else {tissue <- "Root"}



input_file <- paste0("D:/Masters_PhD/Data/Proteomics/2026/DataFiles/R_Exports/",
                     tissue,
                     "Tissue_pg_logFC(min).csv")


og_df <- read.csv(input_file, header = T, check.names = F)

df <- og_df
df[, c(8:19)] <- log(1 + og_df[ , 8:19])  #log transformation


# 2. Prepare data set----------------------------------------

#sample_cols <- colnames(df)[colnames(df) != protein_col]
sample_cols <- colnames(df)[8:19]

long_df <- pivot_longer(
  df,
  cols = all_of(sample_cols),
  names_to = "sample",
  values_to = "values")

long_df <- separate_wider_delim(
  long_df,
  cols = "sample",
  names = c("genotype", "treatment", "sample"),
  delim = "_")

head(long_df, 6)


long_df$values    <- as.numeric(long_df$values)
long_df$genotype  <- as.factor(long_df$genotype)
long_df$treatment <- as.factor(long_df$treatment)


# Anova of full matrix
mdl <- lm(values ~ genotype*treatment, data = long_df)
A <- Anova(mdl, type = "II"); A

r <- residuals(mdl)
#qqPlot(r) #qqnorm(r); qqline(r)

hist(r, freq = FALSE, breaks = 20, xlim = c(min(r), max(r)))
s = seq(min(r), max(r), 1/2)
lines(dnorm(s, sd = sqrt(var(r))) ~ s, lty = 3, lwd = 2.5, col = 'dark red')

