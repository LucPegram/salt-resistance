# Settings that affect reproducible outputs. No measurements belong in this file.
analysis_settings <- list(
  seed = 20261103L,
  proteomics_min_observed_in_one_condition = 2L,
  proteomics_min_sequences = 2L,
  genotype_levels = c("No-0", "NIT2"),
  salt_levels_mM = c(0L, 75L),
  tissue_levels = c("Root", "Shoot"),
  qpcr_genotype_labels = c("No-0", "35S::Nit2"),
  qpcr_treatment_labels = c("Untreated", "75 mM NaCl")
)

