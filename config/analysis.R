# Settings that affect reproducible outputs. No measurements belong in this file.
analysis_settings <- list(
  seed = 20261103L,
  toprank = list(array1_endpoint_mM=100L,factorial_contrast="salt_WT",BH_cutoff=.05,
    table_per_direction=20L,figure_per_direction=10L),
  proteomics_min_observed_in_one_condition = 2L,
  proteomics_min_sequences = 2L,
  genotype_levels = c("No-0", "NIT2"),
  salt_levels_mM = c(0L, 75L),
  tissue_levels = c("Root", "Shoot"),
  qpcr_genotype_labels = c("No-0", "35S::Nit2"),
  qpcr_treatment_labels = c("Untreated", "75 mM NaCl")
)

