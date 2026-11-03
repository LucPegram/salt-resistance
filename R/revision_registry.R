revision_registry <- function() {
  genotype <- "microarray_genotype"
  dose <- c("microarray_root", "microarray_shoot")
  protein <- c("proteomics_root", "proteomics_shoot")
  legacy <- c("saved_proteomics_root", "saved_proteomics_shoot")
  entry <- function(inputs,builder,filename,status="revised analysis")
    list(inputs=inputs,builder=builder,filename=filename,status=status)
  entries <- list(
    table_revision_genotype = entry(genotype,"build_revision_microarray","genotype_contrasts.csv"),
    table_revision_dose = entry(dose,"build_revision_dose","dose_joint_tests.csv"),
    table_revision_protein = entry(protein,"build_revision_protein","protein_contrasts.csv"),
    table_revision_summary = entry(c(genotype,protein),"build_revision_summary","discovery_summary.csv"),
    table_revision_candidates = entry(protein,"build_revision_candidates","candidate_contrasts.csv"),
    table_revision_missingness = entry(protein,"build_revision_missingness","condition_missingness.csv"),
    table_revision_identity = entry(protein,"build_revision_identity","candidate_identity.csv"),
    table_revision_prestress = entry(c(genotype,protein),"build_revision_prestress","prestress_similarity.csv","exploratory resemblance"),
    table_revision_qc = entry(c(dose,genotype,protein),"build_revision_qc","sample_qc.csv"),
    table_revision_model_qc = entry(c(dose,protein),"build_revision_model_qc","model_diagnostics.csv"),
    table_revision_primary_summary = entry(c(genotype,protein),"build_revision_primary_summary","primary_summary.csv"),
    table_revision_primary_candidates = entry(protein,"build_revision_primary_candidates","primary_candidates.csv"),
    table_revision_similarity_summary = entry(c(genotype,protein),"build_revision_similarity_summary","similarity_summary.csv"),
    table_revision_sensitivity_summary = entry(protein,"build_revision_sensitivity_summary","sensitivity_summary.csv"),
    table_revision_legacy_sensitivity = entry(c(protein,legacy),"build_revision_legacy_sensitivity","legacy_sensitivity.csv"),
    table_revision_enrichment = entry(c(protein,legacy,"go_associations"),"build_revision_enrichment","directional_GO.csv"),
    figure_revision_qc = entry(c(dose,genotype,protein),"build_revision_qc_figure","sample_PCA.pdf","QC"),
    figure_revision_dose = entry(dose,"build_revision_dose_figure","candidate_dose_trajectories.pdf"),
    figure_revision_protein = entry(protein,"build_revision_protein_figure","candidate_observed_abundance.pdf"),
    figure_revision_contrasts = entry(protein,"build_revision_contrast_figure","candidate_contrast_intervals.pdf"),
    figure_revision_prestress = entry(c(genotype,protein),"build_revision_prestress_figure","prestress_similarity.pdf","exploratory resemblance")
  )
  for (id in names(entries)[startsWith(names(entries), "figure_")])
    entries[[id]]$formats <- c("pdf", "svg", "tiff", "png")
  entries
}
