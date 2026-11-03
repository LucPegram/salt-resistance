# Descriptive IDs are used until the final manuscript labels are agreed.
artifact_registry <- list(
  table_input_inventory = list(
    inputs = names(input_specs)[vapply(input_specs, function(x)
      x$assay %in% c("microarray", "proteomics", "qpcr_summary"), logical(1))], builder = "build_input_inventory",
    filename = "input_inventory.csv", status = "audit"),
  table_microarray_mapping_audit = list(
    inputs = c("microarray_root", "microarray_shoot", "microarray_genotype"),
    builder = "build_mapping_audit", filename = "microarray_mapping_audit.csv",
    status = "audit"),
  table_proteomics_filtering = list(
    inputs = c("proteomics_root", "proteomics_shoot"),
    builder = "build_proteomics_filtering", filename = "proteomics_filtering.csv",
    status = "audit"),
  figure_proteomics_missingness = list(
    inputs = c("proteomics_root", "proteomics_shoot"),
    builder = "build_proteomics_missingness", filename = "proteomics_missingness.pdf",
    status = "qc"),
  table_qpcr_reported_expression = list(
    inputs = "qpcr_summary", builder = "build_qpcr_table",
    filename = "qpcr_reported_expression.csv", status = "descriptive"),
  figure_qpcr_reported_expression = list(
    inputs = "qpcr_summary", builder = "build_qpcr_figure",
    filename = "qpcr_reported_expression.pdf", status = "descriptive")
)
artifact_registry <- c(artifact_registry, list(
  table_manuscript_design = list(
    inputs = c("microarray_root", "microarray_shoot", "microarray_genotype", "proteomics_root", "proteomics_shoot", "qpcr_summary"),
    builder = "build_manuscript_design", filename = "experimental_design.csv", status = "design metadata"),
  table_manuscript_qpcr_candidates = list(
    inputs = "qpcr_summary", builder = "build_manuscript_qpcr_candidates",
    filename = "qpcr_candidates.csv", status = "descriptive"),
  table_manuscript_qpcr_full = list(
    inputs = "qpcr_summary", builder = "build_manuscript_qpcr_full",
    filename = "qpcr_full_panel.csv", status = "descriptive"),
  figure_manuscript_qpcr_candidates = list(
    inputs = "qpcr_summary", builder = "build_manuscript_qpcr_figure",
    filename = "qpcr_candidates.pdf", status = "descriptive"),
  table_existing_microarray_results = list(
    inputs = names(input_specs)[grepl("^saved_(dose|genotype)", names(input_specs))],
    builder = "build_existing_microarray_results", filename = "existing_microarray_results.csv", status = "legacy summary"),
  table_existing_pairwise_results = list(
    inputs = names(input_specs)[startsWith(names(input_specs), "saved_pairwise_")],
    builder = "build_existing_pairwise_results", filename = "existing_pairwise_results.csv", status = "legacy summary"),
  table_existing_proteomics_results = list(
    inputs = c("saved_proteomics_root", "saved_proteomics_shoot"),
    builder = "build_existing_proteomics_results", filename = "existing_proteomics_results.csv", status = "legacy summary"),
  table_proteomics_completecase_audit = list(
    inputs = c("proteomics_root", "proteomics_shoot", "saved_proteomics_root", "saved_proteomics_shoot"),
    builder = "build_completecase_audit", filename = "proteomics_completecase_audit.csv", status = "numerical audit"),
  table_existing_candidate_proteins = list(
    inputs = c("saved_proteomics_root", "saved_proteomics_shoot"),
    builder = "build_existing_candidates", filename = "existing_candidate_proteins.csv", status = "legacy summary"),
  table_manuscript_candidate_proteins = list(
    inputs = c("saved_proteomics_root", "saved_proteomics_shoot"),
    builder = "build_manuscript_candidates", filename = "candidate_protein_table.csv", status = "legacy summary"),
  table_existing_go_inventory = list(
    inputs = names(input_specs)[startsWith(names(input_specs), "saved_go_")],
    builder = "build_existing_go_inventory", filename = "existing_go_inventory.csv", status = "legacy summary"),
  table_existing_go_terms = list(
    inputs = names(input_specs)[startsWith(names(input_specs), "saved_go_")],
    builder = "build_existing_go_terms", filename = "existing_go_terms.csv", status = "legacy summary"),
  table_existing_phosphopeptides = list(
    inputs = c(names(input_specs)[startsWith(names(input_specs), "saved_phospho_")], "aha2_reference"),
    builder = "build_existing_phosphopeptides", filename = "existing_phosphopeptides.csv", status = "selected subset"),
  figure_existing_dose_response_counts = list(
    inputs = names(input_specs)[startsWith(names(input_specs), "saved_dose_")],
    builder = "build_existing_dose_figure", filename = "existing_dose_response_counts.pdf", status = "legacy summary"),
  figure_existing_proteomic_effect_counts = list(
    inputs = c("saved_proteomics_root", "saved_proteomics_shoot"),
    builder = "build_existing_proteomics_figure", filename = "existing_proteomic_effect_counts.pdf", status = "legacy summary")
))

artifact_registry <- c(artifact_registry, revision_registry())

