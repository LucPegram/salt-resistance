integration_registry <- function() {
  base<-c("microarray_genotype","microarray_root","microarray_shoot","proteomics_root",
    "proteomics_shoot","saved_proteomics_root","saved_proteomics_shoot")
  entry<-function(builder,filename,go=FALSE,figure=FALSE){
    z<-list(inputs=c(base,if(go)"go_associations"),builder=builder,filename=filename,
      status="cross-experiment molecular integration")
    if(figure)z$formats<-c("pdf","svg","tiff","png")
    z
  }
  list(
    table_integration_mapping=entry("build_integration_mapping","mapping_audit.csv"),
    table_integration_genes=entry("build_integration_genes","gene_associations.csv"),
    table_integration_summary=entry("build_integration_summary","integration_summary.csv"),
    table_integration_candidates=entry("build_integration_candidates","integration_candidates.csv"),
    table_integration_pathways=entry("build_integration_pathways","joint_pathways.csv",go=TRUE),
    table_integration_memberships=entry("build_integration_memberships","pathway_memberships.csv",go=TRUE),
    table_integration_influence=entry("build_reviewer_influence","correlation_influence.csv"),
    table_integration_joint_sensitivity=entry("build_reviewer_joint_sensitivity","joint_hit_sensitivities.csv"),
    table_integration_measurement_audit=entry("build_reviewer_measurement_audit","highlighted_measurements.csv"),
    table_integration_candidate_dose=entry("build_reviewer_candidate_dose","candidate_dose_tests.csv"),
    figure_integration_concordance=entry("build_integration_concordance_figure","transcript_protein_effects.pdf",figure=TRUE),
    figure_integration_pathways=entry("build_integration_pathway_figure","pathway_evidence.pdf",go=TRUE,figure=TRUE)
  )
}
