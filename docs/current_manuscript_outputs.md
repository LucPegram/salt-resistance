# Reproducing the current manuscript displays

Statistical output IDs are stable even when manuscript numbering changes.
Manuscript text, measurements and generated files remain outside this code-only
repository. This crosswalk describes the BMC revision with 4 main tables,
6 main figures, 25 supplementary tables and 9 supplementary figures.

Use `scripts/reproduce.R <ID>` for primary output IDs, with `SALT_DATA_ROOT`
pointing to regenerated quantile inputs and `SALT_OUTPUT_ROOT` outside the code.
See `docs/bmc_revision.md` for new audit entry points and external-root variables;
`scripts/reproduce_extensions.R <ID>` for extension outputs; and
`scripts/reproduce_extension_manuscript.R <ID>` for older formatted summaries.
Every script also accepts its documented listing option where available.
Captions are separate from figure artwork. Each output has a provenance record.

| Manuscript table | R output or complete source |
| --- | --- |
| Table 1 | `table_manuscript_design` / experimental_design.csv |
| Table 2 | `table_revision_primary_summary` / primary_summary.csv |
| Table 3 | `table_revision_primary_candidates` / primary_candidates.csv |
| Table 4 | `table_revision_similarity_summary` / similarity_summary.csv |
| S1 | `reproduce_bmc_manuscript_tables.R table_probe_audit` |
| S2 | `reproduce_bmc_manuscript_tables.R table_pair_influence` |
| S3 | `reproduce_extensions.R dose_discoveries` |
| S4 | `table_integration_candidate_dose` / candidate_dose_tests.csv |
| S5 | `table_revision_identity` / candidate_identity.csv |
| S6 | `table_revision_missingness` / condition_missingness.csv |
| S7 | `table_revision_sensitivity_summary` / sensitivity_summary.csv |
| S8 | `table_manuscript_qpcr_full` / qpcr_full_panel.csv |
| S9 | `table_integration_summary` / integration_summary.csv, primary rows |
| S10 | `table_integration_genes` / gene_associations.csv, primary global-BH findings |
| S11 | `table_integration_joint_sensitivity` / joint_hit_sensitivities.csv |
| S12 | `table_integration_measurement_audit` / highlighted_measurements.csv |
| S13 | `table_integration_influence` / correlation_influence.csv |
| S14 | `table_integration_pathways` / joint_pathways.csv, fixed memberships/coverage |
| S15 | `table_integration_candidates` / integration_candidates.csv, basal rows |
| S16 | `reproduce_bmc_manuscript_tables.R table_supporting_candidates` |
| S17 | `reproduce_extension_manuscript.R table_external_baseline` |
| S18 | `reproduce_bmc_manuscript_tables.R table_strict_probe_baseline` |
| S19 | `reproduce_score_inference_sensitivity.R` / table_basal_score_inference.csv |
| S20 | `reproduce_bmc_manuscript_tables.R table_TF_correlations` |
| S21 | `reproduce_extension_manuscript.R table_marker_baseline`; ordinary sensitivity in marker_score_inference_sensitivity.csv |
| S22 | `table_existing_microarray_results` / existing_microarray_results.csv |
| S23 | `table_existing_proteomics_results` / existing_proteomics_results.csv |
| S24 | `table_existing_pairwise_results` / existing_pairwise_results.csv |
| S25 | `table_existing_go_inventory` / existing_go_inventory.csv; current full directional results in `table_revision_enrichment` |

| Manuscript figure | R output |
| --- | --- |
| Figure 1 | `figure_revision_dose` / candidate_dose_trajectories.pdf |
| Figure 2 | `figure_revision_protein` / candidate_observed_abundance.pdf |
| Figure 3 | `figure_revision_contrasts` / candidate_contrast_intervals.pdf |
| Figure 4 | `figure_revision_qc` / sample_PCA.pdf |
| Figure 5 | `figure_integration_concordance` / transcript_protein_effects.pdf |
| Figure 6 | `reproduce_score_inference_sensitivity.R figure_external_score_inference` |
| S1 | `reproduce_bmc_robustness.R figure_dose_omission` |
| S2 | `reproduce_extensions.R figure_dose_programs` |
| S3 | `reproduce_extensions.R figure_dose_candidates` |
| S4 | `figure_manuscript_qpcr_candidates` / qpcr_candidates.pdf |
| S5 | `figure_integration_pathways` / pathway_evidence.pdf |
| S6 | `reproduce_extensions.R figure_external_signatures` |
| S7 | `reproduce_extensions.R figure_external_concordance` |
| S8 | `reproduce_bmc_manuscript_tables.R figure_TF_correlation_sensitivity`; all settings and global BY in S20 |
| S9 | `reproduce_extensions.R figure_marker_signatures` |

The full eligible testing families, excluded features and nonsignificant results
remain in the complete CSV outputs. Display subsets retain their source-family
adjustments; no multiple-testing correction is recalculated only on displayed
hits. Gene-level dose analyses, transcript-ID display analyses and protein-group
analyses have explicitly different units. Strict probe selection, harvested-pair
omission, TF correlation and score-variance analyses are sensitivities, not
additional independent replications.
