# Manuscript output commands

Set external input and output roots as explained in README. Each command works
without another output directory or a pre-existing model cache. It may fit the
required models, including slower dropout sensitivities, before writing the output.

| Manuscript item | Command |
|---|---|
| Main Table 1 | `Rscript --vanilla scripts/reproduce.R table_manuscript_design` |
| Main Table 2 | `Rscript --vanilla scripts/reproduce.R table_revision_primary_summary` |
| Main Table 3 | `Rscript --vanilla scripts/reproduce.R table_revision_primary_candidates` |
| Main Table 4 | `Rscript --vanilla scripts/reproduce.R table_revision_similarity_summary` |
| Main Figure 1 | `Rscript --vanilla scripts/reproduce.R figure_revision_dose` |
| Main Figure 2 | `Rscript --vanilla scripts/reproduce.R figure_revision_protein` |
| Main Figure 3 | `Rscript --vanilla scripts/reproduce.R figure_revision_contrasts` |
| Main Figure 5 | `Rscript --vanilla scripts/reproduce.R figure_integration_concordance` |
| Main Figure 4 | `Rscript --vanilla scripts/reproduce.R figure_revision_qc` |
| Supplementary Table S1 | `Rscript --vanilla scripts/reproduce.R table_existing_microarray_results` |
| Supplementary Table S2 | `Rscript --vanilla scripts/reproduce.R table_existing_proteomics_results` |
| Supplementary Table S3 | `Rscript --vanilla scripts/reproduce.R table_existing_pairwise_results` |
| Supplementary Table S4 historical inventory | `Rscript --vanilla scripts/reproduce.R table_existing_go_inventory` |
| Supplementary Table S4 complete revised enrichment | `Rscript --vanilla scripts/reproduce.R table_revision_enrichment` |
| Supplementary Table S5 | `Rscript --vanilla scripts/reproduce.R table_revision_identity` |
| Supplementary Table S6 | `Rscript --vanilla scripts/reproduce.R table_revision_missingness` |
| Supplementary Table S7 | `Rscript --vanilla scripts/reproduce.R table_revision_sensitivity_summary` |
| Supplementary Table S8 | `Rscript --vanilla scripts/reproduce.R table_manuscript_qpcr_full` |
| Supplementary Table S9 | `Rscript --vanilla scripts/reproduce.R table_integration_summary` |
| Supplementary Table S10 | `Rscript --vanilla scripts/reproduce.R table_integration_genes` |
| Supplementary Table S11 | `Rscript --vanilla scripts/reproduce.R table_integration_pathways` |
| Supplementary Table S11 memberships | `Rscript --vanilla scripts/reproduce.R table_integration_memberships` |
| Supplementary Table S13 | `Rscript --vanilla scripts/reproduce.R table_integration_candidate_dose` |
| Supplementary Table S14 | `Rscript --vanilla scripts/reproduce.R table_integration_influence` |
| Supplementary Table S15 | `Rscript --vanilla scripts/reproduce.R table_integration_joint_sensitivity` |
| Supplementary Table S16 | `Rscript --vanilla scripts/reproduce.R table_integration_measurement_audit` |
| Supplementary Table S12 | `Rscript --vanilla scripts/reproduce.R table_integration_candidates` |
| Supplementary mapping audit | `Rscript --vanilla scripts/reproduce.R table_integration_mapping` |
| Supplementary Figure S2 | `Rscript --vanilla scripts/reproduce.R figure_integration_pathways` |
| Supplementary Figure S1 | `Rscript --vanilla scripts/reproduce.R figure_manuscript_qpcr_candidates` |

Numeric CSV columns are formatted for Word presentation without changing their
values. Manuscript Table 3 explicitly includes ineligible candidate rows, whose
estimates are unavailable. Interval and correction labels must remain attached.

## Complete result tables and audits

| Result | Output ID |
|---|---|
| Every genotype-array contrast | `table_revision_genotype` |
| Every paired dose joint test | `table_revision_dose` |
| Every protein contrast and sensitivity | `table_revision_protein` |
| All method/contrast discovery summaries | `table_revision_summary` |
| Candidate contrasts in all methods | `table_revision_candidates` |
| Every control-split similarity calculation | `table_revision_prestress` |
| Sample distributions and original PCA coordinates | `table_revision_qc` |
| Figure 4 coordinates and treatment/genotype labels | `table_manuscript_pca` |
| Pair correlation, dropout convergence and warnings | `table_revision_model_qc` |
| Legacy/revised discovery overlap | `table_revision_legacy_sensitivity` |
| Selected phosphopeptide identity and canonical mapping | `table_existing_phosphopeptides` |
| Complete-case ordinary-ANOVA numerical check | `table_proteomics_completecase_audit` |

Use the same single-output command with any registered ID above.
The additional requested three top-20 tables and three top-10 treatment figures
are documented in [top_rankings.md](top_rankings.md). They use distinct gene-locus
rankings and explicit ranking-specific correction families; they do not replace
the original manuscript tables or statistics. Run
`Rscript --vanilla scripts/reproduce_toprank.R` to rebuild all fifteen associated
table, figure, audit and signal artifacts.
The seven current manuscript figure commands also export SVG, 600-dpi TIFF and
PNG previews, each recorded in provenance. To rebuild just these five non-integration figures,
run `Rscript --vanilla scripts/reproduce_figures.R`; rebuild the ten integration and reviewer-diagnostic
tables and two additional figures with `Rscript --vanilla scripts/reproduce_integration.R`. See `docs/figure_design.md`.
`Rscript --vanilla scripts/inventory_sources.R` creates source and historical
script inventories locally. `legacy/README.md` lists original script hashes.
No inventory of result values is stored in this repository.

Figure 4 now shows sample PCA for all treatment/genotype combinations. The earlier
`figure_revision_prestress` remains a reproducible historical output and is not
part of the current manuscript. PCA is no longer duplicated in the supplement;
qPCR and pathways are now Supplementary Figures S1 and S2.

Additional candidate prioritization of significant, same-direction transcript/protein
trends is documented in [supporting_trends.md](supporting_trends.md). These eight
outputs supplement the manuscript analyses and do not replace their correction
families or primary joint-evidence criteria.
