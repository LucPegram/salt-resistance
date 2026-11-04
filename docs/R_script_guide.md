# Guide to every R script

This guide explains the purpose, main operations, inputs and outputs of every R
script in this repository. It covers **60 current workflow, configuration, maintenance and test files**.
Descriptions cover the source shipped with this guide, including the dataset-ranking
module added on 3 November 2026 and the eight computational-extension files added on 4 November 2026. Original module descriptions were checked against
[`c1c3847`](https://github.com/LucPegram/salt-resistance/tree/c1c3847f630b87c36ee3d0f9f7adc0974460c869);
[dataset-ranking methods](top_rankings.md) explain the additional analysis.

The current assays are **gene-expression microarrays**, not mature-microRNA
measurements. Statistical replication is at the independently grown plant-pool
level. Root/shoot pairing in the dose experiment does not establish pairing
between genotype-treatment cells or between assays.

## How the current workflow fits together

1. Run the command-line entry points in `scripts/`; they source `bootstrap.R`.
2. `config/inputs.R` names external inputs, `config/analysis.R` holds settings,
   and `config/artifacts.R` combines the output registries.
3. `R/inputs.R` validates measurements. The model modules fit revised analyses;
   historical-audit modules read saved results for explicitly labelled summaries.
4. Table/figure builders write one registered artifact at a time through
   `R/reproduce.R`, recording hashes, versions and a completion marker.
5. Keep data, caches, generated outputs, manuscripts and run records local,
   outside the tracked code/documentation files.

A **builder** is a function that writes a registered table or figure. Most accept
`inputs`, `settings` and `destination`. An **artifact ID** selects a registry
entry, not a script filename. Function lists below include named local helpers;
those helpers belong to their enclosing function and are not standalone APIs.
Files that contain only top-level commands are identified explicitly.

Inactive historical source and its importer were removed during cleanup. They
remain recoverable from Git history; see [cleanup and recovery](cleanup.md).

## Starting commands

Run from the repository root after restoring the locked R environment. Set
`SALT_DATA_ROOT` and `SALT_OUTPUT_ROOT` to separate external input/output locations;
`SALT_R_LIB` and `SALT_MODEL_CACHE` are optional as described in [README.md](../README.md).

```sh
Rscript --vanilla scripts/check_environment.R
Rscript --vanilla scripts/reproduce.R --list
Rscript --vanilla scripts/reproduce.R table_revision_primary_candidates
Rscript --vanilla scripts/reproduce_figures.R
Rscript --vanilla scripts/reproduce_integration.R
```

`reproduce_figures.R` builds five current non-integration manuscript figures;
`reproduce_integration.R` adds the two integration figures and integration tables.
For each final manuscript table/figure and its exact command, see
[current_manuscript_outputs.md](current_manuscript_outputs.md). Model definitions and qualifications
are explained in [analysis_decisions.md](analysis_decisions.md) and
[integration_analysis.md](integration_analysis.md).

## Script catalogue

## Command-line scripts and maintenance

### [`scripts/bootstrap.R`](../scripts/bootstrap.R)

Initializes the current workflow. Finds the repository from the running Rscript path, uses UTC
and C collation, optionally prepends the library named by SALT_R_LIB, then sources the R
modules and configuration.

**Inputs and outputs:** Reads repository source/configuration; creates shared objects in the R
session. Writes no analysis files.

**Operations and interpretation:** This is sourced by the command-line scripts and tests. It
does not source legacy scripts or install packages. The modules define functions before
configuration is loaded; analyses run when a builder is called.

**Structure:** Top-level commands; no named functions are defined.

### [`scripts/check_code_only.R`](../scripts/check_code_only.R)

Checks whether repository paths comply with the code-and-documentation-only publication
policy.

**Inputs and outputs:** Checks paths from an optional text-file argument, otherwise tracked
Git files, or a directory listing when Git metadata is absent. Prints success or stops with
prohibited paths.

**Operations and interpretation:** Rejects data/results/figures/tables directories,
manuscripts and other extensions, local configuration and session/library files. The default
Git mode checks tracked files, so stage intended additions before using it as a publication
check. It is a path/extension allowlist, not a scan for measurements embedded in allowed text
files.

**Structure:** Top-level commands; no named functions are defined.

### [`scripts/check_environment.R`](../scripts/check_environment.R)

Compares the running R version and every package version in renv.lock with the installed
environment.

**Inputs and outputs:** Reads renv.lock and installed package metadata; prints success or
stops with version/missing-package problems.

**Operations and interpretation:** Does not install dependencies or prove that a clean-machine
restoration works. Run after renv::restore() and before analyses.

**Structure:** Top-level commands; no named functions are defined.

### [`scripts/inventory_sources.R`](../scripts/inventory_sources.R)

Inventories external source files and inspects supplied R scripts without executing them.

**Inputs and outputs:** Reads SALT_DATA_ROOT and writes source_inventory.csv and
legacy_analysis_inventory.csv under SALT_OUTPUT_ROOT.

**Operations and interpretation:** Records relative paths, sizes, SHA256 hashes, heuristic
source roles, parse status, referenced packages and Windows-path flags. File presence or a
successful parse is not evidence that a historical analysis was executed.

**Structure:** Top-level commands; no named functions are defined.

### [`scripts/reproduce_figures.R`](../scripts/reproduce_figures.R)

Rebuilds five current non-integration manuscript figure artifacts: dose trajectories,
observed proteins, contrast intervals, treatment/genotype PCA and candidate qPCR.
The former pre-stress graphic is historical and is not rebuilt by this entry point.

**Inputs and outputs:** Requires SALT_DATA_ROOT and SALT_OUTPUT_ROOT. Produces PDF, SVG, TIFF
and PNG for each of these five figures, with provenance records.

**Operations and interpretation:** The local inside() helper checks directory containment.
This script does not build the two additional integration figures; use reproduce_integration.R
for those, or reproduce.R with their IDs.

**Named functions, including local helpers:** `inside()`.

### [`scripts/reproduce_integration.R`](../scripts/reproduce_integration.R)

Rebuilds all twelve registered integration artifacts, including the four diagnostics added
after visual review.

**Inputs and outputs:** Requires external input/output roots; writes ten table artifacts and
two figure artifacts with provenance.

**Operations and interpretation:** Uses integration_registry() to determine the outputs, the
local inside() helper to enforce directory separation, and run_artifact() to build each
independently. Includes the two integration figures needed to complete the seven manuscript
figures.

**Named functions, including local helpers:** `inside()`.

### [`scripts/reproduce.R`](../scripts/reproduce.R)

Main command-line entry point for listing or rebuilding registered outputs. Accepts one
artifact ID, --list, or all.

**Inputs and outputs:** Reads external files selected by the registry and writes one
subdirectory per artifact under SALT_OUTPUT_ROOT; the default is the ignored outputs/
directory.

**Operations and interpretation:** Checks that input and output roots are separate and that
experimental inputs are outside the repository. Its local inside() helper checks directory
containment. Calls run_artifact() with a shared import cache. The all option includes
historical summaries and audits as well as revised analyses.

**Named functions, including local helpers:** `inside()`.

### [`scripts/snapshot_dependencies.R`](../scripts/snapshot_dependencies.R)

Maintainer utility that records the installed direct and transitive dependency versions in
renv.lock.

**Inputs and outputs:** Reads installed package metadata and writes renv.lock, including R,
CRAN and Bioconductor source records.

**Operations and interpretation:** Use only for an intentional dependency update. This changes
the environment specification; it does not restore packages or run analyses.

**Structure:** Top-level commands; no named functions are defined.

## Configuration

### [`config/analysis.R`](../config/analysis.R)

Defines analysis_settings: random seed, original protein-retention thresholds, genotype/tissue
levels, qPCR group labels and dataset-ranking settings.

**Inputs and outputs:** Creates an R list in memory; contains settings rather than
observations.

**Operations and interpretation:** Changing relevant settings can change outputs and cache
fingerprints. Protein retention here means enough observations in at least one cell plus
sequence support; the stricter primary model eligibility is implemented in revision_models.R.

**Structure:** Top-level commands; no named functions are defined.

### [`config/artifacts.R`](../config/artifacts.R)

Assembles artifact_registry, the dispatch catalogue for independently reproducible outputs.

**Inputs and outputs:** Each entry specifies input IDs, a builder function, destination
filename, status and, where appropriate, graphic formats.

**Operations and interpretation:** Defines the initial import/QC, qPCR and historical-audit
entries, then appends revision_registry(), integration_registry() and toprank_registry().
There are 78 artifact IDs, including fifteen dataset-ranking outputs. This is
registration/configuration code, not a model fit.

**Structure:** Top-level commands; no named functions are defined.

### [`config/inputs.R`](../config/inputs.R)

Defines input_specs, mapping logical input IDs to external relative paths, assay types and
relevant metadata. Adds historical result specifications, the AHA2 FASTA reference, GO
associations and explicitly enumerated historical GO files.

**Inputs and outputs:** Paths are resolved against SALT_DATA_ROOT by read_input(); this file
reads no measurements itself.

**Operations and interpretation:** Also declares the supplied gene-name annotation input for
dataset rankings. Calls legacy_result_specs() from manuscript_evidence.R.
Merely declaring a path does not require its file for every artifact: each registry entry
specifies which inputs it needs.

**Structure:** Top-level commands; no named functions are defined.

## Active analysis and output modules

### [`R/artifacts.R`](../R/artifacts.R)

Builds initial input audits, missingness plots and descriptive qPCR outputs.

**Inputs and outputs:** Uses validated import objects to write input_inventory.csv,
microarray_mapping_audit.csv, proteomics_filtering.csv, proteomics_missingness.pdf,
qpcr_reported_expression.csv and qpcr_reported_expression.pdf.

**Operations and interpretation:** write_table() exports CSV with explicit NA.
build_input_inventory(), build_mapping_audit() and build_proteomics_filtering() report
row/ID/filter counts. build_proteomics_missingness() plots missing measurements before
filtering/imputation. build_qpcr_table() copies supplied summaries; build_qpcr_figure() plots
ordered group means and available STE using a local key() helper. These are descriptive/audit
outputs, not new qPCR inference or the redesigned manuscript qPCR figure.

**Named functions, including local helpers:** `write_table()`, `build_input_inventory()`,
`build_mapping_audit()`, `build_proteomics_filtering()`, `build_proteomics_missingness()`,
`build_qpcr_table()`, `build_qpcr_figure()`, `key()`.

### [`R/inputs.R`](../R/inputs.R)

Validates and imports microarray, protein and qPCR summary inputs while preserving assay
identities and missing observations.

**Inputs and outputs:** Reads external CSV/TSV/XLSX inputs and returns structured R lists
containing measurements, annotations and sample metadata; writes no files.

**Operations and interpretation:** assert_columns() checks schemas; as_numeric_checked()
validates numeric values. read_microarray() checks sample names, positive signals and
corrected transcript IDs; collapse_microarray() averages repeated exact transcript IDs on the
log2 scale and retains transcript suffixes. read_proteomics() preserves full protein-group
identities and constructs metadata; proteomics_keep() applies the original retention rule.
read_qpcr_summary() checks group means/STE without inventing replication. read_input()
dispatches by assay, including historical, GO and gene-name annotation readers defined in
other modules.

**Named functions, including local helpers:** `assert_columns()`, `as_numeric_checked()`,
`read_microarray()`, `collapse_microarray()`, `read_proteomics()`, `proteomics_keep()`,
`read_qpcr_summary()`, `read_input()`.

### [`R/integration_figures.R`](../R/integration_figures.R)

Draws the two integration figures using the shared publication style.

**Inputs and outputs:** Produces transcript_protein_effects and pathway_evidence in PDF, SVG,
TIFF and PNG.

**Operations and interpretation:** build_integration_concordance_figure() plots all eligible
effects for basal genotype, No-0 salt and Col-0 context in each protein tissue, highlighting
candidates and joint-BH evidence. build_integration_pathway_figure() shows all nine planned
pathways across twelve families, marking joint support, direction and insufficient coverage.
No cross-assay sample pairing or gene-as-replicate significance test is introduced.

**Named functions, including local helpers:** `build_integration_concordance_figure()`,
`build_integration_pathway_figure()`.

### [`R/integration_models.R`](../R/integration_models.R)

Matches transcript and protein evidence at AGI locus level across experiments and tests
same-direction joint evidence.

**Inputs and outputs:** Uses the three array matrices, two protein matrices and saved protein
AGI annotation workbooks; writes mapping_audit.csv, gene_associations.csv,
integration_summary.csv and integration_candidates.csv.

**Operations and interpretation:** integration_settings() defines contrasts, candidates and
pathway settings; integration_cached() adds integration-source hashes to model caching.
integration_loci() gives each unique transcript ID equal weight within a locus after
repeated-spot collapse. integration_array() refits locus-level models.
integration_protein_mapping() requires a single accession, reported proteotypic support, valid
unique AGI, shared array coverage and primary observation eligibility.
integration_conjunction() computes conservative directional intersection-union P values;
integration_adjust() applies family/global BH and global BY separately by sensitivity.
integration_join() joins one contrast; integration_results() assembles twelve comparison
families with primary, complete-case and single-transcript sensitivities.
integration_summary() calculates descriptive correlations, omission of NIT2 and evidence
counts. The four build_integration_* functions export the corresponding objects. Independent
experiments, unequal salt doses/backgrounds and tissue aggregation remain limitations; genes
are not biological replicates.

**Named functions, including local helpers:** `integration_settings()`,
`integration_cached()`, `integration_loci()`, `integration_array()`,
`integration_protein_mapping()`, `integration_conjunction()`, `integration_adjust()`,
`integration_join()`, `integration_results()`, `build_integration_genes()`,
`build_integration_mapping()`, `integration_summary()`, `build_integration_summary()`,
`build_integration_candidates()`.

### [`R/integration_pathways.R`](../R/integration_pathways.R)

Tests nine focused GO pathways for joint directional transcript/protein enrichment.

**Inputs and outputs:** Adds external GO associations to the integration inputs; writes
joint_pathways.csv and pathway_memberships.csv, including eligible universes and untestable
sets.

**Operations and interpretation:** integration_canonical_associations() resolves secondary GO
IDs using pinned GO.db before ancestor propagation. integration_pathway_inputs() aligns both
assays on a common uniquely mapped, fully observed universe. integration_pathways() runs
CAMERA with estimated nonnegative inter-gene correlation and a separate fixed-0.01
sensitivity; requires at least five members and at most 80% of the universe; combines
directions with the conjunction test and adjusts across families.
build_integration_pathways()/build_integration_memberships() export results and exact
membership. Enrichment direction is not measured transport or ROS activity.

**Named functions, including local helpers:** `integration_pathway_inputs()`,
`integration_canonical_associations()`, `integration_pathways()`,
`build_integration_pathways()`, `build_integration_memberships()`.

### [`R/integration_registry.R`](../R/integration_registry.R)

Defines integration_registry(), the twelve independently reproducible integration/diagnostic
outputs.

**Inputs and outputs:** Returns ten table and two figure registry records; writes no files
itself.

**Operations and interpretation:** The local entry() helper attaches common required inputs,
optional GO associations, builder/filename and figure formats. The candidate-dose diagnostic
shares the integration runner but extracts tests from the original paired dose model.

**Named functions, including local helpers:** `integration_registry()`, `entry()`.

### [`R/manuscript_evidence.R`](../R/manuscript_evidence.R)

Imports and summarizes historical result files, documents the experimental design and qPCR
summaries, and performs narrowly defined numerical/peptide audits.

**Inputs and outputs:** Reads saved CSV/XLSX/GO/FASTA inputs; writes historical summary/count
tables and plots, design/qPCR tables, candidate summaries, the complete-case ANOVA audit and
selected-phosphopeptide audit.

**Operations and interpretation:** legacy_result_specs() declares historical inputs;
read_legacy(), with local read_sheet(), imports them; checked_q() checks saved adjusted P
values. build_existing_microarray_results(), build_existing_pairwise_results() and
build_existing_proteomics_results() summarize supplied results. balanced_protein_anova() and
build_completecase_audit() independently recalculate ordinary ANOVA only for fully observed
proteins. build_existing_candidates() and build_manuscript_candidates() extract historical
candidates. build_existing_go_inventory()/build_existing_go_terms() summarize saved enrichment
reports. build_existing_phosphopeptides() deduplicates supplied selected precursors and calls
map_modified_peptide() to locate modification positions in the supplied reference, without
establishing site occupancy or activation. build_manuscript_design() describes assay designs;
qpcr_summary_table() and its candidate/full builders format descriptive summaries;
build_manuscript_qpcr_figure() delegates to the publication renderer. The two
build_existing_*_figure() count plots display historical discoveries. Summarizing saved
statistics does not validate historical models or replace revised inference.

**Named functions, including local helpers:** `legacy_result_specs()`, `read_legacy()`,
`read_sheet()`, `checked_q()`, `build_existing_microarray_results()`,
`build_manuscript_design()`, `build_manuscript_qpcr_candidates()`,
`build_manuscript_qpcr_full()`, `qpcr_summary_table()`, `build_manuscript_qpcr_figure()`,
`build_existing_pairwise_results()`, `build_existing_proteomics_results()`,
`balanced_protein_anova()`, `build_completecase_audit()`, `build_existing_candidates()`,
`build_manuscript_candidates()`, `build_existing_go_inventory()`, `build_existing_go_terms()`,
`build_existing_phosphopeptides()`, `map_modified_peptide()`, `build_existing_dose_figure()`,
`build_existing_proteomics_figure()`.

### [`R/reproduce.R`](../R/reproduce.R)

Runs one registered artifact and records reproducibility evidence.

**Inputs and outputs:** Reads the entry's required external inputs; writes the artifact,
provenance.R, sessionInfo.txt and RUN_COMPLETE in its output subdirectory.

**Operations and interpretation:** source_files() identifies source/documentation files to
hash. run_artifact() validates registration, caches imports, fixes the seed, calls the named
builder, verifies all expected files, and checks input/source hashes before and after
generation. Its local hash() computes SHA256. Provenance records settings, versions, Git
revision and input/code/output hashes. RUN_COMPLETE is written only after successful
completion. Generated provenance.R files live in local outputs and are not repository analysis
scripts.

**Named functions, including local helpers:** `source_files()`, `run_artifact()`, `hash()`.

### [`R/reviewer_diagnostics.R`](../R/reviewer_diagnostics.R)

Adds the four transparent diagnostics requested during visual review without redefining
primary inference.

**Inputs and outputs:** Writes correlation_influence.csv, joint_hit_sensitivities.csv,
highlighted_measurements.csv and candidate_dose_tests.csv.

**Operations and interpretation:** reviewer_loo_cor() efficiently calculates every
leave-one-locus Pearson correlation. build_reviewer_influence() reports correlations with and
without the post hoc highlighted locus set and every single-locus omission range.
build_reviewer_joint_sensitivity(), with local value(), places complete-case/single-transcript
BH and primary BY beside each primary joint hit. build_reviewer_measurement_audit(), with
local append_cells()/means(), reports observed cell counts/ranges and independently checks
highlighted basal/saline effects against cell means. build_reviewer_candidate_dose() extracts
original root/shoot/tissue-by-dose tests for displayed transcripts with their original
correction universe. These are processed-input/statistical diagnostics, not independent
molecular validation or biological sample exclusions.

**Named functions, including local helpers:** `reviewer_loo_cor()`,
`build_reviewer_influence()`, `build_reviewer_joint_sensitivity()`, `value()`,
`build_reviewer_measurement_audit()`, `append_cells()`, `means()`,
`build_reviewer_candidate_dose()`.

### [`R/revision_enrichment.R`](../R/revision_enrichment.R)

Performs direction-separated protein GO over-representation analysis using an eligible
measured universe.

**Inputs and outputs:** Reads the supplied gene-GO association resource and uses revised
protein effects plus saved AGI annotations; writes directional_GO.csv.

**Operations and interpretation:** read_revision_associations() validates the four-column
resource. revision_go_sets() propagates membership to ancestors in the pinned GO.db.
build_revision_enrichment() excludes ambiguous groups/loci, defines an annotated measured
universe, separates higher/lower queries for three contrasts, calculates upper-tail
hypergeometric P values, and applies family/global BH correction. Integration-specific
secondary-ID resolution is implemented separately in integration_pathways.R; it is not
silently applied to this original analysis.

**Named functions, including local helpers:** `read_revision_associations()`,
`revision_go_sets()`, `build_revision_enrichment()`.

### [`R/revision_figures.R`](../R/revision_figures.R)

Defines shared publication styling and six non-integration figure renderers (five current
outputs plus the historical pre-stress graphic), including
descriptive qPCR.

**Inputs and outputs:** Uses revised models/validated summaries; produces PDF/SVG vectors, 600
dpi TIFF and 200 dpi preview PNG. The model data are not altered for display.

**Operations and interpretation:** publication_figure() handles the four graphic devices;
publication_caption() collects explanatory text in a companion caption file;
publication_range(), publication_axes(), publication_title(), publication_key() and
publication_pool_labels() control ranges, axes, headings, legends and deterministic label
placement. publication_colors stores shared colours. revision_candidate_transcript() chooses
.1 when available, otherwise the first sorted transcript ID, independently of significance.
publication_pca_data() adds dataset, genotype, salt and pool labels to the unchanged
PCA coordinates; build_manuscript_pca_coordinates() exports this sample table.
build_revision_qc_figure() plots the five PCA panels in current Figure 4, using a local
key() helper for explicit treatment/genotype keys; build_revision_dose_figure() records direct tissue-by-dose
q values in its caption companion; build_revision_protein_figure() shows observed pools and cell counts;
build_revision_contrast_figure() shows contrasts and intervals;
build_revision_prestress_figure() leads with the control-sharing sensitivity.
build_publication_qpcr() plots candidate group means/STE using a local key() helper, without
reconstructing replicate-level tests.

**Named functions, including local helpers:** `publication_caption()`, `publication_figure()`, `publication_range()`,
`publication_axes()`, `publication_title()`, `publication_key()`, `publication_pool_labels()`,
`publication_pca_data()`, `build_manuscript_pca_coordinates()`,
`build_revision_qc_figure()`, `revision_candidate_transcript()`,
`build_revision_dose_figure()`, `build_revision_protein_figure()`,
`build_revision_contrast_figure()`, `build_revision_prestress_figure()`,
`build_publication_qpcr()`, `key()`.

### [`R/revision_manuscript_tables.R`](../R/revision_manuscript_tables.R)

Creates compact tables used to assemble the revised manuscript and supplement.

**Inputs and outputs:** Writes primary_summary.csv, primary_candidates.csv,
similarity_summary.csv and sensitivity_summary.csv from revised model objects.

**Operations and interpretation:** build_revision_primary_summary() selects primary discovery
contrasts. build_revision_primary_candidates() reports four contrasts for NIT2, SOS1/NHX7 and
the AHA2-labelled group, retaining ineligible entries and identity qualifications.
build_revision_similarity_summary() reports shared-control correlation and the entire
disjoint-control range. build_revision_sensitivity_summary() compares basal and interaction
discovery counts across protein methods. These tables summarize the fitted models rather than
introduce another fitting procedure.

**Named functions, including local helpers:** `build_revision_primary_summary()`,
`build_revision_primary_candidates()`, `build_revision_similarity_summary()`,
`build_revision_sensitivity_summary()`.

### [`R/revision_models.R`](../R/revision_models.R)

Implements the revised statistical models, missingness sensitivities, sample QC and
exploratory pre-stress comparisons.

**Inputs and outputs:** Consumes validated measurement matrices and, for historical
comparisons only, saved result objects. Returns fitted objects and writes registered contrast,
summary, identity, missingness, diagnostic, QC and resemblance CSVs through its builders.

**Operations and interpretation:** revision_cached() manages hash-keyed memory/optional disk
caches, checks cached payloads and keeps caches outside code and original inputs.
factorial_metadata() constructs four genotype-by-treatment cells; factorial_contrasts()
defines seven explicit contrasts. moderated_contrasts() fits robust, trend-moderated limma
models with pointwise intervals, within-contrast BH and an additional across-contrast
criterion. revision_genotype() fits combined-tissue arrays. revision_dose() fits categorical
tissue-by-dose means, represents root/shoot pool pairing with duplicateCorrelation(), and
tests root, shoot and tissue-by-dose responses. protein_cell_counts() counts observations.
revision_protein() fits the primary observed-data protein model plus available-observation,
complete-case, median-centering and proDA dropout sensitivities. prestress_splits(), with
local summarize(), compares all shared/disjoint control assignments descriptively.
revision_discoveries() summarizes discoveries; revision_qc() calculates distributions,
missingness and centered, unscaled PCA on fully observed features. The build_revision_*
functions export these model/QC objects; build_revision_legacy_sensitivity() compares revised
and historical discovery sets without using old P values as new model inputs.

**Named functions, including local helpers:** `revision_cached()`, `factorial_metadata()`,
`factorial_contrasts()`, `moderated_contrasts()`, `revision_genotype()`, `revision_dose()`,
`protein_cell_counts()`, `revision_protein()`, `prestress_splits()`, `summarize()`,
`build_revision_microarray()`, `build_revision_dose()`, `build_revision_protein()`,
`revision_discoveries()`, `build_revision_summary()`, `build_revision_candidates()`,
`build_revision_missingness()`, `build_revision_identity()`, `build_revision_model_qc()`,
`build_revision_prestress()`, `revision_qc()`, `build_revision_qc()`,
`build_revision_legacy_sensitivity()`.

### [`R/revision_registry.R`](../R/revision_registry.R)

Defines revision_registry(), registering the revised model tables, audits, compact manuscript
tables and five revised figure artifacts.

**Inputs and outputs:** Returns registry entries in memory; writes no files itself.

**Operations and interpretation:** Its local entry() helper constructs
input/builder/filename/status records. Figure records request four graphic formats. The sixth
original manuscript figure, candidate qPCR, is registered in config/artifacts.R.

**Named functions, including local helpers:** `revision_registry()`, `entry()`.

## Synthetic verification

### [`tests/integration_tests.R`](../tests/integration_tests.R)

Runs synthetic integration and reviewer-diagnostic checks.

**Inputs and outputs:** Uses artificial measurements/mappings plus the pinned GO ontology;
prints pass messages or stops on a failed assertion.

**Operations and interpretation:** Checks equal transcript weighting despite unequal duplicate
spots, directional conjunction behaviour, method-specific global BH/BY and family correction,
strict protein-mapping exclusions, secondary GO-ID resolution/de-duplication and
leave-one-locus correlations against direct recomputation. No research measurements are
embedded.

**Structure:** Top-level commands; no named functions are defined.

### [`tests/run_tests.R`](../tests/run_tests.R)

Runs synthetic checks for imports, filtering, registry consistency, peptide mapping and
revised statistical calculations.

**Inputs and outputs:** Uses artificial values and loaded code/configuration; prints pass
messages or stops on a failed assertion. No experimental data are required.

**Operations and interpretation:** expect_error() checks intended failures. Independent
lm()/ANOVA calculations check balanced ANOVA and factorial interaction estimates. Other
assertions check transcript-suffix handling, numeric validation, protein retention, registered
builders, interval ordering, across-contrast correction and exhaustive disjoint-control
partitions. This does not establish biological validity or validate missing acquisition
records.

**Named functions, including local helpers:** `expect_error()`.

## Dataset-ranking additions

### [`R/toprank.R`](../R/toprank.R)

Builds the requested three dataset-specific top-20 up/down tables and top-10
up/down signal-versus-treatment figures, plus complete ranking and signal audits.

**Inputs and outputs:** Uses external array/protein matrices, supplied gene-name
annotations and saved protein AGI mappings. Registers fifteen independently
reproducible CSV, table-PDF and four-format figure artifacts.

**Operations and interpretation:** `toprank_contrast_label()` writes readable
comparison labels. `read_toprank_annotation()`, `toprank_annotation()` and
`toprank_add_names()` import and attach supplied symbols, aliases and functional
names without inventing missing symbols. `toprank_dose()` fits locus-level paired
categorical dose contrasts; `toprank_protein()` selects observed-model protein
statistics with strict identity eligibility. `toprank_select()` corrects before
ranking and chooses one representative tissue per locus. `toprank_result()`
assembles/caches a dataset result; `toprank_signal()` exports plotted observations.
`toprank_table_pdf()` and `toprank_figure()` render the tables and trajectories.
`build_toprank()` dispatches output types, `toprank_registry()` declares artifacts,
and fifteen `build_toprank_<dataset>_<kind>()` wrappers expose independent
builders. The wrappers cover each of `array1`, `array2`, `proteomics` and each of
`table`, `audit`, `signal`, `print`, `figure`. These descriptive rankings have
explicit correction families and do not replace original manuscript inference.

**Named functions:** `toprank_contrast_label()`, `read_toprank_annotation()`,
`toprank_annotation()`, `toprank_add_names()`, `toprank_dose()`, `toprank_protein()`,
`toprank_select()`, `toprank_result()`, `toprank_signal()`, `toprank_table_pdf()`,
`toprank_figure()`, `build_toprank()`, `toprank_registry()`, `build_toprank_array1_table()`,
`build_toprank_array1_audit()`, `build_toprank_array1_signal()`,
`build_toprank_array1_print()`, `build_toprank_array1_figure()`,
`build_toprank_array2_table()`, `build_toprank_array2_audit()`,
`build_toprank_array2_signal()`, `build_toprank_array2_print()`,
`build_toprank_array2_figure()`, `build_toprank_proteomics_table()`,
`build_toprank_proteomics_audit()`, `build_toprank_proteomics_signal()`,
`build_toprank_proteomics_print()`, `build_toprank_proteomics_figure()`.

### [`scripts/reproduce_toprank.R`](../scripts/reproduce_toprank.R)

Rebuilds all fifteen dataset-ranking artifacts using the normal bootstrap,
validated imports and provenance runner.

**Inputs and outputs:** Requires separate external input/output roots and writes
all three tables, table PDFs, treatment figures, full candidate audits and plotted
signal exports. Its named local `inside()` helper checks directory containment.
No data are written into tracked repository paths.

### [`tests/toprank_tests.R`](../tests/toprank_tests.R)

Runs synthetic checks of ranking-specific multiplicity, direction, effect ordering,
distinct-locus selection, symbol/functional-name parsing and paired dose effects.

**Inputs and outputs:** Uses invented values and annotations, printing success or
stopping on failure; no experimental data are used. A local `make()` helper
constructs synthetic root/shoot input objects. Expected dose effects are checked
against independently calculated cell-mean differences.

## Maintaining this guide

Update the relevant entry whenever an R file is added, removed or changes purpose.
Distinguish intended behaviour
in comments from actual executable statements and from independently validated
experimental results.

There is one catalogue heading for each tracked R source file. Generated local
`provenance.R` run records are described under `R/reproduce.R` and are excluded
from the source-file count. This guide contains no experimental measurements or
generated result tables.


## Supporting trend prioritization added after investigator clarification

### [`R/supporting_trends.R`](../R/supporting_trends.R)

Prioritizes significant, same-direction RNA/protein effects using the unchanged
integration fits and positive four-condition mean-profile Pearson correlation.
Inputs are external arrays, protein matrices, supplied mappings and gene names.
Outputs include all-result audits, ranked and same-tissue dual-category candidates,
three-dataset salt agreement, counts, plotted observations, a PDF report and a
four-format profile figure. This is a post hoc candidate screen, not a replacement
for the original joint global BH/BY inference.

**Functions:** `supporting_settings()` defines thresholds and ranking;
`supporting_cells()` computes four observed-value condition means;
`supporting_profile()` reports descriptive Pearson/Spearman correlations;
`supporting_flag()` applies assay-wise significance, sign and profile criteria;
`supporting_order()` ranks without changing q values; `supporting_dual()` requires
salt and genotype support in the same tissue; `supporting_three_datasets()` adds
same-direction Array 1 salt support. `supporting_results()` orchestrates matching,
names, sensitivities and caching; `supporting_signal()` exports plotted pools.
`supporting_label()`, `supporting_page()`, `supporting_text()` and
`supporting_print_table()` format the five-page `supporting_report()`.
`supporting_figure()` shows raw pools and four-mean scatterplots.
`build_supporting_audit()`, `build_supporting_ranked()`, `build_supporting_dual()`,
`build_supporting_three()`, `build_supporting_summary()`, `build_supporting_signal()`,
`build_supporting_report()` and `build_supporting_profiles()` are independent
artifact wrappers; `supporting_registry()` registers all eight outputs. Local
`key()` functions build deterministic matching keys and `entry()` builds registry
records. See [supporting_trends.md](supporting_trends.md) for inference and naming.

### [`scripts/reproduce_supporting.R`](../scripts/reproduce_supporting.R)

Rebuilds all eight supporting-trend outputs from external inputs. Uses bootstrap,
`inside()` containment checks, an input cache and `run_artifact()` provenance.
Every output can also be reproduced individually with `scripts/reproduce.R`.

### [`tests/supporting_tests.R`](../tests/supporting_tests.R)

Uses synthetic values only to check dual-assay significance, effect direction,
positive profiles, missing and boundary cases, four-condition means/correlation,
same-tissue salt/genotype overlap, deterministic ordering and preserved q values.
No named helper functions are defined. No experimental measurements are embedded.

## Microarray preprocessing additions

### [`R/normalization.R`](../R/normalization.R)

Assesses scaling, quantile and cyclic-loess methods on processed Agilent signals;
writes a separate quantile input root without changing original measurements.
Inputs are two processed-signal workbooks, 32 sample exports, corrected mappings
and the three supplied normalized matrices. Outputs include feature-QC records,
method diagnostics, distribution figures and three regenerated RNA matrices.
Controls are held out of fitting. Protein and qPCR files are linked unchanged.

**Functions:** `normalization_settings()` defines the method contract;
`normalization_paths()` enumerates source files; `normalization_feature_file()`
reads complete feature exports; `normalization_sources()` validates workbook and
export identity/signals and collects QC; `normalization_duplicate_mad()` measures
repeat-spot dispersion; `normalization_variants()` fits all methods and metrics;
`normalization_comparison()` audits supplied scaling and new values;
`normalization_diagnostics_figure()` exports distributions;
`normalization_prepare_root()` creates the separate analysis input root.
Local helpers are `canonical()`, `control_rms()` and `inside()`.

### [`scripts/renormalize_microarrays.R`](../scripts/renormalize_microarrays.R)

Command-line preprocessing entry point. Requires distinct external original,
normalized-input and diagnostics roots. Writes preprocessing provenance, verifies
input/code hashes and marks completion. It has no positional arguments. The local
`inside()` helper checks containment. See [normalization methods](microarray_normalization.md).

### [`tests/normalization_tests.R`](../tests/normalization_tests.R)

Synthetic independent rank-mean verification, control-exclusion and separate-tissue
fitting checks, plus tie handling. No experimental measurements are included.
Top-level test commands; no named functions.

### [`scripts/compare_normalization_runs.R`](../scripts/compare_normalization_runs.R)

Compares completed scaling and quantile analysis releases. Reads genotype contrast
tables, dose joint tests and supporting rankings. Writes significant-set overlap,
lost/gained counts, fold-change agreement and candidate-level changes, without
fitting a new model or using discovery counts to select normalization. Distinct
external roots and completed input artifacts are required. Local helpers are
`inside()` for containment, `read()` for completed-file imports, `compare()` for
matched-universe sensitivity summaries, and `key()` for candidate identity.
Provenance records inputs, code, outputs and the Git revision.

## Computational extensions added 4 November 2026

The eight additional R files and the functions of each are described in the [extension script/function guide](computational_extensions.md#scriptfunction-responsibilities). The new entry point is `scripts/reproduce_extensions.R`, with 30 individually reproducible output IDs. Public reference downloads are separate from the experimental input root. These additions do not replace or relabel historical analyses.


### [`scripts/reproduce_extension_manuscript.R`](../scripts/reproduce_extension_manuscript.R)

Formats Supplementary Tables S18-S20 from the complete registered external
signature, TF-target and marker CSVs. It retains all eligible basal salt-signature
tests, and all primary global-BH findings in the stated TF and basal marker
families, with their matched sensitivities. It never refits a model or adjusts
a selected subset again. `SALT_EXTENSION_ROOT` selects the external source
folder; `SALT_OUTPUT_ROOT` selects a separate output folder. Call with `all`,
`--list`, or one table ID. Local `inside()`, `hash()` and `key()` helpers check
folder separation, record SHA-256 and match complete dataset/contrast/TF-or-marker
identities. Outputs are three CSVs with original statistics, provenance and
completion records. Supplementary Table S17 directly uses the existing
`dose_discoveries` CSV.


### [`R/bmc_robustness.R`](../R/bmc_robustness.R)

`bmc_dose_metadata()` parses doses and harvested-pair labels.
`bmc_collapse_logspots()` averages technical spots within transcript IDs and
unique transcript signals within loci. `bmc_fit_dose()` fits the paired eight-cell
model and complete fifteen-contrast/three-joint-test families.
`bmc_dose_omissions()` checks normalization reconstruction and runs all twelve
paired omissions under frozen/refitted normalization; its local `collect()`
records every focal estimate, interval, joint result and mean-shape count.
`bmc_omission_figure()` renders four focal-gene sensitivity panels with all
omissions. No function selects an omission to replace the primary data.

### [`scripts/reproduce_bmc_robustness.R`](../scripts/reproduce_bmc_robustness.R)

Validates separate roots, records input/code/output hashes with `inside()` and
`hash()` helpers, verifies or computes a fingerprinted external cache, and
reproduces six tables and one figure individually or together. See
`docs/bmc_revision.md` for the exact normalization and testing-family scope.

### [`tests/bmc_robustness_tests.R`](../tests/bmc_robustness_tests.R)

Uses synthetic paired pools to check complete-pair omission, independently
calculated cell-mean contrasts, interval construction and technical spot collapse.
The local `mu()` helper computes means separately from the model implementation.


### [`R/bmc_probe_mapping.R`](../R/bmc_probe_mapping.R)

`bmc_probe_evidence()` reconstructs complete consistent query sequences from the
original alignment workbook and audits valid/invalid rows and reported near
matches. `bmc_exact_cds_map()` maps full sequences on both strands to all pinned
CDS, rejecting cross-sequence boundary matches. `bmc_probe_audit()` joins
physical features, corrected IDs, reconstructed queries and exact/near targets.
`bmc_probe_sensitivity()` refits RNA contrasts and the complete eligible external
signature family using the declared strict probe subset. Local `num()` parses
the source number convention. No sequence gap is filled with guessed bases.

### [`scripts/fetch_probe_reference.R`](../scripts/fetch_probe_reference.R)

Downloads or verifies the exact Ensembl Plants release 63 TAIR10 CDS FASTA
against its pinned SHA-256 and writes an external reference manifest.

### [`scripts/reproduce_probe_mapping.R`](../scripts/reproduce_probe_mapping.R)

Validates external roots with `inside()`, records all input/code/output hashes
with `hash()`, runs reconstruction, independent mapping and strict sensitivity,
and exports thirteen tables plus the external recovered-probe FASTA.

### [`tests/probe_mapping_tests.R`](../tests/probe_mapping_tests.R)

Checks independently constructed exact, reverse-strand and duplicate-gene
examples and rejects an apparent match crossing concatenated CDS boundaries.


### [`R/bmc_tf_correlation.R`](../R/bmc_tf_correlation.R)

`bmc_residual_coordinates()` whitens the harvested-pair covariance and obtains
normalized residual coordinates in the residual subspace of the original design.
`bmc_mean_residual_correlation()` calculates the mean off-diagonal correlation
without constructing a gene-by-gene matrix. `bmc_tf_correlation()` applies these
estimates to the fixed DAP-seq target sets and the existing moderated statistics,
then adjusts the full test family by BH and BY. Negative correlation estimates
are set to zero for the test; raw estimates are retained in the audit.

### [`scripts/reproduce_tf_residual_correlation.R`](../scripts/reproduce_tf_residual_correlation.R)

Refits the primary gene-level assay objects, reads the pinned DAP-seq resource,
and exports the correlation audit and complete enrichment sensitivity table.
`inside()` enforces separate input/output roots and `hash()` records provenance.

### [`tests/tf_residual_correlation_tests.R`](../tests/tf_residual_correlation_tests.R)

Checks the correlation calculation against limma on synthetic unpaired and
paired designs, and verifies that the covariance whitening produces identity.


### [`scripts/reproduce_bmc_manuscript_tables.R`](../scripts/reproduce_bmc_manuscript_tables.R)

Formats five reviewer-response display tables from complete registered results:
probe reconstruction and focal mapping, flagged-pair influence under both
normalization strategies, both exploratory candidate screens, every eligible
strict-probe basal external-signature test, and every all-setting TF comparison
for the union of global-BH findings under any setting. It also writes the
comparison figure with a separate caption. `inside()` separates
roots and `hash()` records inputs/code/outputs. No model, source P value or
source-family correction is recomputed after display selection. Run one ID,
`all` or `--list`; the root variables are documented in `docs/bmc_revision.md`.


### [`R/bmc_score_sensitivity.R`](../R/bmc_score_sensitivity.R)

`bmc_score_design()` reconstructs the original categorical pool design and
contrasts from the registered sample table. `bmc_unmoderated_scores()` refits
ordinary residual-variance score models with the original paired-score covariance,
checks equality of point estimates, and adjusts every eligible test by BH/BY.
`bmc_score_inference_figure()` displays every focal-reference Array 2 pool and both basal intervals, with explanations in its caption. The module compares uncertainty with empirical-Bayes moderation without changing scores,
reference selection or target eligibility. Negative outcomes remain in the export.

### [`scripts/reproduce_score_inference_sensitivity.R`](../scripts/reproduce_score_inference_sensitivity.R)

Reads completed sample, external-score and marker-score results, exports all 141
external and 781-per-mode marker tests plus every basal external display row and the external-score inference figure in four formats. Each output ID can be reproduced separately.
`inside()` separates roots and `hash()` records source, code and output hashes.

### [`tests/score_sensitivity_tests.R`](../tests/score_sensitivity_tests.R)

Verifies score contrasts, standard errors, residual degrees of freedom and
ordinary t P values against an independent synthetic `lm` calculation.

The TF comparison helpers in `R/bmc_tf_correlation.R` also support the current display revision. `bmc_tf_comparison()` selects the union of globally adjusted BH findings across fixed 0.01, fixed 0.05 and estimated residual-correlation families and copies every original BH/BY value and direction for all three settings. `bmc_tf_comparison_figure()` plots those BH values alongside estimated correlations, writing separate explanatory captions. `scripts/reproduce_bmc_manuscript_tables.R` exposes this as `figure_TF_correlation_sensitivity`; its probe audit now distinguishes physical spots, query identifiers and unique reconstructed sequences. The existing TF test entry point checks selection against deliberately discordant settings without recomputing q values.

`bmc_target_overlap()` computes all pairwise intersections, unions, Jaccard indices and fractions of the smaller set from unique target memberships. The existing TF test checks these against explicit overlapping and disjoint toy sets. The display entry point requires the pinned AGRIS reference and the complete Array 1 universe; it exports `table_TF_overlap` and `table_TF_target_members` and checks target counts against original CAMERA output. Panel C of `bmc_tf_comparison_figure()` plots the same measured-universe overlaps. These are descriptive set calculations, without enrichment P values or additional biological replication.
