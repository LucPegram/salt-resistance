# Reproducible analysis code

This repository contains R source code, an environment lock, historical source
scripts and reproduction instructions. Experimental data, source archives,
manuscripts, statistical results, generated figures and local run records are
excluded. Every registered output can be rebuilt separately from external inputs.

## Organization

- `R/inputs.R`: validated imports and transcript/protein filtering.
- `R/revision_models.R`: factorial contrasts, paired dose models, missingness
  sensitivities and disjoint-control resemblance calculations.
- `R/revision_enrichment.R`: directional enrichment with eligible universes.
- `R/revision_figures.R` and `R/revision_manuscript_tables.R`: manuscript outputs.
- `R/integration_*.R`: locus-matched transcript/protein comparisons, joint
  directional evidence, focused CAMERA pathways and integration figures.
- `R/reviewer_diagnostics.R`: candidate dose tests, locus influence, joint-hit
  sensitivities and processed-measurement audits requested during review.
- `R/toprank.R`: three dataset-specific top-gene/protein rankings, naming audits
  and treatment-response figures; methods in `docs/top_rankings.md`.
- `R/manuscript_evidence.R`: historical summaries and numerical audits.
- `config/`: external input paths, analysis settings and output registry.
- `scripts/`: reproduction, dependency checks and code-publication checks.
- `tests/`: synthetic examples and independent statistical calculations.
- `legacy/`: original source scripts, inactive and retained for inspection.
- `docs/`: analysis decisions, output commands and limitations.

The revised workflow fits models from supplied quantitative matrices. It does not
claim to reconstruct undocumented upstream acquisition or normalization. Archived
scripts are not executed automatically or represented as validated analyses.

## Environment

Use R 4.5.2 and the package versions in `renv.lock`, including Bioconductor 3.22.
After installing `renv`, run `renv::restore()` from this repository, then check:

```sh
Rscript --vanilla scripts/check_environment.R
Rscript --vanilla tests/run_tests.R
Rscript --vanilla tests/integration_tests.R
Rscript --vanilla tests/toprank_tests.R
```

The lock records exact installed sources and versions. These checks do not prove
a clean-machine restoration. Archived scripts may reference packages outside the
active workflow. Packages are never installed silently during analysis.
If using a separate local R library, set `SALT_R_LIB` to its absolute path.

## External inputs and outputs

Keep inputs outside this repository. The input root contains the extracted
`Microarray_for_Stephano` and `Data_from_Luc_2` directories and the separately
obtained `references/P19456.fasta`. Relative paths are in `config/inputs.R`.
The analysis does not upload measurements or download substitute experiments.

```sh
export SALT_DATA_ROOT=/absolute/path/to/extracted/source_data
export SALT_OUTPUT_ROOT=/absolute/path/to/local/analysis_outputs
Rscript --vanilla scripts/reproduce.R --list
Rscript --vanilla scripts/reproduce.R table_revision_primary_candidates
Rscript --vanilla scripts/reproduce.R figure_revision_protein
Rscript --vanilla scripts/reproduce_figures.R
Rscript --vanilla scripts/reproduce_integration.R
Rscript --vanilla scripts/reproduce_toprank.R
Rscript --vanilla scripts/reproduce.R all
```

The output root defaults to the ignored `outputs/` directory. Each output records
input/code/output SHA256 hashes, effective settings, package versions, a Git
revision when available, session information and a completion marker. A missing
completion marker means that the run did not finish successfully.

An optional `SALT_MODEL_CACHE` directory outside code and source inputs can reuse
model results. Keys include input, relevant source, setting and package-version
hashes. Cached payloads are hash-checked. A fresh run without a cache is supported.
Caches are derived data and must never be pushed.

The seven current manuscript figures use consistent publication styling and export
vector PDF/SVG, 600-dpi TIFF and PNG previews. All formats are hash-recorded.
`docs/figure_design.md` describes visual encodings and display qualifications.

## Reproduction and audit documentation

[`docs/R_script_guide.md`](docs/R_script_guide.md) explains the purpose, inputs,
outputs and functions of every current and archived R script.

`docs/manuscript_outputs.md` maps every main/supplementary table and figure to an
independent command. `config/artifacts.R` and `R/revision_registry.R` provide the
registry. Builders never require a pre-existing output directory.
`docs/analysis_decisions.md` explains model choices and unresolved constraints.
`docs/integration_analysis.md` specifies gene mapping, contrast alignment,
directional conjunction, multiplicity, pathway universes and sensitivities.

```sh
Rscript --vanilla scripts/inventory_sources.R
```

This inventories external files and historical analyses by path, hash and role,
without executing old scripts. `legacy/README.md` maps archived code to original
paths and hashes. The explicit maintenance command `scripts/import_legacy.R`
imports R files only and refuses to overwrite an archive.

Before committing or uploading, stage the intended source files and run:

```sh
Rscript --vanilla scripts/check_code_only.R
```

The allowlist rejects measurements, generated results, libraries, session state
and manuscripts. Public data deposition must be arranged separately.

The supporting-trend candidate screen is described in [docs/supporting_trends.md](docs/supporting_trends.md).
Run `Rscript --vanilla scripts/reproduce_supporting.R` for its eight independent
report, table, audit and profile-figure artifacts. Synthetic selection and
correlation checks: `Rscript --vanilla tests/supporting_tests.R`.
