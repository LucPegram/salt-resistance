# Dose-response, external benchmarking, TF targets and tissue/development signatures

These supplementary analyses extend the primary quantile-normalized release without replacing its transcript-ID tables, available-observation protein models, or registered manuscript figures. No current-project measurements, public reference data, generated results, or manuscripts are included in this repository. The user explicitly excluded the suggested 10x dataset; it is not used in this workflow.

## Run and reproduce a single output

Use the R environment in `renv.lock`. Set these environment variables to directories **outside this repository**:

```sh
export SALT_DATA_ROOT=/absolute/path/to/analysis_inputs/quantile
export SALT_REFERENCE_ROOT=/absolute/path/to/external_references/computational_extension
export SALT_OUTPUT_ROOT=/absolute/path/to/analysis_outputs/computational_extensions
export SALT_MODEL_CACHE=/absolute/path/to/existing/revision_model_cache
# Optional: existing library restored from renv.lock
export SALT_R_LIB=/absolute/path/to/R-library

Rscript --vanilla scripts/fetch_extension_references.R
Rscript --vanilla tests/extension_tests.R
Rscript --vanilla scripts/reproduce_extensions.R --list
Rscript --vanilla scripts/reproduce_extensions.R all
# Examples: reproduce one specific table or figure
Rscript --vanilla scripts/reproduce_extensions.R dose_contrasts
Rscript --vanilla scripts/reproduce_extensions.R figure_external_signatures
Rscript --vanilla scripts/reproduce_extensions.R marker_tests
```

Each of the 30 output IDs has its own directory, original/reference/code/output SHA-256 hashes, embedded microarray-normalization record, statistical settings, package/session information, and completion marker. Figures have PDF, SVG, PNG and 600-dpi TIFF versions. A local derived-data cache is keyed by code, settings, R/package versions and all inputs and verified by a content hash. It is never uploaded. A single-output command computes or verifies the shared cache but writes only the selected output.

Reference retrieval records URLs, byte sizes, timestamps, and SHA-256 hashes in an external `download_manifest.csv`. Archive this manifest and the exact eight downloaded resources with the analysis release. Future versions of public resources are not assumed equivalent: compare their hashes with the release provenance. GEO annotations are interpreted from the retrieved version, not silently updated to current gene mappings.

## Model and testing families

### Root-shoot dose patterns

Use all mapped Array 1 loci, averaging unique transcript log2 signals within each AGI locus. This is a gene-level counterpart of the established transcript-ID analysis, so discovery counts should not be substituted for the original transcript-ID counts. Categorical tissue-by-dose means preserve the four concentrations. Root and shoot from each harvested pool share a block; the 12 pools are defined by dose and replicate, with no cross-treatment pairing. Limma estimates the paired-pool consensus correlation and uses robust, intensity-trend moderation.

There are five comparisons per tissue: 50-0, 75-0, 100-0, 75-50 and 100-75 mM; five direct tissue-difference comparisons are **shoot response minus root response** for those same increments. BH is calculated across all loci and all 15 contrasts (`q_global`), with per-contrast BH also retained (`q`). Three joint three-degree-of-freedom tests (root, shoot, tissue response difference) form a separate global BH family.

Increasing/decreasing/non-monotonic labels classify the four observed group means. They are descriptive, not statistically established monotonic programs. The first globally significant **tested** concentration is reported without interpreting it as a response threshold. No EC50 or elaborate nonlinear curve is estimated. Quantile normalization was performed separately by organ, so these results address relative expression responses, not absolute RNA abundance or total transcription differences between organs. A difference in significance between two organs is not a significant organ difference; use the direct Difference tests.

### Independent salt-response references

References were selected by species, experimental design and available processed measurements before inspecting current-project concordance:

- [GSE5623](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5623) with its matched-time [GSE5620 controls](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5620): AtGenExpress Col-0, 16-day plants, roots and shoots, 150 mM NaCl, 0.5, 1, 3, 6, 12 and 24 h. The actual sample titles identify **two** biological arrays per tissue/time/treatment cell, despite generic series text discussing triplicate designs. Deposited positive MAS5 signals are log2 transformed and quantile normalized jointly across control/salt times separately within tissue, before unambiguous probe-to-AGI collapse using GPL198. No current-project values enter this normalization.
- [GSE223207](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223207): independent hydroponic Col-0 roots, 16-day plants, 25/50/75/100 mM, 15 min, two arrays per condition. The deposited GCRMA log2 matrix is used without a second logarithm or additional normalization. All 16 deposited solute/control groups are retained as separate categorical means for residual variance estimation; only the four NaCl-control contrasts define salt signatures. GPL17133 assignments containing multiple distinct AGI loci are excluded, rather than choosing a first hit.

External models use categorical means and limma robust/trend moderation. No undocumented replicate-number block is inferred. Public differential-expression tables retain their own within-contrast BH values. A public Up or Down signature requires BH q < 0.05 and |log2FC| >= log2(1.5), selected entirely in the external experiment. Both halves must cover at least ten measured genes in a target assay.

For each biological target sample, the signed salt score is mean fractional expression rank of external Up genes minus the corresponding mean of external Down genes. Ranks are computed within each sample over that assay's fixed eligible universe. Scores are fitted with the appropriate project design; the paired correlation is re-estimated from scores for Array 1. External root signatures are used in separated root samples and shoot signatures in separated shoot samples; both are eligible in combined-tissue Array 2. All eligible signature-by-selected-contrast tests across assays form one global BH family, with BY sensitivity (`q_BY`). Correlations of gene-level effect vectors and sign agreement are descriptive, without treating genes as independent biological replicates.

These references represent early stress and differ from the two-week experiment. They test resemblance or persistence of response direction, not equivalent treatment, stress memory, causation, or protection. GSE5620 is the control arm of the same AtGenExpress experiment, not a second independent salt study. Cackett et al. 2022 is not used as independent validation.

### TF target enrichment

[AGRIS AtRegNet](https://agris-knowledgebase.org/Downloads/AtRegNet.zip) provides an accessible reference when PlantRegMap downloads are unavailable. Use only rows with valid TF/target AGI identifiers, explicit DAP annotation, direct target type, and PMID 27203113 ([O'Malley et al. 2016](https://doi.org/10.1016/j.cell.2016.04.038)). Historical CSV rows with extra fields are excluded and counted; no column-shifting repair or evidence-type substitution occurs. DAP rows marked "Unconfirmed" by AGRIS are not claimed to be confirmed in-vivo regulation.

Use at least ten measured targets and exclude sets larger than half the assay universe. Rank-based `limma::cameraPR` tests the paired-model moderated t statistics for Array 1 and the relevant factorial-model statistics for Array 2/proteins. It is a competitive test, accounting for an assumed positive inter-gene correlation. Primary correlation 0.01 and sensitivity 0.05 are analyzed separately; each has one global BH/BY family across all TFs, assays and selected contrasts. Per-assay/contrast BH is also retained. Do not select a preferred correlation after seeing significance. The published target-to-gene assignment is accepted as a reference annotation; binding distances were not rederived from genomic peaks.

Target enrichment is a candidate-regulator hypothesis, not TF activity, binding in the experimental tissue, or activation/repression by that TF. Closely related TFs may have highly overlapping targets. The focal-target table is a reference annotation table, not evidence for new interactions.

### Development and tissue-marker signatures

Use independent [Lee et al. 2025 atlas](https://doi.org/10.1038/s41477-025-02072-z) Tables S1-S2. Table S1 includes published root cell/developmental markers from Shahan et al. 2022, literature mesophyll/vascular/epidermal/SAM markers, and dividing/G1-S markers. Preserve the published membership and original reference. Pool stage-prefixed root markers into Meristem, Elongation and Maturation sets. For each annotated 21-day rosette cluster, select at most 50 markers with published adjusted P < 0.05, log2FC >= 0.5 and pct.1 > pct.2, sorting by decreasing log2FC and AGI; union clusters with the same published cell-type label. Exclude Unknown clusters. These rules do not use the salt project's outcomes.

Require at least five measured genes per set. Mean within-sample fractional rank is a marker-expression score, **not a cell proportion**. Test root-specific signatures in root samples, shoot-specific signatures in shoot samples, proliferation signatures in both, and all organ signatures in combined-tissue Array 2. All anatomically eligible marker-by-selected-contrast tests across assays form one global BH/BY family. A separate sensitivity excludes the prespecified transport candidates, NIT2 and the three previously selected cross-assay loci; its BH/BY family is independent of the primary family. Sets falling below minimum coverage remain in the coverage audit but are not tested. Marker scores are not used as correction covariates in the primary expression models.

All signature/TF protein analyses use only uniquely mapped, single-accession, proteotypically supported, fully observed protein loci. Fixed membership avoids sample-dependent scores from missing values and introduces no imputation. These are supplementary complete-observation analyses, not replacements for the primary protein contrasts that admit two observed pools per cell. Protein representation of a marker set can be sparse; inspect coverage. Combined RNA cannot localize a genotype effect to an organ.

## Output guide and interpretation

- `dose_*`: complete contrasts with confidence intervals and gene names, joint tests, mean shapes, means and discovery summaries. Figures show all mean-shape categories, the top 20 joint loci per tissue, and adjacent-dose estimates for focal candidates. Selection and display ranges are disclosed.
- `external_*`: sample metadata, probe-mapping audit, all external salt contrasts, frozen signature membership, coverage, per-sample scores, score tests and descriptive effect concordance. Signatures failing coverage are retained in the audit.
- `TF_*`: provenance filtering audit, all TF coverage, complete primary/sensitivity enrichment tests, and original DAP annotations linking TFs to focal genes. The figure ranks 20 TFs by their minimum primary global q; it does not imply all 20 are significant.
- `marker_*`: external memberships and source references, all coverage, per-sample scores, primary and candidate-exclusion score tests. The figure includes all eligible stage/cycle sets and 20 strongest other sets, selected by minimum primary global q.
- `assay_universes`, `project_samples`: the exact eligible dimensions and biological sample identities used by these extensions.
- `analysis_report`: automatically generated methods/results/limitations and primary-source references. It is a supplementary report, not a replacement manuscript.

Signed evidence heatmaps use direction times -log10(global BH q), capped at six for display, with stars at q < 0.05 and grey for untested comparisons. Scores, effects and uncapped P/q values remain in tables. TF directions describe the target set's expression shift, not the regulator's activity. Every displayed finding should be checked against its sensitivity and measurement coverage before being included in mechanistic claims.

## Script/function responsibilities

| File | Purpose |
| --- | --- |
| `config/extensions.R` | Reference URLs and prespecified coverage, effect, selection and sensitivity settings. No gene-expression values or derived gene lists. |
| `R/extension_inputs.R` | Validates GEO matrices/platforms; preserves ambiguous-mapping audits; reads the DAP reference safely to EOF; extracts independently defined atlas memberships. |
| `R/extension_models.R` | Fits paired categorical dose contrasts, independent public salt references, sample-level salt/marker scores and competitive TF target enrichment. |
| `R/extension_figures.R` | Generates the six publication-style multi-format figures from computed results with explicit selection and legends. |
| `R/extension_release.R` | Builds the 23 table objects, five-page report and 30-output registry; all reported counts are calculated from results. |
| `scripts/fetch_extension_references.R` | Downloads/validates eight external resources or verifies existing copies; records URL and hash provenance. |
| `scripts/reproduce_extensions.R` | Runs all extensions or one output ID; verifies cache and unchanged inputs/code; writes individual artifact provenance. |
| `tests/extension_tests.R` | Tests contrast signs/pairing, ambiguous mapping, fixed-membership ranks, untested heatmap cells, duplicate rejection and DAP streaming/field validation using synthetic fixtures. |

Additional primary-method reference: [limma User's Guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf). No machine-learning classifier, inferred experimental batch, cell-fraction deconvolution, or cross-assay sample pairing is introduced.
