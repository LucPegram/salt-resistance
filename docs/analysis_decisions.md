# Analysis decisions

## Units and input processing

One independently grown plant pool is a biological replicate. Plants within a
pool are not additional replicates. Root–shoot measurements from one harvested
pool are paired. Matching numbers do not establish cross-treatment pairing.

Primary microarray analyses begin with quantile-normalized matrices regenerated
from supplied processed signals; the original normalized matrices remain a
scaling baseline. See [normalization assessment](microarray_normalization.md). Remove unmapped rows,
log2-transform positive values without a pseudocount, and average repeated exact
transcript IDs. Retain transcript suffixes. Do not call transcript IDs gene loci
or resolved isoforms. Vendor background processing and probe remapping specificity are not reconstructed.

Retain complete protein-group identities. Historical retention requires more
than one sequence and at least two observed values in at least one condition.
Primary inference additionally requires at least two observations in every cell.
No missing observation is replaced by an artificial value in primary inference.

## Models and multiplicity

Four-cell genotype/salt models use seven prespecified contrasts: basal and saline
genotype, wild-type and overexpressor salt responses, interaction and averaged
genotype/salt effects. limma uses robust moderation with an abundance trend.
Report pointwise moderated intervals and BH correction across features within
each contrast.

`q_all_contrasts` is within-contrast BH multiplied by seven, capped at one. This
is conservative Bonferroni correction across seven contrast families and is
nested within per-contrast discoveries. It is deliberately distinct from pooling
all P values in a single BH family, whose discoveries need not be nested within
separate-contrast results. Never interchange these correction labels.

Dose response uses eight tissue-by-concentration means and a harvested-pool
correlation estimated by duplicateCorrelation. Separate three-degree-of-freedom
tests evaluate root, shoot and tissue-by-dose effects. No polynomial shape is
imposed. Do not pool ecotypes, platforms or unmatched biological samples.

## Protein sensitivity analyses

Sensitivities include complete cases; at least one observation per cell and eight
overall; sample-median centering of the primary eligible set; and proDA dropout
models on the broader retained set. proDA uses log2 data and default location and
variance moderation. Convergence and warnings are recorded. Its assumptions are
not established by model execution and absent conditions remain unobserved.

The historical pipeline filled fully missing conditions with a minimum, then used
transposed KNN. Its ordinary ANOVA and ratio estimates are audit material.
Complete-case calculation agreement does not validate imputation-dependent rows.
Counts across different eligibility sets and assumptions are not directly
comparable estimates of biological response magnitude.

## Identity and pre-stress resemblance

Inspect full accession groups, gene labels and proteotypic sequence support.
A single accession does not establish specificity when unique support is absent.
Candidates failing the primary observation rule remain explicitly ineligible.
Label secondary estimates by their actual model.

Compare basal genotype and wild-type salt contrasts within an experiment. Shared
controls induce correlated measurement error. Evaluate every nonempty
complementary untreated-control split to remove shared measurements between axes.
Use complete proteins to avoid split-dependent missing-value handling. Ranges are
descriptive sensitivities, not confidence intervals or significance tests treating
molecules as biological replicates. Protection requires physiological evidence.

## Enrichment

Use supplied gene–GO associations with ancestor propagation through pinned GO.db.
Their original release is undocumented. Eligible genes require a single accession,
at least one reported proteotypic sequence,
unambiguous unique AGI assignment, primary-model eligibility and annotation.
Exclude duplicate locus mappings rather than choosing the first group.

Separate positive and negative basal, saline-genotype and interaction queries at
within-contrast q below 0.05. Test terms with at least five eligible genes and at
most 80 percent of the universe, including zero-hit terms. Report BH within each
family and across all twelve tissue/contrast/direction families. Historical
agriGO backgrounds and thresholds remain distinct audit analyses.

## Unavailable records

Do not invent line identities, treatment protocols, container independence,
qPCR replicates, reference genes, localization, acquisition parameters, phenotypes,
fluxes or blots. qPCR means remain descriptive. Canonical sequence mapping does not
establish site localization, occupancy, differential phosphorylation or activity.
Reproduction from processed matrices is not full raw-assay reproduction.

## Computational extensions, 4 November 2026

Added gene-level paired categorical dose contrasts, two independent public salt experiments (GSE5623 with its GSE5620 control arm; GSE223207), original-paper DAP-seq target enrichment via AGRIS, and independent atlas developmental/tissue signatures. Full methods, multiple-testing families, complete-observation protein eligibility, coverage audits, candidate-exclusion and correlation sensitivities are in [computational_extensions.md](computational_extensions.md). The 10x dataset is excluded following the user's explicit instruction. Main-release outputs and original measurements are not overwritten.
