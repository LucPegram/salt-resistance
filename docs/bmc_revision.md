# Computational response to the Plant Physiology review

The full-data quantile analysis remains primary. No additional wet experiments
are introduced, and no biological sample is removed from that analysis.

## Harvested pair influence

`scripts/reproduce_bmc_robustness.R` performs all twelve root-shoot harvested-pair
omissions in Array 1 under both frozen and refitted quantile normalization.
Refitting uses the same ordinary physical spots and normalizes each organ
separately before exact transcript collapse and locus averaging. The full
reconstruction must agree with the primary locus matrix within 1e-10 log2 units.
Each omission retains 22 profiles, with at least two harvested pools per dose.
No cross-treatment batch pairing is inferred.

Each fit uses the same 27,089 loci and fifteen planned contrasts. Global BH is
within all locus-contrast tests of each fit; the three joint tests have a separate
complete global family. The fit-specific families are sensitivity analyses, not
independent replications or opportunities to select a favorable P value. Effect
ranges include the primary fit. Sign agreement among effects of at least 0.5
absolute log2 units is a prespecified descriptive stability summary. Discovery
counts and mean-shape counts are provided but threshold retention is not treated
as proof of reliability. The flagged root belongs to pair `75_2`; every other
pair is omitted on the same basis.

Set `SALT_R_LIB`, `SALT_DATA_ROOT` (primary quantile inputs), `SALT_RAW_ROOT`
(original extracted inputs), `SALT_OUTPUT_ROOT`, and `SALT_QC_CACHE` to separate
external locations. Run `Rscript --vanilla scripts/reproduce_bmc_robustness.R all`
or replace `all` with one ID returned by `--list`. Six CSV outputs cover genome
metrics, all nine focal genes and their intervals/tests, joint-test counts,
mean-shape counts, and per-gene/contrast ranges across every omission. The figure
shows all twelve omissions for four focal 75-versus-0 mM root effects. Its
explanatory text is in a caption companion. Input, code and output hashes are
recorded; unchanged derived results can be retrieved from a verified cache.

The computational check does not reconstruct scanner images, establish probe
specificity, or resolve the experimental unit of salt application.

## Probe reconstruction and mapping sensitivity

The supplied alignment workbook contains full query sequence evidence. Complete,
consistent A/C/G/T query sequences are recovered only when alignment coordinates
span the entire stated query length. Gaps introduced by alignment are removed;
reverse query orientation is handled explicitly. Partial or conflicting queries
are not filled or silently promoted to complete sequences.

`fetch_probe_reference.R` pins the Ensembl Plants release 63 TAIR10 CDS FASTA by
URL and SHA-256. `reproduce_probe_mapping.R` independently searches every
recovered complete query, on both strands, against every reference CDS using
Biostrings. Matches crossing concatenated CDS boundaries are rejected. Isoforms
are collapsed to gene loci only for the locus-uniqueness diagnostic.

The strict sensitivity uses probes with one exact CDS locus matching the supplied
corrected mapping and no alternative locus among supplied alignments with at
least 95% identity and 95% query coverage. These are stated annotation criteria,
not a physical cross-hybridization model. Exact-match uniqueness does not prove
hybridization specificity or isoform resolution, and CDS-only searching does not
exclude off-target binding elsewhere in the transcriptome/genome. All failures,
missing sequences, alternative loci and original corrected identifiers remain in
the complete feature audit. The recovered FASTA is an external output, never
committed as project data.

After fitting normalization on the same complete ordinary-spot universe, the
sensitivity restricts downstream RNA models to these probes. It refits all dose
and genotype contrasts and the complete eligible independent salt-signature
family. Different eligible loci and score universes remain explicit; the
sensitivity does not replace the declared primary results or compare discovery
counts as if they were biological effect sizes.

Set `SALT_DATA_ROOT`, `SALT_RAW_ROOT`, `SALT_REFERENCE_ROOT` (the existing eight
extension resources), `SALT_PROBE_REFERENCE_ROOT` and `SALT_OUTPUT_ROOT` to
separate external roots. Run `Rscript --vanilla scripts/fetch_probe_reference.R`,
then `Rscript --vanilla scripts/reproduce_probe_mapping.R`. Thirteen CSV outputs,
the complete recovered probe FASTA and per-output provenance are generated.


## TF target residual correlation sensitivity

`reproduce_tf_residual_correlation.R` preserves the original moderated feature
statistics, target memberships and eligible test family. It estimates average
target residual correlation from the original assay design. The paired dose
assay first whitens expression and design with the fitted harvested-pair
covariance; all sets use residual coordinates, not correlation across condition
means. Negative estimates are truncated to zero for rank CAMERA; both raw and
used estimates are exported. Global BH and BY cover all 5,935 tests separately
from the fixed 0.01 and 0.05 analyses. Residual degrees of freedom are 16 for
Array 1, four for Array 2 and eight for each protein assay, so estimated
correlations remain uncertain. This sensitivity does not establish TF activity.

Set `SALT_R_LIB`, `SALT_DATA_ROOT`, `SALT_REFERENCE_ROOT` and `SALT_OUTPUT_ROOT`
to external locations and run
`Rscript --vanilla scripts/reproduce_tf_residual_correlation.R`. The script
exports `TF_residual_correlation.csv` and `TF_estimated_correlation_tests.csv`
with complete input/code/output hashes and package versions.
