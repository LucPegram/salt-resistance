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
