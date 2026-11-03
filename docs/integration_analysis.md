# Transcript and protein integration plan

The plan is specified before inspecting integrated results. These experiments
use independent plant pools. No sample-level transcript/protein pairing, tissue
mixture reconstruction, or joint latent-factor model is justified.

## Questions and comparisons

1. Compare basal genotype effects in the No-0-background arrays and each tissue
   proteome. This has the closest condition alignment, although the arrays combine
   tissues and the protein fractions are plasma-membrane enriched.
2. Compare No-0 salt, NIT2 salt, saline genotype and genotype-by-salt effects.
   Array salt is 100 mM and protein salt is 75 mM; these are related contrasts,
   not identical treatment effects or a test of translational regulation.
3. Separately compare the 75-versus-0 mM Col-0 array effect within each tissue
   with the 75-versus-0 mM No-0 protein effect. This matches dose and tissue but
   changes background and does not establish matched growth/harvest protocols.

The complete family comprises twelve tissue/comparison combinations: five
genotype-array contrasts in two proteomes plus the two contextual Col-0 pairs.

## Mapping and models

For integration only, average log-scale measurements of unique exact transcript
IDs within each AGI locus, with equal weight per transcript ID; repeated probe
rows have already been averaged. Refit locus-level limma models using the same
robust/trend settings. These locus tests do not replace the manuscript's original
transcript-level tests. The Col-0 context fits four dose means separately in each
tissue and tests 75 minus 0; no cross-treatment matching is assumed.

Match complete protein-group identities to supplied AGI annotations. Retain
single-accession groups with reported proteotypic support, a valid unique AGI,
and at least two measured pools in every condition. Exclude duplicate locus
assignments rather than select the most significant group. Audit every exclusion,
including absent array coverage and the candidate proteins.

Report all eligible genes, per-assay effects, pointwise intervals, P values and
within-assay BH values. Describe Pearson/Spearman effect correlations and sign
agreement without treating genes as independent biological replicates. Do not
assign a correlation P value or a gene-bootstrap confidence interval.
After the first integrated summary was inspected, add a diagnostic sensitivity
omitting NIT2 from correlations to assess transgene influence. This diagnostic
does not change the joint-test universe and is not a new primary hypothesis.

## Evidence in both assays

Use a directional intersection-union test: convert each two-sided P value into
left/right one-sided values using its effect sign; take the maximum within each
direction; multiply the smaller directional maximum by two, capped at one.
This requires both assays to support the same direction and prevents a single
strong assay from dominating a combined P value. It makes no claim of a common
effect size or of a replicated physiological mechanism.

Report BH within each of twelve families and BH/BY across all gene/family tests.
BY is a dependence-robust sensitivity. The primary joint discovery criterion is
global BH below 0.05. Retain nonsignificant effects and intervals. Separate
opposite point estimates from evidence of opposite changes in both assays.

Sensitivities use the existing complete-case protein model and, separately,
restrict primary comparisons to loci with one mapped transcript ID. Correct each
sensitivity across its own complete set of gene/family tests; do not select the
most favourable sensitivity as the primary result.

## Focused pathway integration

Use nine GO sets chosen from the project's salt, transport, auxin and redox
questions: GO:0009651, GO:0006814, GO:0006813, GO:0015992, GO:0009733,
GO:0000302, GO:0006979, GO:0042744 and GO:0016209. During final table QA,
GO:0015992 was identified as a secondary ID of GO:1902600 in pinned GO.db.
The planned proton-transport concept is retained under GO:1902600. Resolve
all supplied secondary GO IDs using the pinned GOTERM Secondary slots before
ancestor propagation, and de-duplicate gene/term assignments. This correction
was triggered by identifier quality, not statistical significance. Recover term names from pinned
GO.db, not guessed labels. The supplied association release remains unknown.

Use the same uniquely mapped, fully observed gene universe for both assays in
each comparison. Test sets with at least five measured members and no more than
80 percent of that universe; preserve untested sets and their exclusion reasons.
Use CAMERA with inter-gene correlation estimated separately for each set and
negative correlation disallowed, plus a fixed-correlation 0.01 sensitivity.
Apply the same directional conjunction and global BH/BY correction across all
testable set/family combinations separately for each correlation setting.
Supply exact gene universes and memberships. CAMERA tests relative enrichment,
not pathway activation or ROS/ion flux, and does not prove every member agrees.

## References

- Wu D, Smyth GK. 2012. Camera: a competitive gene set test accounting for
  inter-gene correlation. Nucleic Acids Research 40:e133.
  https://doi.org/10.1093/nar/gks461
- Djordjilovic V, Sofer T, Dreyfuss JM. 2026. Directional replicability: When can
  the factor of two be omitted. Statistics and Probability Letters 233:110662.
  https://doi.org/10.1016/j.spl.2026.110662

Here the conservative factor of two is retained. The scientific terminology is
joint molecular association because the assays do not sample identical material
or treatment regimes. All data and generated results remain outside GitHub.

## Additional diagnostics after visual review

The response to visual review adds four reproducible tables without changing
primary tests or eligibility: selected transcript tissue-by-dose tests from the
existing paired model; descriptive influence of omitting the three loci selected
by primary joint BH, plus every leave-one-locus Pearson correlation; complete-case
and single-transcript correction values beside each primary joint hit; and
per-cell processed-measurement ranges and reported identity support for all
highlighted loci. The omission set is explicitly post hoc and is not a new
discovery analysis. No samples or genes are excluded from the primary figures
or correlations. Gene counts are not biological replication. The measurement
audit independently checks that cell means reproduce highlighted model contrasts;
it does not validate raw processing, probe specificity or peptide spectra.
