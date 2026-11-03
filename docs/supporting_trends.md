# Supporting transcript and protein trends

This additional post hoc candidate screen answers which genes show significant,
same-direction salt or genotype effects in both assays. Explaining every
cross-assay disagreement is not a prerequisite for selecting supported candidates.
The original integration models, eligible mappings and all correction families
are preserved; the manuscript's primary analyses are not replaced.

## Criteria and comparisons

Use the closest background comparison: Array 2 (No-0 and NIT2, combined tissue,
0/100 mM NaCl) with each tissue proteome (0/75 mM). Evaluate five contrasts:
No-0 salt, NIT2 salt, basal genotype, genotype under salt, and genotype by salt.
Require the original within-assay BH q < 0.05 in **both** assays and matching
nonzero effect signs. Keep all original P values, pointwise intervals, joint
family/global BH and global BY values, rather than recomputing significance
within a selected list.

For each gene/protein-tissue pair, calculate four arithmetic condition means
on the log2 scale, in the order No-0 control, No-0 salt, NIT2 control, NIT2 salt.
Use observed protein values; do not impute missing pools. All eligible protein
cells already have at least two observations. Calculate Pearson and Spearman
correlations between the array and protein four-mean profiles. The screen also
requires a positive Pearson correlation. No correlation P value is calculated:
four condition means are not four matched biological replicates, and selecting
on significance and correlation is descriptive prioritization.

The RNA signal is a locus summary with equal weights for unique transcript IDs,
after repeated exact transcript spots are collapsed on the log2 scale. Protein
assignments require original retention, a single accession, at least one reported
proteotypic sequence, valid unique AGI, no duplicate retained mapping, array
coverage and at least two observed pools in every cell. See
[integration_analysis.md](integration_analysis.md) for the unchanged models.

## Ranking and evidence levels

Order supported comparisons by the larger of the RNA and protein q values,
ascending, then profile Pearson r descending, smaller absolute log2 effect
descending, AGI, tissue and contrast. This balances support rather than allowing
one exceptionally strong assay to dominate. Display ten salt candidates per
protein tissue and ten distinct genes per genotype comparison in the report;
all qualifying results remain in CSV.

The salt-plus-genotype priority table requires at least one supported salt
contrast **and** one supported genotype contrast for the same gene and the
same protein tissue. Select the contrast pair with the smallest worst q across
its four assay tests, with deterministic contrast-name ties. A root salt hit
and shoot genotype hit alone do not qualify. These two significant comparisons
are not a significant genotype-by-salt interaction.

A separate three-dataset salt table requires the No-0 salt screen plus a
significant, same-direction Array 1 effect at 75 versus 0 mM in the corresponding
tissue. Array 1 uses Col-0; Array 2 combines tissues at 100 mM. The same protein
contrast enters both array comparisons and is not counted as two independent
protein replications. No joint-FDR claim is assigned to the three-list overlap.

The primary joint criterion remains the original directional-conjunction global
BH < 0.05 across all twelve comparison/tissue families. Family BH and global BY,
complete-case and single-transcript sensitivities are retained and distinguished
from the broader **assay-wise supported** candidate screen. A significant assay
intersection is not relabelled as a new joint-FDR discovery set. A missing
sensitivity estimate is unavailable evidence, not a zero effect or a failed test.

Names use supplied array symbols and aliases. If no array symbol exists, a
meaningful supplied protein gene label is used, excluding clone/ORF-style names
and bare AGI labels. Otherwise the supplied functional description and AGI remain.
`Name_source` identifies the source; current official nomenclature is not claimed.

## Outputs and reproduction

Set external roots and restore the locked environment as in the README, then run:

```sh
Rscript --vanilla scripts/reproduce_supporting.R
```

| Output ID | Contents |
|---|---|
| `table_supporting_audit` | All primary comparisons, names, profile metrics, original statistics, sensitivity results and selection flags |
| `table_supporting_ranked` | Every supported Array 2/protein comparison, deterministically ordered |
| `table_supporting_dual` | Same-tissue genes supported in both salt and genotype comparisons; full effects and uncertainty for each contrast |
| `table_supporting_three` | Salt support across both arrays and the corresponding tissue proteome |
| `table_supporting_summary` | Eligible and supported counts by original comparison/tissue family, with stricter joint-evidence counts |
| `table_supporting_signal` | Every plotted biological-pool value and observed/missing flag for the leading dual candidates |
| `report_supporting_trends` | Five-page candidate report with methods, ranked tables and original global-joint results |
| `figure_supporting_profiles` | Leading dual candidates: Array 2 pools, corresponding protein pools, and four-condition mean scatterplots; PDF/SVG/TIFF/PNG |

Each ID also works through `Rscript --vanilla scripts/reproduce.R <ID>`. Provenance
records source/input/output hashes, settings, versions and a completion marker.
`supporting_settings()` records the thresholds, ordering and four profile cells.
All data, candidate result values, report files and plots remain local.

`tests/supporting_tests.R` checks that significance in one assay, opposite signs,
nonpositive profiles or threshold equality do not enter the list; verifies
four-mean correlations and equal condition weighting; prevents cross-tissue
salt/genotype intersections; and checks that ranking does not change q values.
