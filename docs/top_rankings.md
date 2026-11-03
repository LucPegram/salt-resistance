# Dataset rankings and treatment figures

These additional descriptive outputs answer the request for three tables of the
20 largest increases and 20 largest decreases, and three figures showing the
10 largest increases and 10 largest decreases. They do not replace the original
transcript-level manuscript analyses or their correction families.

## Comparisons and units

| Dataset | Ranking contrast | Units and displayed treatments |
|---|---|---|
| Array 1 | Col-0, 100 mM minus 0 mM NaCl | Distinct AGI loci; roots and shoots at 0, 50, 75 and 100 mM; three pools per tissue/concentration |
| Array 2 | No-0, 100 mM minus 0 mM NaCl | Distinct AGI loci, combined tissues; both No-0 and NIT2 shown at 0/100 mM; two pools per cell |
| Proteomics | No-0, 75 mM minus 0 mM NaCl | Unambiguously assigned protein groups with reported identifying support; selected tissue shown for both genotypes at 0/75 mM; up to three observed pools per cell |

The salt-response comparison in No-0 was explicitly confirmed by the investigator.
The ranking contrast is recorded in `analysis_settings$toprank`, the result tables
and every provenance record. The second genotype is plotted as context; a selected
No-0 salt response does not establish a genotype-by-salt interaction.

Array gene-level signals first average repeated exact transcript IDs on the log2
scale, then give unique transcript IDs equal weight within each AGI locus. Array 1
uses a locus-level counterpart of the paired categorical tissue-by-dose model,
with `duplicateCorrelation` for roots/shoots from the same harvested pool. The
endpoint contrasts are estimated within that model. Array 2 uses the existing
locus-level four-cell genotype model. Thus the tables contain distinct genes
rather than repeated transcript IDs, and array signals are processed log2 locus
summaries, not raw scanner intensities or resolved isoforms.

Proteomics preserves full group IDs and uses the existing primary observed-data
effect estimates and P values. Ranking additionally requires a single accession,
at least one reported proteotypic sequence, a valid unique AGI assignment and at
least two observed pools in each genotype-treatment cell. Duplicate loci among
retained groups are excluded. No array coverage requirement is imposed on this
standalone protein ranking. The audit retains all originally retained groups and
their identity/observation exclusions; excluded rows have no ranking BH value.

## Ranking and correction

For each dataset, BH correction is applied before ranking. Array 1 corrects both
tissue endpoint contrasts together across all tested loci; Array 2 corrects the
selected contrast across all loci; proteomics corrects both tissues together
across the eligible identified groups. These are explicit ranking-specific
families, distinct from the original manuscript or integration correction
families. Their adjusted P values must not be substituted into those analyses.

Only comparisons with BH q below 0.05 enter the ranking. For a gene observed in
both tissues, the significant comparison with the largest absolute log2 change
is its representative. Ties use q, AGI and tissue deterministically. Each locus
then appears once, in the direction of that representative comparison. Within
each direction the 20 largest absolute log2 effects form the table, and the first
10 form the figure. The code stops if there are fewer than 20 significant eligible
genes in either direction; it never fills a list with nonsignificant genes.

The selected tissue is reported explicitly. `opposite_sign_other_tissue` flags
opposite estimated signs across tissues, without requiring the other estimate to
be significant. A selected tissue-specific response is not whole-plant regulation
or a direct test of root/shoot differences. The top lists and pointwise intervals
are descriptive selections from the same data; they are not an independent
validation or simultaneous inference for the selected genes.

## Gene labels and signal plots

Names and aliases come from the supplied `Stu_MappingDescription.csv`. The first
reported symbol is used when available; all supplied aliases and full functional
descriptions are retained. If no symbol is supplied, the functional name and AGI
are retained rather than inventing a symbol or dropping the gene from the ranking.
`Name_source` records that distinction. Historical annotation is not a new
independent genome-wide remapping or a claim of current official nomenclature.

Each figure contains twenty panels, ten up and ten down. Array 1 shows both
tissues for the selected loci; Array 2 shows both genotypes in combined tissue;
proteomics shows both genotypes in the selected tissue. Points are observed plant
pools; lines connect condition means, not longitudinal or cross-treatment paired
measurements. Signal axes vary by panel. Protein measurements are log2 abundance,
not absolute concentrations, and missing values are never filled for plotting.
The title identifies the selected comparison and its BH value. Very long
functional names can wrap or shorten in graphics; full names remain in CSV.

## Reproduction

Set the external input/output roots and restore the locked environment as in
the repository README, then run:

```sh
Rscript --vanilla scripts/reproduce_toprank.R
```

Each of the fifteen artifacts also works independently:

| Dataset | Forty-row CSV | Two-page table PDF | Twenty-panel figure |
|---|---|---|---|
| Array 1 | `toprank_array1_table` | `toprank_array1_print` | `toprank_array1_figure` |
| Array 2 | `toprank_array2_table` | `toprank_array2_print` | `toprank_array2_figure` |
| Proteomics | `toprank_proteomics_table` | `toprank_proteomics_print` | `toprank_proteomics_figure` |

Use an ID with `Rscript --vanilla scripts/reproduce.R <ID>`. Additional
`toprank_<dataset>_audit` artifacts contain all tested/retained candidates with
selection flags, and `toprank_<dataset>_signal` artifacts contain every plotted
signal, sample label and observed/missing flag. Figures export PDF/SVG, 600 dpi
TIFF and 200 dpi PNG. The table PDFs use rounded display values; CSV retains full
precision. Each artifact has independent input/code/output hashes, session
information and a completion marker. Outputs and annotation data remain outside
the GitHub repository.

`tests/toprank_tests.R` checks BH correction before selection, distinct-gene
ranking, direction/order, missing-symbol handling and paired-model effect
estimates against independently calculated cell means on synthetic inputs.
