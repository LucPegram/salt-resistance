# Manuscript figure presentation

The seven current manuscript figures use one set of base-R drawing helpers in
`R/revision_figures.R`. Models, biological replicates, feature selection,
observations, uncertainty estimates and correction levels are unchanged by
presentation choices. Figure 1D was subsequently changed at the investigator's
request from HKT1 to KUP4/POT3 (AT4G23640), using the same prespecified
transcript selection rule and original paired dose model. Supplementary Table
S13 follows the four currently displayed transcripts; the full genome-wide
test results retain HKT1. Historical summary figures remain separate outputs.

- Roots are blue and shoots vermilion. Genotype panels use blue for No-0 and
  purple for the overexpressor; the relevant key is shown outside the data area.
- Shape and line style also distinguish groups, providing cues without colour.
- Dose axes mark only the actual sampled concentrations. Small deterministic
  horizontal offsets separate observations; mean lines use the true dose values.
- Dose-panel subtitles also show the selected transcript tissue-by-dose BH
  value from correction across all 29,391 transcript IDs.
- Figure 2 overlays all available pool measurements on boxplots of observed
  log2 abundance. Boxes show medians and Tukey hinges; whiskers reach the most
  extreme observations within 1.5 times the hinge spread. Outlying values are
  displayed once as pool points. A cell with one observed pool has a point only,
  without a box; empty cells remain empty. The prior mean bars are replaced by
  these distribution summaries. Observed pool counts remain below each cell.
  Missing values are never drawn as artificial observations. AHA2 group
  specificity and root NIT2 model eligibility remain qualified in figure legends.
  This presentation change does not alter observations, fitted models or tests.
- Contrast panels retain pointwise moderated intervals. BH values sit in their
  own aligned column. Neither intervals nor display colour represent an
  additional significance test.
- Figure 4 replaces the former pre-stress resemblance graphic with five PCA
  panels: Array 1 roots and shoots, Array 2 combined tissue, and proteomics roots
  and shoots. All treatment/genotype combinations are fitted together within
  each panel. The existing centered, unscaled complete-feature coordinates are
  retained; no missing values are filled. PCA is not fitted separately within
  two- or three-pool cells. The axes give explained variance; subtitles give
  feature and pool counts. Every pool is labelled. The key explicitly identifies
  all four Col-0 doses and all four No-0/NIT2 by control/salt combinations.
  No confidence ellipse, sample exclusion or significance claim is added.
  `table_manuscript_pca` exports every coordinate with dataset, genotype, dose,
  tissue, sample and pool labels. The former resemblance output stays available
  as a historical registered artifact; its numeric sensitivity results remain
  in Table 4. PCA is promoted from the supplement, where qPCR and pathway figures
  are renumbered S1 and S2 to avoid duplication.
- qPCR uses separate root/shoot panels with independently labelled y axes.
  Points remain reported group means and error bars remain the supplied STE.
  Short error bars can be covered by the mean symbols; no missing error is
  interpreted as zero uncertainty. Every selected candidate has a reported STE.
  Displaying them in separate panels does not reconstruct biological replication.
- Integration scatterplots retain every eligible locus and use the same axis
  limits for root and shoot within each comparison column. Candidate labels and
  joint-BH symbols are separate encodings. No fitted common biological slope or
  sample pairing is implied.
- Focused pathway panels show all planned terms and contrasts, including
  untestable terms. Symbols describe joint relative enrichment, with direction,
  evidence and multiplicity explained in the legend; none measures activation.

## Export and reproduction

Every current figure command writes vector PDF and SVG, a losslessly compressed
600-dpi TIFF, and a 200-dpi PNG preview beside its primary PDF. The completion
marker requires every registered format, and provenance hashes all four files.
Use PDF for scalable figure assembly and TIFF where a raster submission is needed.
Cairo graphics support is required; no extra plotting packages are installed.

```sh
Rscript --vanilla scripts/reproduce_figures.R
```

Set `SALT_DATA_ROOT` and `SALT_OUTPUT_ROOT` first; model caching remains optional.
Individual commands in `docs/manuscript_outputs.md` still work independently.
All rendered images and figures are generated locally and excluded from Git.
