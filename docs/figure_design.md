# Manuscript figure presentation

The eight current manuscript figures use one set of base-R drawing helpers in
`R/revision_figures.R`. Models, biological replicates, feature selection,
observations, uncertainty estimates and correction levels are unchanged by
presentation choices. Historical summary figures remain separate outputs.

- Roots are blue and shoots vermilion. Genotype panels use blue for No-0 and
  purple for the overexpressor; the relevant key is shown outside the data area.
- Shape and line style also distinguish groups, providing cues without colour.
- Dose axes mark only the actual sampled concentrations. Small deterministic
  horizontal offsets separate observations; mean lines use the true dose values.
- Protein panels retain all available pool measurements and observed log-scale
  means. Missing values are never drawn as artificial observations. AHA2 group
  specificity and root NIT2 model eligibility remain qualified in figure legends.
- Contrast panels retain pointwise moderated intervals. BH values sit in their
  own aligned column. Neither intervals nor display colour represent an
  additional significance test.
- The resemblance figure retains every feature and extreme value. Its lower
  panel displays the already defined shared-control and all disjoint-control
  correlations. The connecting line spans the disjoint-split range; it is not a
  confidence interval. Small vertical offsets separate the individual splits.
- PCA uses unchanged coordinates. Deterministic label placement reduces
  collisions, and sample keys occupy a separate panel.
- qPCR uses separate root/shoot panels with independently labelled y axes.
  Points remain reported group means and error bars remain the supplied STE.
  Asterisks mark missing STE rather than implying a zero error bar.
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
