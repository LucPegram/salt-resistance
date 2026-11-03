# Microarray normalization and regeneration

The primary between-array method is **quantile normalization of log2 Agilent
`gProcessedSignal`**, fitted to noncontrol physical features before transcript
collapse. This is performed independently for Array 1 roots (12 arrays), Array 1
shoots (12 arrays), and Array 2 combined tissue (8 arrays). Every treatment and
genotype in a group participates in the same fit. Conditions are not normalized
separately. Ecotypes and experiments are never pooled for normalization.

## Assessment and decision

The supplied normalized matrices use one multiplicative factor per array. The
source notes identify 75th-percentile Excel scaling. The audit verifies the
per-array scaling directly against the pre-normalization workbooks, without
claiming to reconstruct Excel's exact percentile range, interpolation or rounding.

The assessment compares unnormalized processed signals, noncontrol 75th-percentile
scaling (R type 7, target 1000), quantile normalization and cyclic loess (limma's
fast default, three iterations). It evaluates within-condition residual RMS,
upper-signal residual RMS, repeated-probe dispersion, intensity-dependent
within-condition differences, and positive/negative control dispersion. Biological
replicates contain biological variation, so residual RMS is not pure technical
noise. Equal distributions after quantile normalization are imposed by the method
and are not, by themselves, evidence of its superiority. Discovery counts and
candidate significance do not determine the choice.

Quantile normalization is the primary choice because the combined intensity-bias
and control diagnostics improve over scaling, with broadly maintained or improved
within-condition consistency. It is not superior on every diagnostic: duplicate
spot dispersion can rise slightly, and cyclic loess may have smaller residuals in
some groups. Seventy-fifth-percentile scaling remains a useful sensitivity.
Quantile normalization assumes most of the signal distribution is comparable
within each fitting group; a true global transcriptional shift could violate that
assumption. Absolute RNA content and spike-in concentration records are unavailable,
so these data cannot establish a uniquely optimal biological ground truth.

Controls are excluded from fitting. For held-out diagnostics only, their signals
are mapped through each sample's fitted monotone transformation by interpolation.
Only controls within the noncontrol intensity range in every group sample are
used; corner controls are omitted. These checks assess consistency, not known
concentration accuracy, because spike-in mix records are absent.

The [Agilent one-color manual](https://www.agilent.com/cs/library/usermanuals/Public/G4140-90040_GeneExpression_OneColor_6.7.pdf)
describes noncontrol 75th-percentile scaling. The [limma User's Guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)
provides a single-channel Agilent quantile-normalization workflow. Our input is
vendor-processed signal, not raw foreground; no second background correction is
applied. Two exports retain processed values but lack complete feature/pixel QC.
The preprocessing does not claim to reproduce scanner acquisition, vendor
background estimation or probe remapping specificity.

## Reproduce from immutable original inputs

```sh
export SALT_DATA_ROOT=/absolute/path/to/original/source_data
export SALT_NORMALIZED_ROOT=/absolute/path/to/analysis_inputs/quantile
export SALT_NORMALIZATION_OUTPUT_ROOT=/absolute/path/to/normalization_diagnostics
Rscript --vanilla tests/normalization_tests.R
Rscript --vanilla scripts/renormalize_microarrays.R
export SALT_DATA_ROOT="$SALT_NORMALIZED_ROOT"
export SALT_OUTPUT_ROOT=/absolute/path/to/quantile_reanalysis
Rscript --vanilla scripts/reproduce.R all
```

All three roots must be separate from the code and from one another. The new input
root links unchanged files to their original absolute paths and replaces only the
three active RNA matrices. Links require the original root to remain available.
Do not edit through those links. A nonempty destination without this workflow's
provenance marker is rejected. Original inputs are hashed before/after processing.

Ordinary spots, including unmapped spots, define the normalization distribution;
controls are then omitted from regenerated matrices. Corrected identifiers,
sample labels and downstream transcript/locus collapse rules are preserved. No
additional expression filter, exclusion, missing-value imputation, protein
normalization or qPCR processing is introduced. Regenerated matrix values are
`2^normalized_log2_signal`, because current imports apply log2 themselves.

The diagnostics root contains source-alignment/feature-QC records, method metrics,
comparison with supplied matrices, sample order, distribution figures and hashes.
`NORMALIZATION_PROVENANCE.R` in the new input root identifies source hashes,
normalization settings, exact code hashes, software versions and normalized matrix
hashes. Every registered output includes this normalization record in its own
provenance. All measurements, diagnostics, caches and generated outputs stay
outside this repository. Historical audit modules continue to read explicitly
historical saved summaries; they are not silently relabelled as new inference.

`R/normalization.R` implements this workflow, and
`scripts/renormalize_microarrays.R` is its external-input entry point. This
preprocessing run has its own completion/provenance records and is not an extra
registered manuscript artifact. The existing artifact IDs remain independently
reproducible from the new input root.
