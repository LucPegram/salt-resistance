# Adding a manuscript table, figure, or panel

1. Define the scientific comparison, experimental unit, preprocessing, required
   inputs, and statistical method. Resolve missing metadata rather than guessing.
2. Add any new input specification in `config/inputs.R`. Keep paths relative to
   `SALT_DATA_ROOT`; keep the measurements outside the repository.
3. Put shared calculations in `R/`. Add a builder accepting `inputs`, `settings`,
   and `destination`, with no dependence on objects left in an interactive session.
4. Add a registry entry in `config/artifacts.R` with a stable manuscript ID,
   required input IDs, builder, exact output filename, and status. Register
   individual panels separately if they must be independently regenerated.
5. Document the single-output command and the method in the README and decision
   record. Describe the figure's error bars, transformations, tested contrasts,
   correction for multiple testing, and filtering in its final manuscript legend.
6. Add meaningful checks using artificial examples. Run the output independently
   from a fresh R session and verify its local provenance and visual/numeric output.
7. Update the dependency lock deliberately if new packages are required. Verify
   restoration in a clean environment before claiming a reproducible release.
8. Run the code-only check on the exact files to be committed. Commit code and
   documentation, never generated tables/figures, source measurements, notebooks
   with embedded results, or local run/provenance files.

Manuscript releases should be tagged with a Git version after the complete
registry reproduces every final table, figure, supplementary item, and required
panel from the documented external inputs.
