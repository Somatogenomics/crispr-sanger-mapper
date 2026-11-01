# Validation status — CRISPR Sanger Mapper v0.2.1

## Completed in the artifact-building environment

- Rebased v0.2.0's affine-gap aligner / batch CLI / repository workflow.
- Restored and adapted v0.1.1's PAM-cut phase-change decision logic.
- Generated a wholly synthetic 360-bp reference with an inserted synthetic
  protospacer/PAM (guide starts at 160; predicted SpCas9 cut after base 176).
- Added explicit unit-level regression tests for near-cut versus isolated,
  diffuse, far-downstream and delayed mixed peaks, missing baseline/post coverage,
  poor primary signal and low base quality; plus forward/reverse cut mapping and
  synthetic AB1 end-to-end scenarios.
- Checked that expected R source files are present, obvious bracket/quote nesting
  is balanced, and no locus-specific target sequences are present in the repository.
- ZIP archive integrity will be checked before distribution.

## Not yet completed

**Rscript was unavailable in the build environment.** Consequently the
`Rscript tests/test_core.R` suite has *not* been executed against the merged
code. Browser/Shinylive smoke tests and independent real-data validation have
also **not** been performed. This is a release candidate, not a validated
CRISPR genotyping tool. GitHub Actions will run the synthetic R test suite on
push and block deployment if that suite fails; success still must be checked
before public claims of functioning.

## Acceptance steps after checkout

```bash
Rscript tests/test_core.R
```

Check GitHub **Actions → Core R tests** after pushing. Do not claim passing
execution until that workflow reports success. For a browser run:

```r
install.packages(c("shiny", "shinylive"))
shiny::runApp(".")
shinylive::export(".", "site")
```

If tests fail, inspect the error and revise before publishing scientific claims;
never lower quality gates simply to make the fixture pass.

See `SANGER_REVIEW_STANDARD.md` for exact decision definitions and
`DEVELOPMENT.md` for required benchmark design and limitations.

## v0.2.2 browser reporting

Static checks completed in the build environment; R/webR unavailable for direct
execution. A figure showing mixed signal does not
upgrade a QC failure to edited genotype.
