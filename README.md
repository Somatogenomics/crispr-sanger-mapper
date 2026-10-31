# CRISPR Sanger Mapper (v0.2.1)

A browser-local [Shiny/Shinylive](https://posit-dev.github.io/r-shinylive/) prototype
for checking AB1 chromatograms against a user-supplied FASTA, sgRNA(s), and the
strand-correct PAM. **Screening hypotheses, not definitive edit/genotype calls.**

## What is in v0.2.1?

- Dependency-free native R ABIF parser: PBAS, PCON, PLOC, FWO_, DATA9–12.
- Mott high-quality trimming and mean-Phred gating; low-quality fallback never
  upgrades a sample to genotype-like PASS.
- Local alignment in both orientations using an affine-gap Gotoh aligner.
- One unique strand-valid IUPAC PAM+protospacer required before locating cut.
- Exact **two-flank** cut mapping for mixed-phase calls, plus a separate
  gap-bracketed evaluation for cut-spanning clean deletions/insertions.
- Quality-filtered **pre-cut baseline -35..-9**, **onset -5..+8**, and sustained
  **post-cut +3..+28** evidence, all in original sequencing direction.
- An extrapolated/one-sided cut is a *visual review marker only*; never an
  automated mixed indel or WT-like call.
- Chromatogram with guide/PAM/cut markers, secondary/primary peak ratios,
  aligned mismatches, annotated gaps and batch CSV download.
- Command-line runner (`tools/run_batch.R`) outputs TSV and run parameters.
- Synthetic (nonexperimental) AB1 fixtures and tests; GitHub Actions workflow.

## Calls and their meaning

| Call | Meaning | Status |
|---|---|---|
| `GUIDE_OR_PAM_NOT_FOUND` / `AMBIGUOUS_GUIDE` | PAM invalid or multiple guide sites | FAIL |
| `INSUFFICIENT_READ_QUALITY` | Short/low-Q read or fallback to untrimmed low quality | FAIL |
| `LOW_ALIGNMENT_COVERAGE` / `LOW_ALIGNMENT_SCORE` | Mapping insufficient | FAIL |
| `SGRNA_NOT_COVERED` | Cut not mapped or gap bracket unavailable | FAIL |
| `INSUFFICIENT_CUT_WINDOW_COVERAGE` | Missing quality-filtered baseline, onset or post window | FAIL |
| `MIXED_SIGNAL_NOT_CUT_LOCALIZED` | Secondary peaks fail near-cut phase-change criterion | FAIL |
| `CUT_BOUNDARY_NOT_RESOLVED` | One-sided cut or unresolved clean-gap bracketing | FAIL |
| `POSSIBLE_HET_MOSAIC_INDEL` | Near-cut onset + sustained downstream mixture after valid PAM and QC | PASS (provisional) |
| `POSSIBLE_CLEAN_INDEL_NEAR_CUT` | Single high-quality gap near cut, bracketing anchors | PASS (provisional) |
| `WT_LIKE_NO_HET_INDEL_SIGNAL` | No cut-associated indel signature under tested conditions | PASS (screen only) |

*An automatic `PASS` is a QC/screening status, not validated biological genotype.*
The classifier does not infer exact mutant allele(s), allelic fractions, or
homozygous versus mosaic states from an unresolved chromatogram.

## Run locally

```r
install.packages("shiny")
shiny::runApp(".")
```

## CLI / HPC

```bash
Rscript tools/run_batch.R --ref example/synthetic_reference.fa \
  --guides ATGCTACTCGTACGATTCGA --outdir results --figures reads/*.ab1
```

`--param name=value` overrides any `analyze_sample_guide()` numeric argument;
for example `--param min_mean_q=20`, `--param het_threshold=0.35`.
The CLI logs all used parameters and R session info. Set `CSM_REPO` when invoking
it from outside the repo.

## Tests

```bash
Rscript tests/test_core.R
```

Covers PAM strands and incorrect PAM rejection, 60 randomized affine-gap scores
and traceback assertions, synthetic ABIF parsing, cut-window positive/negative
controls, original read orientation and conservative end-to-end invariants.
Synthetic tests **are not proof of accuracy on real chromatograms**. Before
relying on calls, benchmark against verified edits and matched WT traces.

## Shinylive export and GitHub Pages

```r
install.packages("shinylive")
shinylive::export(".", "site")
```

`.github/workflows/deploy-pages.yml` contains a reusable GitHub Pages workflow.
After pushing to `main`, enable **Settings → Pages → Source: GitHub Actions**.
The app executes under webR in the user's browser; data are processed locally
when used through Shinylive (not when running Shiny in hosted server mode).

## Research-disclosure guard

Only synthetic data/examples are distributed. **Do not commit live experimental
AB1, sgRNA, unpublished phenotype/genotype, or sample metadata** to a public repo.
The synthetic reference is not a biological locus.

## Limitations

- **Matched WT subtraction**, allele deconvolution and automatic allele-fraction
  estimation are not implemented.
- Sensitivity to post-cut mixture is reduced when the quality-trimmed alignment
  fails to cover the cut or downstream trace; this yields `FAIL`/manual review.
- Thresholds and false-positive/negative rates remain uncalibrated.
- The Python reference pipeline is a comparison target; parity is not established.
- Shinylive/webR deployment must be validated on the deployed GitHub Pages app.

MIT License. See `SANGER_REVIEW_STANDARD.md` and `DEVELOPMENT.md`.
