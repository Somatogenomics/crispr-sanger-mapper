# CRISPR Sanger Mapper

A browser-local **Shiny/Shinylive** application for QC and screening of Sanger `.ab1` chromatograms around CRISPR cut sites.

The app is designed for the workflow that motivated it: screen candidate edited individuals, align forward/reverse Sanger reads to a strain-matched reference, validate the sgRNA **together with its PAM**, map the expected SpCas9 cleavage position, inspect the raw four-channel chromatogram, and flag evidence consistent with clean or mixed indels.

## Why this tool exists

Simple consensus-sequence alignment can miss the main signal in heterozygous or mosaic Sanger traces. After a length-changing indel, the two alleles become phase-shifted and the chromatogram may turn into sustained mixed peaks after the cut. Conversely, a clean SNP away from the cut may simply reflect a strain/reference difference. This app keeps those cases separate and deliberately avoids calling an exact allele from an unresolved mixed trace.

## Current MVP features

- Reads Applied Biosystems ABIF/`.ab1` files directly in R (no Bioconductor dependency).
- Extracts base calls, Phred qualities, peak locations, channel order, and analyzed `DATA9–DATA12` traces.
- Mott quality trimming with a fallback high-quality window.
- Automatic forward vs reverse-complement reference mapping.
- Multiple FASTA records and multiple sgRNAs per run.
- Strand-aware IUPAC PAM validation (default `NGG`).
- SpCas9 cut-site mapping only after a valid guide+PAM match.
- Per-base primary/secondary trace metrics.
- **Cut-anchored mixed-peak screening** requiring a unique guide+PAM match, mapped cut *boundary*,
  a near-cut onset, a sustained post-cut signal in sequencing direction, and an increase
  over an adjacent pre-cut baseline.
- Independently reports alignment gaps near the expected cut; an alignment gap does not
  establish an exact allele when traces are phase-shifted.
- Batch summary plus interactive per-sample chromatogram/alignment review.
- CSV export.
- Runs with only `{shiny}` plus base R, which keeps the Shinylive/webR dependency surface small.

## Conservative interpretation

Automated calls are **screening/QC**, not proof of genotype.

- `NO_CUT_LOCALIZED_MIXED_SIGNAL`: no cut-centered mixed-phase signature was found in a read that passed QC.
  This is not proof of an unedited genotype or exclusion of a clean indel.
- `POSSIBLE_CUT_LOCALIZED_MIXED_INDEL`: sustained secondary peaks occur after a
  *PAM-validated, directly mapped* cut boundary, begin near that boundary, and are clearly elevated
  relative to the nearby pre-cut trace. This is a screening candidate, not a confirmed edit.
- `POSSIBLE_CLEAN_INDEL_NEAR_CUT`: the alignment contains a gap close to the expected cut.
- `LOW_ALIGNMENT_SCORE`, `LOW_ALIGNMENT_COVERAGE`, `INSUFFICIENT_READ_QUALITY`,
  `CUT_BOUNDARY_NOT_RESOLVED`, or `INSUFFICIENT_CUT_WINDOW_COVERAGE`:
  no interpretable automated indel call from that read.

A low-quality trace may still receive a **manual-review flag** if mixed signal is concentrated at the cut,
  but that flag never upgrades a failed automated QC call to a positive result.
  A read with only 14 high-quality bases would fail the default 20-base minimum and remain **inconclusive**. Exact mosaic/heterozygous alleles should be resolved with independently verified clones or deep amplicon sequencing.

## Cut-anchored phase-change rule (v0.1.1)

For canonical SpCas9, the app first requires **one** exact reference protospacer and an adjacent
strand-correct PAM matching the supplied IUPAC motif (e.g. NGG). The expected cut is
between the reference bases three nucleotides upstream of the PAM. The cut is not
assumed to coincide with a visually mixed base, and the PAM alone is not an indel signal.

Both bases flanking the predicted cut must map to adjacent raw sequencing basecalls.
For reverse chromatograms, raw sequencing order runs from high to low reference
coordinates; all offsets below follow **raw sequencing direction**, not increasing
reference coordinate. The candidate pattern requires:

- A quality-filtered **pre-cut baseline**, raw offsets **-35 to -9** from the cut boundary.
- Evidence that mixed peaks begin **near the cut**, within offsets **-5 to +8**.
- A **sustained** increase across post-cut offsets **+3 to +28**, not a lone mixed peak.
- An increased secondary/primary peak ratio after the cut compared with the baseline.
- Adequate mapping, peak signal, local base quality and usable reads. If any QC fails,
  the automated call remains inconclusive even if the pattern looks compelling.

These are **provisional, manually reviewable thresholds**. Mixed signal in primer regions,
unmapped tails, far downstream of the expected cut, or throughout the entire chromatogram
is not automatically evidence of a Cas9-induced indel. Widening the event window may be
reasonable for larger on-target rearrangements but must be explicitly validated.

The exported CSV reports the guide strand, PAM sequence/coordinates, reference cut
boundary, mapped raw pre/post-cut base indices, baseline and post-cut statistics,
onset offset, and final QC status. Compare edited specimens to genuine unedited
strain-matched controls; variant allele identity and editing mechanism require
orthogonal confirmation.

## Evidence requirements and limitations

The scientific decision tree and publication-readiness validation plan are in
[`SANGER_REVIEW_STANDARD.md`](SANGER_REVIEW_STANDARD.md). **This MVP does not use a
matched WT `.ab1` to calibrate background and does not implement TIDE/ICE/DECODR-style
allele decomposition.** It reports cut-localized *screening flags*, not validated
indel proportions, editing efficiencies, heterozygosity, mosaicism, or knockouts.
A control-aware comparator, affine-gap alignment, and synthetic-mixture performance
benchmarks are necessary before claims of genotype discrimination.

If the cut cannot be bracketed by the alignment because an indel spans it, the
mixed-peak classifier deliberately abstains. A nearby alignment gap may still
be reported as a **manual review flag**. A negative automated result must never
be interpreted as proof of no edit.

## Run locally with Shiny

```r
install.packages("shiny")
shiny::runApp(".")
```

## Export with Shinylive

```r
install.packages("shinylive")
shinylive::export(".", "site")
httpuv::runStaticServer("site")
```

Shinylive runs R in the browser through webR. The repository intentionally avoids `sangerseqR`, `Biostrings`, and other compiled analysis dependencies so uploaded AB1 files can be processed with a small, portable dependency set.

## GitHub Pages

The repository contains `.github/workflows/deploy-pages.yml`, using Posit's official reusable Shinylive deployment workflow. After pushing:

1. Open **Settings → Pages** in the GitHub repository.
2. Set **Source** to **GitHub Actions**.
3. Push to `main` or manually run the workflow.

## Input recommendations

Use an **amplicon- or locus-sized reference FASTA**, not an entire chromosome, when possible. Include an unedited strain-matched WT chromatogram in each experimental batch. Shared high-quality mismatches seen in multiple edited candidates should be checked in WT controls before being interpreted as editing.

For a canonical SpCas9 experiment, enter the 20-nt protospacer without the PAM and leave PAM as `NGG`.

## Example

`example/synthetic_reference.fa` contains a short **artificial sequence**, with an artificial guide and an NGG PAM. No organism-specific guide, unpublished target sequence, or experimental AB1 file is included in the public repository.

## Privacy

When deployed with Shinylive, the app executes in the user's browser. Uploaded chromatograms and FASTA files are processed locally by webR rather than sent to a Shiny analysis server.

## Validation status

This is an early research-tool release. The source includes synthetic regression tests
for the cut-relative phase-change rule (clean trace, true near-cut shift, diffuse noise,
late shift, isolated mixed base, missing baseline and missing cut). The production Python workflow was run separately on private test chromatograms. Those data and findings are not distributed in the repository. **The revised R implementation has not
been executed under R/webR in the build environment**; run `Rscript tests/test_core.R`
plus the browser smoke test prior to releasing it publicly.

The original Python pipeline remains a reference implementation while the browser port is validated across manually reviewed traces and orthogonal amplicon sequencing.

## License

MIT.
