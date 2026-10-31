## v0.2.2 — downloadable evidence figures and reports (2026-10-31)

- Added three on-screen and downloadable per-sample figures: four-channel AB1 chromatogram with sgRNA/PAM/cut markers; secondary/primary peak ratio with pre-cut/onset/post-cut windows; Phred and reference-mapping coverage.
- Added a four-page base-graphics PDF sample report, per-base TSV with raw trace peak metrics/reference mapping, single-sample evidence CSV, batch QC PNG, and complete batch CSV/TSV.
- Added a settings CSV manifest; figure thresholds remain frozen to the values used during the Analyze action, so changing input controls cannot silently redraw an old call with new thresholds.
- Applied Posit Community's documented Chrome/Edge Shinylive downloadButton workaround; use ordinary Shiny downloadHandler with browser-local output.
- Every figure reports the QC PASS/FAIL flag; mixed peaks in QC-failed traces are review evidence, never confirmed edits.
- Added synthetic-only report export regression checks; no experimental chromatograms or guide sites are bundled.

**Validation:** source/static checks only in the assembly environment (no R installed). Core R and deployment browser tests remain necessary; no scientific performance benchmark is claimed.

# Release notes

## v0.2.1 — PAM-anchored phase change (2026-10-31)

Merged v0.2.0's affine-gap aligner, synthetic ABIF regression fixtures, CLI batch runner,
per-sample filenames and GitHub workflows with v0.1.1's stringent PAM-anchored
phase-change logic. Key changes:

- Require both immediately cut-flanking reference bases to map adjacently into the raw
  read before automatically calling a **mixed** indel.
- Evaluate quality-filtered local upstream baseline (-35:-9), onset (-5:+8), and
  sustained downstream signal (+3:+28) in the sequencing direction.
- Mixed peaks beyond the near-cut window, diffuse noise, isolated peaks, a missing
  baseline, weak signal and ambiguous PAM do **not** support a CRISPR mixed-indel call.
- Extrapolated cuts are manual-review/plotting only: they never yield an automated
  mutant or WT-like call.
- Preserve separate **clean-gap** candidate detection when a deletion/insertion
  intersects the immediate cut boundary, with flanking anchors and QC requirements.
- Restore mean-HQ-Phred quality gate and detailed PAM/guide/onset evidence columns.
- Synchronize peak highlighting with classifier's secondary/primary ratio metric.
- Replace the locus-specific public example and tests with wholly synthetic DNA.
- Retain v0.2.0 source-code tests and expand with explicit negative-control scenarios.
- Gate GitHub Pages deployment on passing the core R synthetic test suite.

**Caveat:** This is an uncalibrated research prototype; synthetic tests are not a
clinical/genotyping validation, and no exact alleles are inferred from mixed Sanger traces.

## v0.2.0

Affine-gap alignment, batch CLI, simulated chromatogram fixtures, expanded tests.

## v0.1.1

Cut-anchored phase-change QC, guide/PAM report fields, and conservative screening.
