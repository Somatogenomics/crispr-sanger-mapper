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
