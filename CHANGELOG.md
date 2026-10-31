# Changelog

## v0.1.1 (2026-10-31)

- Remove unpublished gene/guide identifiers and experimental sample details from public examples.
- Add evidence hierarchy, matched-WT design, and benchmarking requirements (`SANGER_REVIEW_STANDARD.md`).

- Fix mixed-peak inference to require unique PAM-validated sgRNA placement and an adjacent mapped cut boundary.
- Evaluate the raw sequencing direction with baseline, cut-onset and near post-cut windows.
- Reject diffuse background, remote peaks, isolated peaks and missing local coverage as cut-edit evidence.
- Enforce high-quality trimmed read gate and minimum mean Phred; preserve failed reads as inconclusive.
- Expose PAM and cleavage positions, onset, and baseline/post-cut metrics in the UI/CSV.
- Annotate guide/PAM/cut in chromatogram.
- Preserve original uploaded sample names instead of Shiny temporary filenames.
- Add regression tests for cut direction, location and non-local noise.
- R/webR execution and genotype sensitivity/specificity remain to be validated.
