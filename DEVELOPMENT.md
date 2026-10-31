# Development and validation notes (v0.2.2)

## Scope

Main branch: v0.2.0 affine-gap aligner / ABIF parser / CLI / GitHub workflow.
Merged: v0.1.1's stricter PAM-cut baseline/onset/sustained signal rule,
cut-flanking base mapping, mean Phred gate and evidence reporting.

## Evidence windows in raw sequencing direction

| Window | Offsets from last pre-cut base |
|---|---|
| Baseline | -35 through -9 |
| Onset | -5 through +8 |
| Sustained after cut | +3 through +28 |

Defaults: min 5 scored baseline bases; min 10 scored downstream bases;
min 4 downstream mixed peaks; min 25% downstream mixed; 2+ in near-cut onset;
post median >= 0.30; post-minus-pre median >= 0.15; pre median <= 0.25;
post >= 1.8x pre. Only mapped bases with primary height >=50 and Phred >=10
are eligible. Ratios are the secondary-to-primary channel peak heights.

These are **starting assumptions**, not validated sensitivity/specificity targets.
Missing mapped cut coverage, baseline, onset, or post coverage produces an
inconclusive `FAIL`, not a WT-like or a mutant call. Extrapolations are visual
only; clean cut-spanning gaps have a separate bracketing rule.

## Benchmark before publication

1. Confirm R syntax and run `Rscript tests/test_core.R` (CI does this).
2. Run a browser smoke test on the exported Shinylive site with synthetic AB1.
3. Compare v0.2.2 output against manually interpreted chromatograms with
   orthogonally verified WT, heterozygous/mosaic indels, homozygous indels,
   SNPs, contamination and failed reads; split by locus and instrument run.
4. Compare against the pinned Python `sanger_crispr_qc.py` implementation;
   record differences in alignment, cut mapping and mixed-peak metrics.
5. Calibrate thresholds on a development set, report sensitivity, specificity,
   ROC/PR with uncertainty, and freeze evaluation on an independent test set.
6. Perform parameter sweeps; register false-positive/false-negative reviews and
   manual reviewer agreement; protect unpublished experimental data.

## R Shiny / native limitations

R's webR browser package availability and performance must be tested after
GitHub Pages deployment. Only base R + `shiny` are required at runtime.
Real ABI traces from multiple basecaller generations and indel mixture levels
are necessary for external validation.

## CLI

`tools/run_batch.R --param min_mean_q=15 --param het_threshold=0.3` writes
summary.tsv, mismatches.tsv, errors.tsv and processing parameters to the
output directory. Runs are separate from public test data.

## v0.2.2 reporting verification

Run `Rscript tests/test_core.R` (including synthetic PDF/PNG/CSV/TSV export
checks), then confirm actual Chrome/Edge and Firefox downloads on the deployed
Shinylive site. Full batch/plot correctness requires a matched WT, known-edited
and no-edit real-chromatogram benchmark. Do not commit experimental AB1 files.
