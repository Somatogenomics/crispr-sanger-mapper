# Sanger / CRISPR screening evidence standard — v0.2.1

## Decision hierarchy

1. **Reference and PAM:** one unambiguous protospacer with valid IUPAC PAM on the
   correct strand. For SpCas9 the cut is a *boundary*, not a single base.
2. **Chromatogram QC:** meaningful called read, valid trace channels, sufficient
   high-quality Mott window and mean Phred; failures remain failures.
3. **Alignment QC:** read direction determined from sequence, adequate alignment
   span and per-base score; local clean gaps require their own review.
4. **Cut coverage:** directly map both bases flanking expected cut into raw base
   indices. Extrapolation is for manual visualization only.
5. **Mixed-phase onset:** baseline at -35..-9 in sequencing direction must be
   sufficiently clean; at least two mixed positions around -5..+8, sustained
   mixing +3..+28, adequate quality/primary signal and increase over baseline.
6. **Clean-gap evidence:** a single cut-near alignment gap bracketed by aligned,
   quality-supported flanks can be provisionally flagged even if a reference
   base immediately at the cut is absent. Never assign exact genotype without
   confirmation.
7. **WT-like:** absence of a *detected* cut-centered indel signature is not
   absence of editing, nor does it rule out pure SNVs or low-frequency clones.
8. **Control comparisons:** check matched strain WT (preferably same sequencing
   primer and PCR context) for shared mismatches and baseline trace characteristics.
   **Note:** automated matched-WT subtraction is NOT implemented in v0.2.1.
9. **Validation:** confirmation by independent PCR, bidirectional Sanger and/or
   allele-resolving deep amplicon sequencing where warranted.

## Coordinate conventions

- Reference positions and cut-after positions are 1-based.
- A cut_after value of `c` means cleavage between bases `c|c+1`.
- Raw read peak indices are 1-based in original sequencing direction.
- Mixed signal is measured as the **second-highest channel / highest channel**
  at the called peak, *not* as 1 minus the fraction of intensity in the called base.
- `cut_mapping=aligned_two_flanks` is required for mixed-indel calls;
  `gap_bracketed` may support a separate clean-gap call;
  `extrapolated` is manual-review only.

All thresholds are provisional until benchmarks on verified WT, clean-indel,
mixture, and failed-trace ground truths establish error rates and uncertainty.
