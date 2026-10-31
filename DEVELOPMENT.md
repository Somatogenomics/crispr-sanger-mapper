# Development notes

The browser implementation was ported from a production Python Sanger/CRISPR QC workflow. Key invariants retained in the port:

1. PAM validation is mandatory and strand-aware; finding a protospacer sequence alone is insufficient.
2. A reverse-complemented guide is accepted only when the upstream reference PAM reverse-complements to the requested guide-strand PAM motif.
3. Read orientation is selected by alignment rather than filename/primer assumptions.
4. Both reference bases flanking the expected cleavage boundary are mapped through the gapped alignment into adjacent raw chromatogram basecalls.
5. Mixed-peak evidence is quantified from the four raw/analyzed channel intensities at called-base peak locations; the call requires a cut-relative onset, a sustained after-cut signal, and a pre-cut contrast.
6. Low-quality/low-alignment reads remain inconclusive even when a visual pattern is interesting; manual-review flags are separate from genotype calls.
7. Exact allele size is not inferred from a phase-shifted mixed Sanger trace.

## Known MVP limitations

- The in-browser aligner currently uses seed-localized Smith-Waterman with a linear gap penalty rather than the affine-gap Biopython `PairwiseAligner` used by the Python reference implementation.
- The ABIF reader targets standard sequencing tags (`PBAS`, `PCON`, `PLOC`, `FWO_`, `DATA9–12`) and should be expanded/tested against more ABI instrument generations.
- The classifier is intentionally conservative and should be benchmarked on validated WT, clean-indel, heterozygous, and mosaic datasets before publication-grade use.

## v0.1.1 cut-site classifier audit

A generic post-cut window can spuriously count mixed peaks far from the Cas9 cut.
The new pure-R `cut_localized_mixed_evidence()` only examines aligned bases around
an sgRNA/PAM-validated cleavage boundary. It scores onset (-5..+8), nearby baseline
(-35..-9) and sustained post-cut (+3..+28) in **raw read direction**, and requires
sample-level and local QC. An unbracketed/gapped cut is explicitly inconclusive.

The thresholds are starting points, not calibrated specificity/sensitivity estimates.
Negative controls, orthogonal edited samples, and alignment-indel simulation remain
necessary before a publication-ready accuracy claim. A read that fails a
minimum high-quality window or the original Python QC gate must not be
promoted to a positive automatic call by relaxing quality criteria.
