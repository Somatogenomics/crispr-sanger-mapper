# Evidence standard for cut-site-anchored CRISPR Sanger screening

**Scope:** canonical SpCas9, amplicon-based capillary Sanger AB1 traces. The current
R/Shinylive browser MVP performs a conservative cut-localized screening analysis.
This standard describes the *target analytical design*, not a claim that every
feature has already been implemented or benchmarked.

## Biological unit and question

An AB1 is one sequencing reaction; a PCR is one DNA preparation; an individual,
clone, or experimental pool is the biological unit. A chromatogram from one reaction
must not be presented as independent biological confirmation. In mosaic G0/G1
specimens, a Sanger peak mixture alone cannot distinguish germline from somatic
editing or reliably distinguish diploid heterozygosity from cellular mosaicism.

**Question:** is the candidate's sequence evidence compatible with an on-target
SpCas9 repair event near a uniquely PAM-validated cleavage site, over and above
quality artifacts, colony polymorphism and WT background?

## Five independent evidence layers

1. **Amplicon/template:** confirm one intended PCR product of expected size,
   appropriate extraction/NTC controls, no obvious off-target products;
   confirm primer binding specificity against the appropriate reference build.
2. **Reference/guide:** confirm accession/build, protospacer orientation,
   *adjacent strand-correct PAM*, unique guide locus, and the predicted SpCas9
   cleavage boundary (3 nt upstream of PAM). A guide match without a valid PAM
   is a configuration failure, **never a cut-site call**. Handle reference SNPs
   in guide and PAM explicitly against a colony-matched WT.
3. **Raw trace:** evaluate chromatogram peak separation, baseline, peak spacing,
   SNR, dye artifacts, beginning/end read decay, Phred distribution, ambiguous
   bases, and both forward and reverse sequencing traces if available. Avoid
   inferring from consensus-only alignment.
4. **Cut-centered evidence:** anchor control and candidate to a high-quality
   *pre-cut* region; align in sequencing direction. Compare normalized A/C/G/T
   signals in the same mapped reference coordinates; inspect the local baseline,
   near-cut onset, and persistence of phase-shifted signal after the cut. A
   single double peak suggests a substitution/heterozygous site, not a shifted
   indel. Diffuse mixed peaks or peaks far from the cut are not cut-specific.
5. **Independent validation:** replicate PCR plus opposite-direction sequencing,
   matched unedited strain/colony WT, and orthogonal deep amplicon sequencing or
   cloned/allele-resolved products for exact edit sequence and proportion.

## Minimum control for interpretable *mixed indel* inference

**Recommended required inputs for a final probabilistic model:** candidate AB1,
matched WT AB1 generated with the same sequencing primer/chemistry (plus preferably
a WT biological replicate), amplicon FASTA, validated guide and PAM, and primer
orientation. The present v0.1.1 MVP does **not** have matched-WT trace subtraction.
Its uncalibrated flag is labeled *possible* and should not be called positive.

Assess the **change**, not just an absolute mixed-peak ratio:

- Lower secondary-to-primary fluorescence before the cut, in both WT and candidate.
- A local candidate-vs-WT divergence at or near the expected cleavage site.
- Sustained downstream phase incoherence extending tens of bases in the actual
  **sequencing direction**, i.e. toward decreasing genomic coordinates for a
  reverse read on a plus-strand reference.
- A repeatable signal in independently generated AB1 reads, ideally from
  both read directions (the deterioration appears on the opposite side of the
  cut in reference coordinates in reverse vs forward reads).

A cut-spanning deletion may prevent the two adjacent reference bases from both
mapping; that is not a negative result. Inspect nearby gaps and use allele-resolved
sequencing. Position tolerance must be sensitivity-tested, since microhomology,
larger deletions, and aligner ambiguity affect apparent onset.

## Evidence classes: do NOT collapse into a binary WT/KO call

| Class | Required evidence | Claim permitted |
| --- | --- | --- |
| INVALID / INCONCLUSIVE | failed PCR/NTC, poor trace, unverified PAM, cut not covered, or insufficient windows | No inference |
| NON-LOCAL OR DIFFUSE MIXTURE | mixed peaks also present pre-cut or in WT, or onset far from cut | Technical/mixed-template issue for review |
| SEQUENCE VARIANT | local high-quality SNP or clean mismatch, not an indel pattern | Variant observed; origin undetermined |
| PUTATIVE CUT-ASSOCIATED MIXTURE | quality gate + strong cut-anchored change and persistent phase mixture; preferably matched-WT contrast | Compatible with mixed indel; exact alleles unresolved |
| CANDIDATE CLEAN INDEL | coherent gapped alignment supported by raw peak quality and independent sequencing | Candidate indel; coordinate and allele provisional |
| CONFIRMED ALLELE | independent allele-resolving sequencing with read-level support and controls | Exact edit, with detection limits and uncertainty |

Never infer knockout solely from a frameshift-looking Sanger read: confirm
location in an exon/transcript, gene model, allele dosage, and germline
transmission/phenotype as applicable.

## Allele deconvolution (future version, not in v0.1.1)

Implement TIDE/ICE-style candidate allele models: align a reliable pre-cut WT
trace, simulate WT and indel variants around the predicted cut, fit normalized
A/C/G/T signals with constrained NNLS, report fit diagnostics, model residuals,
indel-size limits, background error and sensitivity to analysis window. Use
re-sampling or validated calibration to quantify uncertainty. Support multiple
sgRNAs and large edits separately. If poor quality, ambiguous alignment or
model fit failure, abstain rather than force a fractional genotype.

## Experimental validation before public scientific claims

Benchmark against: matched WT replicates; control/reference polymorphisms;
known single-SNP heterozygotes; clean confirmed +/-1 and larger indels;
synthetic WT/edited mixtures at 0/5/10/25/50/75/100% when technically
feasible; diluted/failed reads; opposite-strand Sanger; duplicate PCRs;
misleading double primers and multiple products. Use hold-out samples and
report confusion matrices, false-positive rates, detection limits, calibration,
and confidence interval coverage; predefine thresholds before validation.

**Interpretation constraint:** trace-based Sanger signals alone cannot
establish allele frequencies with calibrated uncertainty without a validated
model and controls. Mosaic, heterozygous and biallelic calls must remain
unresolved until allele-resolving evidence is available.

## References

- Brinkman EK et al. (2014). *Easy quantitative assessment of genome editing
  by sequence trace decomposition*. Nucleic Acids Research 42:e168.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4267669/
- Conant D et al. (2022). *Inference of CRISPR Edits from Sanger Trace Data*.
  https://github.com/synthego-open/ice
- Bloh KM et al. (2021). *Deconvolution of Complex DNA Repair (DECODR)*.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7898406/
- CRISPResso2: allele-resolved analysis of amplicon reads.
  https://crispresso2.pinellolab.org/help
