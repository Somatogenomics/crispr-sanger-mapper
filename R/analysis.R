# CRISPR Sanger Mapper v0.2.1 — conservative, PAM-anchored screening.
# Preserve v0.2.0's affine-gap aligner and batch runner, restoring v0.1.1's
# near-cut onset / baseline / downstream phase-change logic. An extrapolated
# cut is a visual aid only and NEVER supports an automated edit or WT-like call.

empty_summary <- function(rec, ref_name, guide, call, hq, meanq, manual) {
  data.frame(
    sample=rec$name, reference=ref_name, guide=guide, status="FAIL", call=call,
    strand=NA_character_, guide_strand=NA_character_, pam=NA_character_,
    pam_start=NA_integer_, pam_end=NA_integer_,
    hq_method=hq$method, hq_len=hq$len, mean_q=round(meanq, 2),
    aligned_bases=NA_integer_, identity=NA_real_, score_per_base=NA_real_,
    ref_start=NA_integer_, ref_end=NA_integer_, cut_after=NA_integer_,
    cut_between=NA_character_, cut_mapping="not_covered", cut_raw=NA_integer_,
    cut_boundary_resolved=FALSE, pre_scored_n=NA_integer_,
    pre_double_fraction=NA_real_, pre_median_ratio=NA_real_,
    onset_mixed_n=NA_integer_, onset_first_offset=NA_integer_,
    post_double_n=NA_integer_, post_scored_n=NA_integer_,
    post_double_fraction=NA_real_, post_median_ratio=NA_real_,
    post_minus_pre_ratio=NA_real_, cut_localized_pattern=FALSE,
    gaps_near_cut=NA_integer_, gap_detail="", manual_review=manual,
    stringsAsFactors=FALSE)
}

analyze_sample_guide <- function(rec, ref_name, ref, guide, pam="NGG",
                                 mott_cutoff=0.05, min_hq_len=20L, min_qual=20L,
                                 min_mean_q=15, min_aligned_bases=25L,
                                 min_score_per_base=1.2,
                                 het_threshold=0.30, min_qual_het=10L,
                                 min_double_count=4L, min_double_fraction=0.25,
                                 min_baseline_scored=5L, min_post_scored=10L,
                                 min_primary_signal=50, min_ratio_change=0.15,
                                 edit_window=12L, max_extrapolate=30L,
                                 match=2, mismatch=-2, gap_open=-5, gap_extend=-1) {
  hq <- pick_hq_window(rec, mott_cutoff, min_qual, min_hq_len, allow_low=TRUE)
  low_quality <- isTRUE(hq$fallback) || hq$len < min_hq_len
  qseq <- if (hq$len > 0L) substr(rec$seq, hq$start, hq$end) else ""
  meanq <- if (hq$len > 0L) safe_mean(rec$qual[hq$start:hq$end]) else NA_real_
  low_quality <- low_quality || !is.finite(meanq) || meanq < min_mean_q
  metrics <- trace_metrics(rec)
  guides <- find_guide_hits(ref, guide, pam)
  if (length(guides$hits) != 1L) {
    call <- if (!length(guides$hits)) "GUIDE_OR_PAM_NOT_FOUND" else "AMBIGUOUS_GUIDE"
    s <- empty_summary(rec, ref_name, guide, call, hq, meanq,
                       paste(c("Verify guide and strand-correct PAM", guides$rejected), collapse="; "))
    return(list(summary=s, rec=rec, hq=hq, aln=NULL, guide_hit=NULL,
                metrics=metrics, mismatches=data.frame(), gaps=data.frame(),
                cut_raw=NA_integer_, cut_boundary=NULL, cut_evidence=NULL))
  }
  gh <- guides$hits[[1L]]
  if (!nzchar(qseq)) {
    s <- empty_summary(rec, ref_name, guide, "INSUFFICIENT_READ_QUALITY", hq, meanq, "")
    s$guide_strand <- gh$strand; s$pam <- gh$pam
    s$pam_start <- gh$pam_start; s$pam_end <- gh$pam_end
    s$cut_after <- gh$cut_after; s$cut_between <- paste0(gh$cut_after, "|", gh$cut_after + 1L)
    return(list(summary=s, rec=rec, hq=hq, aln=NULL, guide_hit=gh,
                metrics=metrics, mismatches=data.frame(), gaps=data.frame(),
                cut_raw=NA_integer_, cut_boundary=NULL, cut_evidence=NULL))
  }

  aln <- align_orientation(ref, qseq, match=match, mismatch=mismatch,
                           gap_open=gap_open, gap_extend=gap_extend)
  # 'candidate_cut' may be extrapolated for plotting/manual review, but an
  # automated phase-change call requires the stricter two-flank boundary.
  candidate_cut <- map_cut_to_raw(gh$cut_after, aln, hq$start, hq$len,
                                  nchar(rec$seq), max_extrapolate=max_extrapolate)
  cut <- map_cut_boundary(gh$cut_after, aln, hq$start, hq$len)
  cut_resolved <- isTRUE(cut$valid)
  gaps <- alignment_gaps_near(aln, gh$cut_after, edit_window)
  gap_n <- nrow(gaps)
  gap_detail <- if (gap_n) paste(sprintf("%s%d@%d", gaps$kind, gaps$length, gaps$anchor), collapse=";") else ""
  bracket <- if (gap_n) map_cut_gap_bracket(gh$cut_after, aln, hq$start, hq$len, flank=edit_window)
             else list(valid=FALSE, raw_pre=NA_integer_, raw_post=NA_integer_, reason="No cut gap")
  clean_gap_bracket <- !cut_resolved && gap_n > 0L && isTRUE(bracket$valid)
  raw_for_eval <- if (cut_resolved) cut$raw_pre else if (clean_gap_bracket) bracket$raw_pre else NA_integer_
  mapped <- aligned_raw_indices(aln, hq$start, hq$len, gh$cut_after, flank=50L)
  evidence <- cut_localized_mixed_evidence(
    metrics, raw_for_eval, mapped, threshold=het_threshold, min_qual=min_qual_het,
    min_primary=min_primary_signal, min_baseline=min_baseline_scored,
    min_post=min_post_scored, min_double=min_double_count,
    min_double_fraction=min_double_fraction, min_change=min_ratio_change)
  # Mixed-phase evidence must use an EXACT adjacent cut-flank mapping.
  mixed_pattern <- cut_resolved && isTRUE(evidence$pattern)

  # Mismatches are 1-based in both raw-read and reference coordinates.
  mism <- data.frame()
  if (nrow(aln$map)) {
    mm <- aln$map
    rows <- which(is.finite(mm[, "ref"]) & is.finite(mm[, "query"]))
    if (length(rows)) {
      ref_chars <- strsplit(ref, "")[[1]]
      oq_chars <- strsplit(aln$query_oriented, "")[[1]]
      tmp <- lapply(rows, function(k) {
        rp <- as.integer(mm[k, "ref"]); qp <- as.integer(mm[k, "query"])
        if (ref_chars[rp] == oq_chars[qp]) return(NULL)
        rawi <- query_to_raw(qp, aln, hq$start, hq$len)
        data.frame(ref_pos=rp, ref_base=ref_chars[rp], read_base=oq_chars[qp],
                   raw_index=rawi,
                   qual=if (rawi >= 1L && rawi <= length(rec$qual)) rec$qual[rawi] else NA_real_)
      })
      tmp <- Filter(Negate(is.null), tmp)
      if (length(tmp)) mism <- do.call(rbind, tmp)
    }
  }

  # Prefer a single continuous, high-confidence gap aligned close to the cut;
  # multiple gaps / elevated mixed signal require manual review instead.
  gap_confident <- gap_n == 1L && isTRUE(evidence$assessable) &&
    is.finite(evidence$baseline_median) && evidence$baseline_median < het_threshold &&
    is.finite(evidence$post_median) && evidence$post_median < het_threshold
  gap_located <- cut_resolved || clean_gap_bracket
  if (gap_confident && gap_located) {
    anchor <- if (cut_resolved) cut$raw_pre else bracket$raw_pre
    # The aligned bases bracketing a clean indel must themselves have usable Q.
    nearby <- c(anchor, if (cut_resolved) cut$raw_post else bracket$raw_post)
    gap_confident <- all(nearby >= 1L & nearby <= nrow(metrics) &
                         metrics$qual[nearby] >= min_qual_het)
  }

  # First failure wins. An interesting trace never overrides read/align QC.
  if (low_quality) call <- "INSUFFICIENT_READ_QUALITY"
  else if (aln$aligned_bases < min_aligned_bases) call <- "LOW_ALIGNMENT_COVERAGE"
  else if (!is.finite(aln$score_per_base) || aln$score_per_base < min_score_per_base) call <- "LOW_ALIGNMENT_SCORE"
  else if (!cut_resolved && !clean_gap_bracket) call <- "SGRNA_NOT_COVERED"
  else if (!isTRUE(evidence$assessable)) call <- "INSUFFICIENT_CUT_WINDOW_COVERAGE"
  else if (mixed_pattern) call <- "POSSIBLE_HET_MOSAIC_INDEL"
  else if (gap_confident && gap_located) call <- "POSSIBLE_CLEAN_INDEL_NEAR_CUT"
  else if (!cut_resolved) call <- "CUT_BOUNDARY_NOT_RESOLVED"
  else if (is.finite(evidence$baseline_median) && evidence$baseline_median >= het_threshold) call <- "MIXED_SIGNAL_NOT_CUT_LOCALIZED"
  else if (evidence$post_double >= min_double_count && !isTRUE(evidence$pattern)) call <- "MIXED_SIGNAL_NOT_CUT_LOCALIZED"
  else call <- "WT_LIKE_NO_HET_INDEL_SIGNAL"

  status <- if (call %in% c("POSSIBLE_HET_MOSAIC_INDEL", "POSSIBLE_CLEAN_INDEL_NEAR_CUT",
                          "WT_LIKE_NO_HET_INDEL_SIGNAL")) "PASS" else "FAIL"
  notes <- character()
  if (mixed_pattern) notes <- c(notes, "PAM-anchored near-cut onset with sustained downstream phase shift; exact allele unresolved")
  if (!cut_resolved && is.finite(candidate_cut$raw)) notes <- c(notes,
    "Cut coordinate approximate/extrapolated or one-sided; manual visualization only")
  if (!cut_resolved && gap_n) notes <- c(notes, "Cut-flanking bases not adjacent; inspect gap coordinates manually")
  if (gap_n) notes <- c(notes, "Alignment gap(s) near expected cut; inspect manually")
  if (isTRUE(evidence$assessable) && evidence$post_double >= min_double_count &&
      !isTRUE(evidence$pattern) && !mixed_pattern) notes <- c(notes,
    "Mixed peaks exist but no validated near-cut onset/baseline change")
  if (nrow(mism)) notes <- c(notes, "Reference mismatch(es) present; compare strain-matched WT")
  if (!isTRUE(evidence$assessable)) notes <- c(notes, evidence$reason)
  if (low_quality) notes <- c(notes, "Read quality failed; candidate evidence is not an automated genotype call")

  cut_mapping <- if (cut_resolved) "aligned_two_flanks"
                 else if (clean_gap_bracket) "gap_bracketed"
                 else if (identical(candidate_cut$method, "extrapolated")) "extrapolated"
                 else if (identical(candidate_cut$method, "aligned")) "one_sided_only"
                 else "not_covered"
  cut_raw <- if (cut_resolved) cut$raw_pre else if (clean_gap_bracket) bracket$raw_pre
             else candidate_cut$raw
  pre_fraction <- if (evidence$baseline_n) {
    ix <- intersect(seq.int(max(1L, raw_for_eval - 35L), raw_for_eval - 9L), mapped)
    ix <- ix[ix >= 1L & ix <= nrow(metrics) & is.finite(metrics$second_to_primary[ix]) &
               metrics$qual[ix] >= min_qual_het & metrics$primary[ix] >= min_primary_signal]
    if(length(ix)) mean(metrics$second_to_primary[ix] >= het_threshold) else NA_real_
  } else NA_real_

  list(summary=data.frame(
    sample=rec$name, reference=ref_name, guide=guide, status=status, call=call,
    strand=aln$strand, guide_strand=gh$strand, pam=gh$pam,
    pam_start=gh$pam_start, pam_end=gh$pam_end,
    hq_method=hq$method, hq_len=hq$len, mean_q=round(meanq, 2),
    aligned_bases=aln$aligned_bases, identity=round(aln$identity, 4),
    score_per_base=round(aln$score_per_base, 3),
    ref_start=aln$ref_start, ref_end=aln$ref_end,
    cut_after=gh$cut_after, cut_between=paste0(gh$cut_after, "|", gh$cut_after + 1L),
    cut_mapping=cut_mapping, cut_raw=if (is.finite(cut_raw)) as.integer(cut_raw) else NA_integer_,
    cut_boundary_resolved=cut_resolved, pre_scored_n=evidence$baseline_n,
    pre_double_fraction=if (is.finite(pre_fraction)) round(pre_fraction,3) else NA_real_,
    pre_median_ratio=if (is.finite(evidence$baseline_median)) round(evidence$baseline_median,3) else NA_real_,
    onset_mixed_n=evidence$onset_n, onset_first_offset=evidence$onset_first,
    post_double_n=evidence$post_double, post_scored_n=evidence$post_n,
    post_double_fraction=if (is.finite(evidence$post_fraction)) round(evidence$post_fraction,3) else NA_real_,
    post_median_ratio=if (is.finite(evidence$post_median)) round(evidence$post_median,3) else NA_real_,
    post_minus_pre_ratio=if (is.finite(evidence$delta)) round(evidence$delta,3) else NA_real_,
    cut_localized_pattern=mixed_pattern,
    gaps_near_cut=gap_n, gap_detail=gap_detail, manual_review=paste(unique(notes),collapse="; "),
    stringsAsFactors=FALSE),
    rec=rec, hq=hq, aln=aln, guide_hit=gh, metrics=metrics,
    mismatches=mism, gaps=gaps,
    cut_raw=if (is.finite(cut_raw)) as.integer(cut_raw) else NA_integer_,
    cut_boundary=cut, cut_evidence=evidence, gap_bracket=bracket)
}

# Preserve v0.2.0 batch behavior (upload names, duplicate protection, error logging).
analyse_batch <- function(files, ref_name, ref, guides, pam="NGG", sample_names=NULL, ...) {
  if (is.null(sample_names)) sample_names <- basename(files)
  sample_names <- sub("\.[Aa][Bb]1$", "", sample_names)
  if (anyDuplicated(sample_names)) sample_names <- make.unique(sample_names, sep="_dup")
  out <- list(); rows <- list(); errors <- list(); idx <- 1L
  for (fi in seq_along(files)) {
    f <- files[fi]
    rec <- tryCatch(read_ab1(f), error=function(e)e)
    if (inherits(rec,"error")) { errors[[sample_names[fi]]] <- conditionMessage(rec); next }
    rec$name <- sample_names[fi]
    for (g in guides) {
      ans <- tryCatch(analyze_sample_guide(rec, ref_name, ref, g, pam, ...), error=function(e)e)
      if (inherits(ans,"error")) { errors[[paste(rec$name,g)]] <- conditionMessage(ans); next }
      key <- paste(rec$name, g, sep="||")
      out[[key]] <- ans; rows[[idx]] <- ans$summary; idx <- idx + 1L
    }
  }
  list(results=out,
       summary=if (length(rows)) do.call(rbind,rows) else data.frame(), errors=errors)
}
