# This classifier first establishes a unique guide + strand-correct PAM on the
# reference, then calculates the expected SpCas9 cut boundary. Mixed peaks are
# NOT treated as editing evidence unless they are localized to that boundary,
# change from a cleaner pre-cut baseline, and persist downstream in the actual
# sequencing direction. The output is a screening hypothesis, never genotype proof.

empty_result <- function(rec, ref_name, guide, call, meanq=NA_real_, hq_len=0L, reason="") {
  list(summary=data.frame(
    sample=rec$name, reference=ref_name, guide=guide,
    status="FAIL", call=call, read_strand=NA_character_, guide_strand=NA_character_,
    pam=NA_character_, pam_start=NA_integer_, pam_end=NA_integer_,
    hq_len=hq_len, mean_q=meanq, aligned_bases=NA_integer_,
    identity=NA_real_, score_per_base=NA_real_, ref_start=NA_integer_, ref_end=NA_integer_,
    cut_after=NA_integer_, cut_between=NA_character_, cut_covered=FALSE,
    raw_pre_cut=NA_integer_, raw_post_cut=NA_integer_,
    pre_scored_n=NA_integer_, pre_median_ratio=NA_real_,
    onset_mixed_n=NA_integer_, onset_first_offset=NA_integer_,
    post_double_n=NA_integer_, post_scored_n=NA_integer_,
    post_double_fraction=NA_real_, post_median_ratio=NA_real_,
    post_minus_pre_ratio=NA_real_, cut_localized_pattern=FALSE,
    gaps_near_cut=NA_integer_, manual_review=reason,
    stringsAsFactors=FALSE),
    rec=rec, hq=NULL, aln=NULL, guide_hit=NULL, metrics=trace_metrics(rec),
    mismatches=data.frame(), gaps=data.frame(), cut_raw=NA_integer_,
    cut_boundary=NULL, cut_evidence=NULL)
}

analyze_sample_guide <- function(rec, ref_name, ref, guide, pam="NGG",
                                 mott_cutoff=0.05, min_hq_len=20, min_qual=20,
                                 min_aligned_bases=25, min_score_per_base=1.2,
                                 min_mean_q=15,
                                 het_threshold=0.30, min_qual_het=10,
                                 min_double_count=4, min_double_fraction=0.25,
                                 edit_window=12) {
  hq <- pick_hq_window(rec, mott_cutoff, min_qual, min_hq_len, allow_low=TRUE)
  qseq <- if (hq$len>0) substr(rec$seq, hq$start, hq$end) else ""
  meanq <- if (hq$len>0) safe_mean(rec$qual[hq$start:hq$end]) else NA_real_
  guides <- find_guide_hits(ref, guide, pam)
  if (length(guides$hits) != 1L) {
    err <- if (!length(guides$hits)) "GUIDE_OR_PAM_NOT_FOUND" else "AMBIGUOUS_GUIDE"
    ans <- empty_result(rec, ref_name, guide, err, meanq, hq$len,
                        paste(c("Check protospacer and strand-correct PAM", guides$rejected), collapse="; "))
    ans$hq <- hq
    return(ans)
  }
  gh <- guides$hits[[1L]]
  metrics <- trace_metrics(rec)
  aln <- align_orientation(ref, qseq)
  cut <- map_cut_boundary(gh$cut_after, aln, hq$start, hq$len)
  cut_covered <- isTRUE(cut$valid)
  mapped <- aligned_raw_indices(aln, hq$start, hq$len, gh$cut_after)
  ev <- cut_localized_mixed_evidence(metrics, cut$raw_pre, mapped,
                                     threshold=het_threshold, min_qual=min_qual_het,
                                     min_double=min_double_count,
                                     min_double_fraction=min_double_fraction)

  # All reported differences use 1-based reference and original read positions.
  mism <- data.frame()
  if (nrow(aln$map)) {
    rr <- aln$map
    rows <- which(is.finite(rr[,"ref"]) & is.finite(rr[,"query"]))
    if (length(rows)) {
      ref_chars <- strsplit(ref, "")[[1]]
      oq_chars <- strsplit(aln$query_oriented, "")[[1]]
      tmp <- lapply(rows, function(k) {
        rp <- as.integer(rr[k,"ref"]); qp <- as.integer(rr[k,"query"])
        if (ref_chars[rp] == oq_chars[qp]) return(NULL)
        rawi <- if (aln$strand=="+") hq$start+qp-1L else hq$start+(hq$len-qp)
        data.frame(ref_pos=rp, ref_base=ref_chars[rp], read_base=oq_chars[qp], raw_index=rawi,
                   qual=if (rawi>=1 && rawi<=length(rec$qual)) rec$qual[rawi] else NA_real_)
      })
      tmp <- Filter(Negate(is.null), tmp)
      if (length(tmp)) mism <- do.call(rbind,tmp)
    }
  }
  gaps <- alignment_gaps_near(aln, gh$cut_after, edit_window)
  gap_n <- nrow(gaps)

  # QC failures take precedence over any suggestive local trace or gap pattern.
  status <- "PASS"
  call <- "NO_CUT_LOCALIZED_MIXED_SIGNAL"
  manual <- ""
  if (hq$len < min_hq_len || identical(hq$method, "full low-quality fallback") ||
      !is.finite(meanq) || meanq < min_mean_q) {
    status <- "FAIL"; call <- "INSUFFICIENT_READ_QUALITY"
  }
  if (aln$aligned_bases < min_aligned_bases) { status <- "FAIL"; call <- "LOW_ALIGNMENT_COVERAGE" }
  if (aln$score_per_base < min_score_per_base) { status <- "FAIL"; call <- "LOW_ALIGNMENT_SCORE" }
  if (!cut_covered && status=="PASS") { status <- "FAIL"; call <- "CUT_BOUNDARY_NOT_RESOLVED" }
  if (cut_covered && !ev$assessable && status=="PASS") {
    status <- "FAIL"; call <- "INSUFFICIENT_CUT_WINDOW_COVERAGE"
  }
  if (status=="PASS" && gap_n>0L) call <- "POSSIBLE_CLEAN_INDEL_NEAR_CUT"
  if (status=="PASS" && isTRUE(ev$pattern)) call <- "POSSIBLE_CUT_LOCALIZED_MIXED_INDEL"

  if (isTRUE(ev$pattern)) {
    manual <- "PAM-validated cut-localized phase change: inspect traces and confirm with orthogonal sequencing"
  } else if (cut_covered && isTRUE(ev$assessable) && ev$post_double >= min_double_count) {
    manual <- "Mixed peaks observed but NOT cut-localized relative to baseline; do not call indel"
  }
  if (gap_n>0L) manual <- paste(c(manual,"Alignment gap near predicted cut; inspect manually"), collapse=if(nzchar(manual)) "; " else "")
  if (nrow(mism) && !nzchar(manual)) manual <- "Reference mismatch(es); compare unedited strain-matched WT"
  if (!cut_covered) manual <- paste(c(cut$reason, manual), collapse=if(nzchar(manual)) "; " else "")
  if (!ev$assessable && cut_covered && !nzchar(manual)) manual <- ev$reason

  list(summary=data.frame(
    sample=rec$name, reference=ref_name, guide=guide,
    status=status, call=call,
    read_strand=aln$strand, guide_strand=gh$strand,
    pam=gh$pam, pam_start=gh$pam_start, pam_end=gh$pam_end,
    hq_len=hq$len, mean_q=round(meanq,2), aligned_bases=aln$aligned_bases,
    identity=round(aln$identity,4), score_per_base=round(aln$score_per_base,3),
    ref_start=aln$ref_start, ref_end=aln$ref_end,
    cut_after=gh$cut_after, cut_between=paste0(gh$cut_after,"|",gh$cut_after+1L),
    cut_covered=cut_covered, raw_pre_cut=cut$raw_pre, raw_post_cut=cut$raw_post,
    pre_scored_n=ev$baseline_n,
    pre_median_ratio=if (is.finite(ev$baseline_median)) round(ev$baseline_median,3) else NA_real_,
    onset_mixed_n=ev$onset_n, onset_first_offset=ev$onset_first,
    post_double_n=ev$post_double, post_scored_n=ev$post_n,
    post_double_fraction=if (is.finite(ev$post_fraction)) round(ev$post_fraction,3) else NA_real_,
    post_median_ratio=if (is.finite(ev$post_median)) round(ev$post_median,3) else NA_real_,
    post_minus_pre_ratio=if (is.finite(ev$delta)) round(ev$delta,3) else NA_real_,
    cut_localized_pattern=isTRUE(ev$pattern),
    gaps_near_cut=gap_n, manual_review=manual, stringsAsFactors=FALSE),
    rec=rec, hq=hq, aln=aln, guide_hit=gh, metrics=metrics,
    mismatches=mism, gaps=gaps, cut_raw=cut$raw_pre,
    cut_boundary=cut, cut_evidence=ev)
}

analyse_batch <- function(files, ref_name, ref, guides, pam="NGG", sample_names=NULL, ...) {
  out <- list(); rows <- list(); errors <- list()
  idx <- 1L
  for (file_i in seq_along(files)) {
    f <- files[[file_i]]
    rec <- tryCatch(read_ab1(f), error=function(e)e)
    if (inherits(rec,"error")) { errors[[basename(f)]] <- conditionMessage(rec); next }
    if (length(sample_names)>=file_i && nzchar(sample_names[[file_i]]))
      rec$name <- sub("\\.[Aa][Bb]1$", "", sample_names[[file_i]])
    for (g in guides) {
      ans <- tryCatch(analyze_sample_guide(rec, ref_name, ref, g, pam, ...), error=function(e)e)
      if (inherits(ans,"error")) { errors[[paste(basename(f),g)]] <- conditionMessage(ans); next }
      key <- paste(rec$name, g, sep="||")
      out[[key]] <- ans; rows[[idx]] <- ans$summary; idx <- idx+1L
    }
  }
  list(results=out, summary=if(length(rows)) do.call(rbind,rows) else data.frame(), errors=errors)
}
