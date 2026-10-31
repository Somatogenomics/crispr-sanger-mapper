mott_trim <- function(qual, cutoff=0.05) {
  if (!length(qual)) return(c(1L, 0L))
  score <- cutoff - 10^(-as.numeric(qual)/10)
  best <- cur <- 0
  best_s <- cur_s <- 1L
  best_e <- 0L
  for (i in seq_along(score)) {
    cur <- cur + score[i]
    if (cur < 0) { cur <- 0; cur_s <- i + 1L }
    if (cur > best) { best <- cur; best_s <- cur_s; best_e <- i }
  }
  if (best <= 0 || best_e < best_s) c(1L, 0L) else c(best_s, best_e)
}

longest_quality_window <- function(qual, min_q=20, min_len=20) {
  good <- qual >= min_q
  best <- c(1L, 0L)
  i <- 1L
  while (i <= length(good)) {
    if (!good[i]) { i <- i+1L; next }
    j <- i
    while (j <= length(good) && good[j]) j <- j+1L
    if ((j-i) > (best[2]-best[1]+1L)) best <- c(i, j-1L)
    i <- j
  }
  if ((best[2]-best[1]+1L) < min_len) c(1L,0L) else best
}

pick_hq_window <- function(rec, mott_cutoff=0.05, min_q=20, min_len=20, allow_low=TRUE) {
  w <- mott_trim(rec$qual, mott_cutoff)
  method <- sprintf("Mott %.3g", mott_cutoff)
  if ((w[2]-w[1]+1L) < min_len) {
    z <- longest_quality_window(rec$qual, min_q, min_len)
    if ((z[2]-z[1]+1L) > (w[2]-w[1]+1L)) { w <- z; method <- paste0("longest Q", min_q) }
  }
  if ((w[2]-w[1]+1L) < min_len && allow_low && nchar(rec$seq)) {
    w <- c(1L, nchar(rec$seq)); method <- "full low-quality fallback"
  }
  list(start=w[1], end=w[2], method=method, len=max(0L,w[2]-w[1]+1L),
       fallback=identical(method, "full low-quality fallback"))
}

trace_metrics <- function(rec) {
  n <- nchar(rec$seq)
  out <- data.frame(
    read_index=seq_len(n), base=strsplit(rec$seq, "")[[1]], qual=rec$qual,
    trace_pos=rec$ploc, A=NA_real_, C=NA_real_, G=NA_real_, T=NA_real_,
    primary_base=NA_character_, secondary_base=NA_character_,
    primary=NA_real_, secondary=NA_real_, second_to_primary=NA_real_, aberrant=NA_real_,
    stringsAsFactors=FALSE
  )
  # PLOC is a zero-based trace index in ABIF; R vectors are one-based.
  for (i in seq_len(n)) {
    p <- as.integer(rec$ploc[i]) + 1L
    h <- setNames(rep(NA_real_,4), BASES)
    for (b in BASES) if (p >= 1L && p <= length(rec$channels[[b]])) h[b] <- rec$channels[[b]][p]
    out[i, BASES] <- h
    if (all(is.finite(h)) && sum(h)>0) {
      ord <- order(h, decreasing=TRUE)
      out$primary_base[i] <- names(h)[ord[1]]
      out$secondary_base[i] <- names(h)[ord[2]]
      out$primary[i] <- h[ord[1]]
      out$secondary[i] <- h[ord[2]]
      out$second_to_primary[i] <- if (h[ord[1]]>0) h[ord[2]]/h[ord[1]] else NA_real_
      # Match production code: fraction of total signal not assigned to called base.
      called <- out$base[i]
      out$aberrant[i] <- if (called %in% BASES) 1 - h[called]/sum(h) else NA_real_
    }
  }
  out
}


# Evaluate a PHASE CHANGE localized to a PAM-validated Cas9 cut boundary.
# Offsets are measured in actual chromatogram/sequencing direction; 0 is the
# last raw base before the reference cut, +1 the first base after it.
# Windows are intentionally restricted to the aligned cut neighborhood:
# - baseline: -35..-9 bp (away from the cut but nearby)
# - onset:    -5..+8 bp (tolerance for repair variation around predicted cut)
# - sustained signal: +3..+28 bp (after the cut in read direction)
# An isolated mixed base, diffuse noise, high peaks far from the cut, or no
# sufficiently clean baseline cannot constitute a cut-localized mixed-indel call.
cut_localized_mixed_evidence <- function(metrics, raw_pre, mapped_raw,
                                         threshold=0.30,
                                         min_qual=10L, min_primary=50,
                                         min_baseline=5L, min_post=10L,
                                         min_double=4L,
                                         min_double_fraction=0.25,
                                         min_change=0.15) {
  if (!is.finite(raw_pre) || !length(mapped_raw) || !nrow(metrics)) {
    return(list(assessable=FALSE, pattern=FALSE, reason="No mapped cut boundary",
                baseline_n=0L, baseline_median=NA_real_, onset_n=0L,
                onset_first=NA_integer_, post_n=0L, post_double=0L,
                post_fraction=NA_real_, post_median=NA_real_, delta=NA_real_))
  }
  mapped <- unique(as.integer(mapped_raw))
  idx_for <- function(lo, hi) {
    ix <- seq.int(raw_pre+lo, raw_pre+hi)
    ix <- ix[ix>=1L & ix<=nrow(metrics) & ix %in% mapped]
    if (!length(ix)) return(integer())
    ok <- is.finite(metrics$second_to_primary[ix]) &
          is.finite(metrics$primary[ix]) & metrics$primary[ix]>=min_primary &
          is.finite(metrics$qual[ix]) & metrics$qual[ix]>=min_qual &
          metrics$base[ix] %in% BASES
    ix[ok]
  }
  baseline <- idx_for(-35L,-9L)
  onset <- idx_for(-5L,8L)
  downstream <- idx_for(3L,28L)
  npre <- length(baseline); npost <- length(downstream)
  pre_med <- safe_median(metrics$second_to_primary[baseline])
  post_med <- safe_median(metrics$second_to_primary[downstream])
  onset_hit <- onset[metrics$second_to_primary[onset] >= threshold]
  post_hit <- downstream[metrics$second_to_primary[downstream] >= threshold]
  onset_n <- length(onset_hit); post_double <- length(post_hit)
  post_fraction <- if (npost) post_double/npost else NA_real_
  delta <- if (is.finite(pre_med) && is.finite(post_med)) post_med-pre_med else NA_real_
  assessable <- npre>=min_baseline && npost>=min_post && length(onset)>=4L
  pattern <- assessable && onset_n>=2L && post_double>=min_double &&
             post_fraction>=min_double_fraction &&
             is.finite(pre_med) && pre_med<=0.25 &&
             is.finite(post_med) && post_med>=threshold &&
             is.finite(delta) && delta>=min_change &&
             post_med >= 1.8 * max(pre_med,0.01)
  reason <- if (!assessable) {
    "Insufficient mapped, quality-filtered bases around cut"
  } else if (!pattern) {
    "No sustained cut-centered increase above pre-cut baseline"
  } else {
    "Near-cut onset and sustained downstream phase mixture relative to cleaner pre-cut baseline"
  }
  list(assessable=assessable, pattern=pattern, reason=reason,
       baseline_n=npre, baseline_median=pre_med,
       onset_n=onset_n,
       onset_first=if (onset_n) as.integer(min(onset_hit)-raw_pre) else NA_integer_,
       post_n=npost, post_double=post_double, post_fraction=post_fraction,
       post_median=post_med, delta=delta)
}
