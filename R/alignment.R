# Seed-localized local alignment for amplicon/locus references.
# Default scoring: match +2, mismatch -2, gap open -5, gap extend -1 (affine).
# Alignment is performed in both read orientations.

seed_region <- function(ref, query, k=11L, flank=120L) {
  if (nchar(ref) <= 1500L || nchar(query) < k) return(c(1L, nchar(ref)))
  starts <- unique(round(seq(1, max(1,nchar(query)-k+1), length.out=min(15,max(1,nchar(query)-k+1)))))
  offsets <- integer()
  for (q0 in starts) {
    s <- substr(query, q0, q0+k-1L)
    if (grepl("^[ACGT]+$", s)) {
      hits <- find_all_fixed(s, ref)
      if (length(hits)) offsets <- c(offsets, hits-q0)
    }
  }
  if (!length(offsets)) {
    # progressively shorter seed fallback
    if (k > 7L) return(seed_region(ref, query, k=k-2L, flank=flank))
    return(c(1L, min(nchar(ref), 5000L)))
  }
  center0 <- as.integer(round(median(offsets)))
  lo <- max(1L, center0 + 1L - flank)
  hi <- min(nchar(ref), center0 + nchar(query) + flank)
  c(lo, hi)
}

# Local alignment with affine gaps (Gotoh). Gap of length k scores
# gap_open + (k - 1) * gap_extend, matching the Biopython PairwiseAligner
# convention (open_gap_score applies to the first gap position).
# Cells on each anti-diagonal are independent, so each diagonal is computed
# as one vectorised step; this keeps the cost acceptable under webR.
smith_waterman <- function(ref, query, match=2, mismatch=-2, gap_open=-5, gap_extend=-1) {
  r <- strsplit(ref, "")[[1]]; q <- strsplit(query, "")[[1]]
  nr <- length(r); nq <- length(q)
  if (!nr || !nq) stop("Empty sequence supplied to aligner")
  NR <- nr + 1L
  ncell <- NR * (nq + 1L)
  H <- numeric(ncell)                 # best local score ending at (i,j)
  E <- rep(-Inf, ncell)               # ends with gap in reference (insertion in read)
  F <- rep(-Inf, ncell)               # ends with gap in read (deletion in read)
  PH <- integer(ncell)                # 0 stop, 1 diag, 2 from E, 3 from F
  PE <- logical(ncell)                # TRUE: E extended from E; FALSE: opened from H
  PF <- logical(ncell)
  for (d in 2:(nr + nq)) {
    i <- max(1L, d - nq):min(nr, d - 1L)
    j <- d - i
    k  <- i + j * NR + 1L             # linear index of (i, j), 0-based i/j
    kL <- k - NR                      # (i, j-1)
    kU <- k - 1L                      # (i-1, j)
    kD <- kL - 1L                     # (i-1, j-1)
    e_ext <- E[kL] + gap_extend; e_opn <- H[kL] + gap_open
    PE[k] <- e_ext > e_opn
    E[k] <- pmax(e_ext, e_opn)
    f_ext <- F[kU] + gap_extend; f_opn <- H[kU] + gap_open
    PF[k] <- f_ext > f_opn
    F[k] <- pmax(f_ext, f_opn)
    dg <- H[kD] + ifelse(r[i] == q[j], match, mismatch)
    h <- pmax(0, dg, E[k], F[k])
    H[k] <- h
    # Tie order: diagonal, then gap in reference, then gap in read.
    PH[k] <- ifelse(h <= 0, 0L, ifelse(h == dg, 1L, ifelse(h == E[k], 2L, 3L)))
  }
  best <- max(H)
  kb <- which.max(H)
  bi <- (kb - 1L) %% NR                # 0-based row = reference position
  bj <- (kb - 1L) %/% NR
  i <- bi; j <- bj; state <- 1L        # 1 = H, 2 = E, 3 = F
  ar <- aq <- character(); mr <- mq <- mk <- integer()
  while (i > 0L && j > 0L) {
    k <- i + j * NR + 1L
    if (state == 1L) {
      if (H[k] <= 0) break
      p <- PH[k]
      if (p == 0L) break
      if (p == 2L) { state <- 2L; next }
      if (p == 3L) { state <- 3L; next }
      ar <- c(r[i], ar); aq <- c(q[j], aq)
      mr <- c(i, mr); mq <- c(j, mq); mk <- c(1L, mk)
      i <- i - 1L; j <- j - 1L
    } else if (state == 2L) {
      ar <- c("-", ar); aq <- c(q[j], aq)
      mr <- c(NA_integer_, mr); mq <- c(j, mq); mk <- c(3L, mk)
      if (!PE[k]) state <- 1L
      j <- j - 1L
    } else {
      ar <- c(r[i], ar); aq <- c("-", aq)
      mr <- c(i, mr); mq <- c(NA_integer_, mq); mk <- c(2L, mk)
      if (!PF[k]) state <- 1L
      i <- i - 1L
    }
  }
  m <- cbind(ref=mr, query=mq, kind=mk)
  paired <- ar != "-" & aq != "-"
  mism <- sum(paired & ar != aq)
  aligned_bases <- sum(paired)
  identity <- if (aligned_bases) (aligned_bases - mism) / aligned_bases else NA_real_
  list(score=best, aligned_ref=paste(ar, collapse=""), aligned_query=paste(aq, collapse=""),
       map=m, ref_start=i + 1L, ref_end=bi, query_start=j + 1L, query_end=bj,
       aligned_bases=aligned_bases, mismatches=mism, identity=identity,
       score_per_base=if (aligned_bases) best / aligned_bases else 0)
}

align_orientation <- function(ref, query, ...) {
  candidates <- list()
  for (strand in c("+","-")) {
    oq <- if (strand=="+") query else rc_seq(query)
    reg <- seed_region(ref, oq)
    subref <- substr(ref, reg[1], reg[2])
    a <- smith_waterman(subref, oq, ...)
    # Map local reference coordinates to original reference coordinates.
    if (nrow(a$map)) {
      ii <- is.finite(a$map[,"ref"])
      a$map[ii,"ref"] <- a$map[ii,"ref"] + reg[1] - 1L
    }
    a$ref_start <- a$ref_start + reg[1] - 1L
    a$ref_end <- a$ref_end + reg[1] - 1L
    a$strand <- strand
    a$query_oriented <- oq
    candidates[[strand]] <- a
  }
  if (candidates[["+"]]$score >= candidates[["-"]]$score) candidates[["+"]] else candidates[["-"]]
}

query_to_raw <- function(qi, aln, hq_start, hq_len) {
  if (aln$strand == "+") hq_start + qi - 1L else hq_start + (hq_len - qi)
}

# Raw read index of the last base before the expected cut, in sequencing direction.
# + strand reads run toward increasing reference coordinates, so the last pre-cut base is
# reference position cut_after; - strand reads run the other way, so it is cut_after + 1.
# If that reference base is deleted in the read, the nearest aligned base on the pre-cut
# side is used. If the aligned (high-quality) part of the read stops short of the cut in
# sequencing direction by at most max_extrapolate bases, the position is extrapolated
# without assuming any indel; this happens when mixed peaks after a heterozygous cut lower
# base quality and the quality trim ends just before the cut.
map_cut_to_raw <- function(cut_after, aln, hq_start, hq_len, n_read, max_extrapolate=30L) {
  out <- list(raw=NA_integer_, method="not_covered")
  if (is.null(aln) || !nrow(aln$map) || is.na(cut_after)) return(out)
  m <- aln$map
  paired <- is.finite(m[, "ref"]) & is.finite(m[, "query"])
  if (!any(paired)) return(out)
  pr <- m[paired, "ref"]; pq <- m[paired, "query"]
  if (aln$strand == "+") {
    target <- cut_after
    if (target >= min(pr) && target <= max(pr)) {
      k <- max(which(pr <= target))
      return(list(raw=query_to_raw(pq[k], aln, hq_start, hq_len), method="aligned"))
    }
    if (target > max(pr) && target - max(pr) <= max_extrapolate) {
      k <- which.max(pr)
      raw <- query_to_raw(pq[k], aln, hq_start, hq_len) + (target - pr[k])
      if (raw <= n_read) return(list(raw=as.integer(raw), method="extrapolated"))
    }
  } else {
    target <- cut_after + 1L
    if (target >= min(pr) && target <= max(pr)) {
      k <- min(which(pr >= target))
      return(list(raw=query_to_raw(pq[k], aln, hq_start, hq_len), method="aligned"))
    }
    if (target < min(pr) && min(pr) - target <= max_extrapolate) {
      k <- which.min(pr)
      raw <- query_to_raw(pq[k], aln, hq_start, hq_len) + (pr[k] - target)
      if (raw <= n_read) return(list(raw=as.integer(raw), method="extrapolated"))
    }
  }
  out
}

alignment_gaps_near <- function(aln, cut_after, window=12L) {
  if (!nrow(aln$map) || is.na(cut_after)) return(data.frame())
  rows <- list()
  # Consecutive gap runs in aligned strings, anchored by reference position.
  ar <- strsplit(aln$aligned_ref, "")[[1]]; aq <- strsplit(aln$aligned_query, "")[[1]]
  ref_cursor <- aln$ref_start
  k <- 1L
  while (k <= length(ar)) {
    if (ar[k] != "-" && aq[k] != "-") { ref_cursor <- ref_cursor+1L; k <- k+1L; next }
    kind <- if (aq[k]=="-") "D" else "I"
    s <- k; anchor <- ref_cursor
    while (k <= length(ar) && ((kind=="D" && aq[k]=="-") || (kind=="I" && ar[k]=="-"))) {
      if (ar[k] != "-") ref_cursor <- ref_cursor+1L
      k <- k+1L
    }
    len <- k-s
    dist <- abs(anchor-cut_after)
    if (dist <= window) rows[[length(rows)+1L]] <- data.frame(kind=kind, anchor=anchor, length=len, distance=dist)
  }
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

# Strict two-flank boundary mapping from v0.1.1, compatible with affine-gap alignments.
# A single aligned base or an extrapolated coordinate is never enough for an
# automated mixed-phase call. Coordinates are 1-based in reference and read.
map_ref_to_raw <- function(ref_pos, aln, hq_start, hq_len) {
  if (is.null(aln) || !nrow(aln$map)) return(NA_integer_)
  k <- which(is.finite(aln$map[, "ref"]) & aln$map[, "ref"] == ref_pos &
             is.finite(aln$map[, "query"]))
  if (!length(k)) return(NA_integer_)
  as.integer(query_to_raw(as.integer(aln$map[k[1L], "query"]), aln, hq_start, hq_len))
}

map_cut_boundary <- function(cut_after, aln, hq_start, hq_len) {
  left <- map_ref_to_raw(cut_after, aln, hq_start, hq_len)
  right <- map_ref_to_raw(cut_after + 1L, aln, hq_start, hq_len)
  valid <- is.finite(left) && is.finite(right) && abs(left - right) == 1L
  list(valid=valid, ref_left=cut_after, ref_right=cut_after + 1L,
       raw_pre=if (valid) min(left, right) else NA_integer_,
       raw_post=if (valid) max(left, right) else NA_integer_,
       raw_left=left, raw_right=right,
       reason=if (valid) "Adjacent cut-flanking reference bases mapped" else
         "Both immediately cut-flanking reference bases are not mapped adjacently")
}

aligned_raw_indices <- function(aln, hq_start, hq_len, cut_after, flank=50L) {
  if (is.null(aln) || !nrow(aln$map)) return(integer())
  mm <- aln$map
  valid <- is.finite(mm[, "ref"]) & is.finite(mm[, "query"]) &
    abs(mm[, "ref"] - cut_after) <= flank
  if (!any(valid)) return(integer())
  as.integer(query_to_raw(as.integer(mm[valid, "query"]), aln, hq_start, hq_len))
}

# A homozygous/clean indel can remove the reference base that would ordinarily
# anchor the cut boundary. Detect whether high-quality mapped bases bracket the
# expected cut within a limited reference window; this is ONLY for a candidate
# CLEAN gap, not for mixed-peak classification or WT-like evidence.
map_cut_gap_bracket <- function(cut_after, aln, hq_start, hq_len, flank=12L) {
  invalid <- list(valid=FALSE, raw_pre=NA_integer_, raw_post=NA_integer_,
                  left_ref=NA_integer_, right_ref=NA_integer_, reason="No bracketing gap anchors")
  if (is.null(aln) || !nrow(aln$map)) return(invalid)
  mm <- aln$map
  paired <- is.finite(mm[, "ref"]) & is.finite(mm[, "query"])
  if (!any(paired)) return(invalid)
  r <- as.integer(mm[paired, "ref"]); q <- as.integer(mm[paired, "query"])
  l <- which(r <= cut_after & r >= cut_after - flank)
  h <- which(r >= cut_after + 1L & r <= cut_after + 1L + flank)
  if (!length(l) || !length(h)) return(invalid)
  li <- l[which.max(r[l])]; hi <- h[which.min(r[h])]
  rl <- query_to_raw(q[li], aln, hq_start, hq_len)
  rh <- query_to_raw(q[hi], aln, hq_start, hq_len)
  valid <- is.finite(rl) && is.finite(rh) && rl != rh
  if (!valid) return(invalid)
  list(valid=TRUE, raw_pre=as.integer(min(rl,rh)), raw_post=as.integer(max(rl,rh)),
       left_ref=r[li], right_ref=r[hi], reason="Gap-region anchors bracket expected cut")
}
