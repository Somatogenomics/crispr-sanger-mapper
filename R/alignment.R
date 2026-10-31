# Seed-localized Smith-Waterman alignment. Intended for amplicon/locus references.
# Scoring mirrors the production Python tool approximately: match +2, mismatch -2,
# gap -5 (linear in this browser implementation). Alignment is performed in both orientations.

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

smith_waterman <- function(ref, query, match=2, mismatch=-2, gap=-5) {
  r <- strsplit(ref, "")[[1]]; q <- strsplit(query, "")[[1]]
  nr <- length(r); nq <- length(q)
  if (!nr || !nq) stop("Empty sequence supplied to aligner")
  H <- matrix(0, nrow=nr+1L, ncol=nq+1L)
  P <- matrix(0L, nrow=nr+1L, ncol=nq+1L) # 1 diag, 2 up, 3 left
  best <- 0; bi <- bj <- 1L
  for (i in seq_len(nr)) {
    for (j in seq_len(nq)) {
      d <- H[i,j] + if (r[i]==q[j]) match else mismatch
      u <- H[i,j+1L] + gap
      l <- H[i+1L,j] + gap
      vals <- c(0,d,u,l)
      k <- which.max(vals)
      H[i+1L,j+1L] <- vals[k]
      P[i+1L,j+1L] <- k-1L
      if (vals[k] > best) { best <- vals[k]; bi <- i+1L; bj <- j+1L }
    }
  }
  i <- bi; j <- bj; ar <- aq <- character(); map <- list()
  while (i>1L && j>1L && H[i,j] > 0) {
    p <- P[i,j]
    if (p==1L) {
      ar <- c(r[i-1L], ar); aq <- c(q[j-1L], aq)
      map[[length(map)+1L]] <- c(ref=i-1L, query=j-1L, kind=1L)
      i <- i-1L; j <- j-1L
    } else if (p==2L) {
      ar <- c(r[i-1L], ar); aq <- c("-", aq)
      map[[length(map)+1L]] <- c(ref=i-1L, query=NA, kind=2L)
      i <- i-1L
    } else if (p==3L) {
      ar <- c("-", ar); aq <- c(q[j-1L], aq)
      map[[length(map)+1L]] <- c(ref=NA, query=j-1L, kind=3L)
      j <- j-1L
    } else break
  }
  map <- rev(map)
  m <- if (length(map)) do.call(rbind, map) else matrix(numeric(), ncol=3)
  aligned_ref <- paste(ar, collapse="")
  aligned_q <- paste(aq, collapse="")
  paired <- ar!="-" & aq!="-"
  mism <- sum(paired & ar!=aq)
  aligned_bases <- sum(paired)
  identity <- if (aligned_bases) (aligned_bases-mism)/aligned_bases else NA_real_
  list(score=best, aligned_ref=aligned_ref, aligned_query=aligned_q,
       map=m, ref_start=i, ref_end=bi-1L, query_start=j, query_end=bj-1L,
       aligned_bases=aligned_bases, mismatches=mism, identity=identity,
       score_per_base=if (aligned_bases) best/aligned_bases else 0)
}

align_orientation <- function(ref, query) {
  candidates <- list()
  for (strand in c("+","-")) {
    oq <- if (strand=="+") query else rc_seq(query)
    reg <- seed_region(ref, oq)
    subref <- substr(ref, reg[1], reg[2])
    a <- smith_waterman(subref, oq)
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

map_ref_to_raw <- function(ref_pos, aln, hq_start, hq_len) {
  if (!nrow(aln$map)) return(NA_integer_)
  row <- which(is.finite(aln$map[,"ref"]) & aln$map[,"ref"]==ref_pos & is.finite(aln$map[,"query"]))
  if (!length(row)) return(NA_integer_)
  qi <- as.integer(aln$map[row[1],"query"])
  if (aln$strand=="+") hq_start + qi - 1L else hq_start + (hq_len - qi)
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

# A Cas9 cut is a BOUNDARY between reference positions cut_after and cut_after+1.
# Do not use only one mapped base: this is unsafe for a deletion/insertion at
# the cleavage boundary and for reads sequenced in the reverse direction.
# raw_pre/raw_post are relative to the CHROMATOGRAM sequencing direction,
# independent of the reference and guide strands.
map_cut_boundary <- function(cut_after, aln, hq_start, hq_len) {
  left <- map_ref_to_raw(cut_after, aln, hq_start, hq_len)
  right <- map_ref_to_raw(cut_after+1L, aln, hq_start, hq_len)
  valid <- is.finite(left) && is.finite(right) && abs(left-right)==1L
  list(
    valid=valid,
    ref_left=cut_after, ref_right=cut_after+1L,
    raw_pre=if (valid) min(left,right) else NA_integer_,
    raw_post=if (valid) max(left,right) else NA_integer_,
    raw_left=left, raw_right=right,
    reason=if (valid) "Adjacent cut-flanking bases mapped" else "Cut boundary not bracketed by adjacent mapped read bases"
  )
}

aligned_raw_indices <- function(aln, hq_start, hq_len, cut_after, flank=45L) {
  if (!nrow(aln$map)) return(integer())
  mm <- aln$map
  ok <- is.finite(mm[,"ref"]) & is.finite(mm[,"query"]) &
        abs(mm[,"ref"]-cut_after) <= flank
  if (!any(ok)) return(integer())
  qi <- as.integer(mm[ok,"query"])
  as.integer(if (aln$strand=="+") hq_start+qi-1L else hq_start+(hq_len-qi))
}
