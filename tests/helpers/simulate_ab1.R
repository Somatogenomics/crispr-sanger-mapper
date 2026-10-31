# Test helper: simulate Sanger chromatograms and write them as ABIF (.ab1).
# Used only by the test suite. Writes the tags read_ab1() consumes:
# FWO_1, PBAS2, PCON2, PLOC2, DATA9-DATA12.

abif_dir_entry <- function(tag, number, type, esize, n, data_size, offset_or_inline) {
  con <- rawConnection(raw(), "wb")
  on.exit(close(con))
  writeChar(tag, con, nchars=4L, eos=NULL)
  writeBin(as.integer(number), con, size=4L, endian="big")
  writeBin(as.integer(type), con, size=2L, endian="big")
  writeBin(as.integer(esize), con, size=2L, endian="big")
  writeBin(as.integer(n), con, size=4L, endian="big")
  writeBin(as.integer(data_size), con, size=4L, endian="big")
  if (is.raw(offset_or_inline)) {
    writeBin(c(offset_or_inline, raw(4L - length(offset_or_inline))), con)
  } else {
    writeBin(as.integer(offset_or_inline), con, size=4L, endian="big")
  }
  writeBin(0L, con, size=4L, endian="big")
  rawConnectionValue(con)
}

int16_be <- function(x) {
  con <- rawConnection(raw(), "wb"); on.exit(close(con))
  writeBin(as.integer(round(x)), con, size=2L, endian="big")
  rawConnectionValue(con)
}

write_ab1 <- function(path, seq, qual, ploc, channels, order="GATC") {
  stopifnot(nchar(seq) == length(qual), length(qual) == length(ploc))
  items <- list(
    list(tag="FWO_", num=1L, type=2L, esize=1L, data=charToRaw(order)),
    list(tag="PBAS", num=2L, type=2L, esize=1L, data=charToRaw(seq)),
    list(tag="PCON", num=2L, type=2L, esize=1L, data=as.raw(pmin(pmax(round(qual), 0), 99))),
    list(tag="PLOC", num=2L, type=4L, esize=2L, data=int16_be(ploc))
  )
  ord <- strsplit(order, "")[[1]]
  for (i in 0:3) {
    items[[length(items) + 1L]] <- list(tag="DATA", num=9L + i, type=4L, esize=2L,
                                        data=int16_be(pmin(channels[[ord[i + 1L]]], 32767)))
  }
  header_len <- 128L
  body <- raw(); dirs <- raw()
  for (it in items) {
    n_bytes <- length(it$data)
    n_el <- n_bytes / it$esize
    if (n_bytes <= 4L) {
      dirs <- c(dirs, abif_dir_entry(it$tag, it$num, it$type, it$esize, n_el, n_bytes, it$data))
    } else {
      off <- header_len + length(body)
      dirs <- c(dirs, abif_dir_entry(it$tag, it$num, it$type, it$esize, n_el, n_bytes, off))
      body <- c(body, it$data)
    }
  }
  dir_off <- header_len + length(body)
  hdr_con <- rawConnection(raw(), "wb")
  writeChar("ABIF", hdr_con, nchars=4L, eos=NULL)
  writeBin(101L, hdr_con, size=2L, endian="big")
  hdr <- rawConnectionValue(hdr_con); close(hdr_con)
  root <- abif_dir_entry("tdir", 1L, 1023L, 28L, length(items), length(dirs), dir_off)
  out <- c(hdr, root, raw(header_len - length(hdr) - length(root)), body, dirs)
  writeBin(out, path)
  invisible(path)
}

# Gaussian peak trace for one or more aligned alleles.
# alleles: character vector of sequences in sequencing direction, all starting at the
# same template position. weights: relative allele dosage.
simulate_trace <- function(alleles, weights=rep(1, length(alleles)), spacing=12L, sigma=2.4,
                           q_clean=42, ramp_start=15L, ramp_end=40L, q_scale=1, noise=8,
                           seed=1L) {
  set.seed(seed)
  n <- min(nchar(alleles))
  len <- n * spacing + spacing
  ch <- setNames(lapply(1:4, function(i) rep(0, len)), c("A", "C", "G", "T"))
  x <- seq_len(len)
  w <- weights / sum(weights)
  ploc0 <- (seq_len(n) - 1L) * spacing + spacing %/% 2L   # zero-based peak positions
  for (a in seq_along(alleles)) {
    b <- strsplit(alleles[a], "")[[1]][seq_len(n)]
    h <- runif(n, 700, 1300) * w[a]
    for (i in seq_len(n)) {
      if (!b[i] %in% names(ch)) next
      lo <- max(1L, ploc0[i] + 1L - 4L * spacing); hi <- min(len, ploc0[i] + 1L + 4L * spacing)
      xs <- lo:hi
      ch[[b[i]]][xs] <- ch[[b[i]]][xs] + h[i] * exp(-0.5 * ((xs - (ploc0[i] + 1)) / sigma)^2)
    }
  }
  for (k in names(ch)) ch[[k]] <- pmax(0, ch[[k]] + rnorm(len, 0, noise))
  hmat <- sapply(c("A", "C", "G", "T"), function(k) ch[[k]][ploc0 + 1L])
  called <- c("A", "C", "G", "T")[max.col(hmat, ties.method="first")]
  srt <- t(apply(hmat, 1, sort, decreasing=TRUE))
  ratio <- srt[, 2] / pmax(srt[, 1], 1)
  q <- q_clean * (1 - ratio)
  # Low-quality read starts and ends, as on real capillary reads.
  ramp <- pmin(1, seq_len(n) / ramp_start, (n - seq_len(n) + 1) / ramp_end)
  q <- pmax(2, round(q * ramp * q_scale))
  list(seq=paste(called, collapse=""), qual=q, ploc=ploc0, channels=ch)
}

simulate_ab1 <- function(path, alleles, ...) {
  s <- simulate_trace(alleles, ...)
  write_ab1(path, s$seq, s$qual, s$ploc, s$channels)
  invisible(s)
}
