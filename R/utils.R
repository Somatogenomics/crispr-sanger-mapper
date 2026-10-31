BASES <- c("A", "C", "G", "T")
IUPAC_ALLOWED <- list(
  A="A", C="C", G="G", T="T",
  R=c("A","G"), Y=c("C","T"), S=c("G","C"), W=c("A","T"),
  K=c("G","T"), M=c("A","C"), B=c("C","G","T"),
  D=c("A","G","T"), H=c("A","C","T"), V=c("A","C","G"),
  N=c("A","C","G","T")
)

rc_seq <- function(x) {
  x <- toupper(gsub("\\s+", "", x))
  from <- strsplit("ACGTRYKMSWBDHVN", "")[[1]]
  to   <- strsplit("TGCAYRMKSWVHDBN", "")[[1]]
  chars <- strsplit(x, "")[[1]]
  lut <- setNames(to, from)
  out <- lut[chars]
  out[is.na(out)] <- chars[is.na(out)]
  paste(rev(out), collapse="")
}

clean_dna <- function(x) toupper(gsub("[^ACGTRYKMSWBDHVN]", "", x))

read_fasta_base <- function(path) {
  lines <- readLines(path, warn=FALSE)
  hdr <- grep("^>", lines)
  if (!length(hdr)) stop("No FASTA records found")
  ends <- c(hdr[-1] - 1L, length(lines))
  recs <- vector("list", length(hdr))
  for (i in seq_along(hdr)) {
    nm <- sub("^>", "", lines[hdr[i]])
    seq_lines <- if (ends[i] > hdr[i]) lines[(hdr[i]+1L):ends[i]] else character()
    seq <- clean_dna(paste(seq_lines, collapse=""))
    recs[[i]] <- list(name=nm, seq=seq)
  }
  recs
}

pam_ok <- function(pam, motif="NGG") {
  pam <- strsplit(toupper(pam), "")[[1]]
  motif <- strsplit(toupper(gsub("\\s+", "", motif)), "")[[1]]
  if (length(pam) != length(motif) || !length(pam)) return(FALSE)
  for (i in seq_along(pam)) {
    allowed <- IUPAC_ALLOWED[[motif[i]]]
    if (is.null(allowed) || !(pam[i] %in% allowed)) return(FALSE)
  }
  TRUE
}

find_all_fixed <- function(pattern, text) {
  if (!nzchar(pattern) || nchar(pattern) > nchar(text)) return(integer())
  starts <- integer()
  pos <- 1L
  repeat {
    m <- regexpr(pattern, substr(text, pos, nchar(text)), fixed=TRUE)[1]
    if (m < 0) break
    hit <- pos + m - 1L
    starts <- c(starts, hit)
    pos <- hit + 1L
    if (pos > nchar(text)) break
  }
  starts
}

find_guide_hits <- function(ref, guide, pam="NGG") {
  ref <- clean_dna(ref)
  guide <- clean_dna(guide)
  pam <- toupper(gsub("\\s+", "", pam))
  if (nchar(guide) < 15) stop("Guide is unusually short; provide the protospacer sequence")
  if (!nzchar(pam)) stop("PAM motif cannot be empty")
  paml <- nchar(pam)
  hits <- list()
  rejected <- character()
  idx <- find_all_fixed(guide, ref)
  if (length(idx)) for (i in idx) {
    p0 <- i + nchar(guide)
    p1 <- p0 + paml - 1L
    pseq <- if (p1 <= nchar(ref)) substr(ref, p0, p1) else ""
    if (pam_ok(pseq, pam)) {
      hits[[length(hits)+1L]] <- list(
        guide=guide, strand="+", guide_start=i, guide_end=i+nchar(guide)-1L,
        pam=pseq, pam_start=p0, pam_end=p1,
        cut_after=i + nchar(guide) - 4L,
        detail=sprintf("guide %d-%d; downstream PAM=%s", i, i+nchar(guide)-1L, pseq)
      )
    } else rejected <- c(rejected, sprintf("+ hit %d rejected: PAM=%s", i, pseq))
  }
  rg <- rc_seq(guide)
  idx <- find_all_fixed(rg, ref)
  if (length(idx)) for (j in idx) {
    p0 <- j - paml
    p1 <- j - 1L
    pref <- if (p0 >= 1L) substr(ref, p0, p1) else ""
    pguide <- if (nzchar(pref)) rc_seq(pref) else ""
    if (pam_ok(pguide, pam)) {
      hits[[length(hits)+1L]] <- list(
        guide=guide, strand="-", guide_start=j, guide_end=j+nchar(guide)-1L,
        pam=pguide, pam_start=p0, pam_end=p1,
        cut_after=j + 2L,
        detail=sprintf("rc(guide) %d-%d; upstream reference PAM=%s; guide-strand PAM=%s",
                       j, j+nchar(guide)-1L, pref, pguide)
      )
    } else rejected <- c(rejected, sprintf("- hit %d rejected: guide-strand PAM=%s", j, pguide))
  }
  list(hits=hits, rejected=rejected)
}

parse_guides <- function(x) {
  z <- unlist(strsplit(x %||% "", "[\\r\\n,;]+"))
  z <- unique(vapply(z, clean_dna, character(1)))
  z[nzchar(z)]
}

`%||%` <- function(x, y) if (is.null(x) || length(x)==0 || all(is.na(x))) y else x

safe_mean <- function(x) if (!length(x) || all(!is.finite(x))) NA_real_ else mean(x[is.finite(x)])
safe_median <- function(x) if (!length(x) || all(!is.finite(x))) NA_real_ else median(x[is.finite(x)])
