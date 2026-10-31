source("R/utils.R")
source("R/alignment.R")

stopifnot(rc_seq("ACGT") == "ACGT")
stopifnot(pam_ok("AGG","NGG"))
stopifnot(!pam_ok("GTC","NGG"))

ref <- "GGGTTAAACCACGTTGCTACGATCCTGATGAGGTTTGGCTAGACCTAGGCCAT"
guide <- "ACGTTGCTACGATCCTGATG"
g <- find_guide_hits(ref, guide, "NGG")
stopifnot(length(g$hits)==1L)
stopifnot(g$hits[[1]]$strand=="+")
stopifnot(g$hits[[1]]$pam=="AGG")

# Invalid reverse-complement guide without a matching guide-strand PAM must not silently pass.
badref <- paste0("AAA",rc_seq(guide),"GTC", "TTT")
bad <- find_guide_hits(badref, guide, "NGG")
stopifnot(length(bad$hits)==0L)

q <- paste0(guide,"AGG")
a <- align_orientation(ref,q)
stopifnot(a$strand=="+")
stopifnot(a$identity > .95)

q2 <- rc_seq(q)
a2 <- align_orientation(ref,q2)
stopifnot(a2$strand=="-")
stopifnot(a2$identity > .95)

cat("core tests passed\n")

source("R/abif.R")
fx <- "tests/fixtures/synthetic_minimal.ab1"
r <- read_ab1(fx)
stopifnot(r$seq == "ACGTA")
stopifnot(identical(as.numeric(r$qual), c(40,35,30,25,20)))
stopifnot(identical(as.numeric(r$ploc), c(2,5,8,11,14)))
stopifnot(r$channel_order == "GATC")
stopifnot(length(r$channels$G) == 20L)
cat("ABIF parser fixture passed\n")

# Cas9 cleavage is a boundary. For reverse sequencing reads, the raw read
# moves toward lower reference coordinates; raw_pre must still precede raw_post.
ref2 <- paste0("AACCGGTT", "ACGTTGCTACGATCCTGATGAGG", "TTAAACCC")
hit <- find_guide_hits(ref2, guide, "NGG")$hits[[1]]
stopifnot(hit$pam == "AGG")
stopifnot(hit$cut_after == hit$guide_start + nchar(guide)-4L)
al <- align_orientation(ref2,ref2)
cb <- map_cut_boundary(hit$cut_after,al,1L,nchar(ref2))
stopifnot(cb$valid, cb$raw_post-cb$raw_pre==1L)
ar <- align_orientation(ref2,rc_seq(ref2))
cb_rev <- map_cut_boundary(hit$cut_after,ar,1L,nchar(ref2))
stopifnot(cb_rev$valid,cb_rev$raw_post-cb_rev$raw_pre==1L)
stopifnot(cb_rev$raw_left > cb_rev$raw_right)

# Synthetic four-channel-derived peak ratios with a cut boundary at 100|101.
source("R/qc.R")
rawpre <- 100L
mock <- data.frame(
  base=rep("A",200),qual=rep(40,200),primary=rep(600,200),
  second_to_primary=rep(0.05,200),stringsAsFactors=FALSE
)
mapped <- seq_len(200)
call_mock <- function(m) cut_localized_mixed_evidence(m,rawpre,mapped)
stopifnot(call_mock(mock)$assessable,!call_mock(mock)$pattern)

# A sustained phase shift initiated at the cut should be highlighted.
shift <- mock
shift$second_to_primary[rawpre + seq.int(1L,30L)] <- 0.54
stopifnot(call_mock(shift)$assessable,call_mock(shift)$pattern)
stopifnot(call_mock(shift)$onset_first==1L)

# Diffuse high signal across the whole trace is NOT a localized edit.
all_noise <- mock
all_noise$second_to_primary[] <- 0.55
stopifnot(!call_mock(all_noise)$pattern)

# Peaks far downstream, but not near the cut, are NOT cut-localized.
far_noise <- mock
far_noise$second_to_primary[rawpre+seq.int(38L,65L)] <- 0.55
stopifnot(!call_mock(far_noise)$pattern)

# Isolated mixed base at the cut is NOT sufficient.
single <- mock
single$second_to_primary[rawpre+1L] <- 0.65
stopifnot(!call_mock(single)$pattern)

# A sustained mixture starting well after the cut is NOT sufficient.
late_shift <- mock
late_shift$second_to_primary[rawpre+seq.int(14L,30L)] <- 0.55
stopifnot(!call_mock(late_shift)$pattern)

# No pre-cut mapped baseline means inconclusive, not WT-like or an edit.
no_pre <- cut_localized_mixed_evidence(shift,rawpre,101:200)
stopifnot(!no_pre$assessable,!no_pre$pattern)

# A read with no mapped cleavage boundary can never be called.
missing_cut <- cut_localized_mixed_evidence(shift,NA_integer_,mapped)
stopifnot(!missing_cut$assessable,!missing_cut$pattern)

cat("cut-relative mixed-peak regression cases passed\n")

# Valid reverse-strand guide: reverse complement of a plus-strand reference
# must have CCN upstream in reference orientation, verified as NGG on the guide strand.
minus_ref <- rc_seq(ref2)
minus_hits <- find_guide_hits(minus_ref, guide, "NGG")$hits
stopifnot(length(minus_hits)==1L)
stopifnot(minus_hits[[1L]]$strand=="-")
stopifnot(minus_hits[[1L]]$cut_after == nchar(ref2)-hit$cut_after)
stopifnot(minus_hits[[1L]]$pam=="AGG")
cat("strand-aware PAM and cut coordinates passed\n")
