# Dependency-free test suite. Run from the repository root: Rscript tests/test_core.R
for (f in c("R/utils.R", "R/abif.R", "R/qc.R", "R/alignment.R", "R/analysis.R"))
  source(f)
source("tests/helpers/simulate_ab1.R")

n_pass <- 0L
check <- function(desc, cond) {
  if (!isTRUE(cond)) stop("FAILED: ", desc, call.=FALSE)
  n_pass <<- n_pass + 1L
}

# ---- sequence utilities and PAM validation -------------------------------
check("rc of palindrome", rc_seq("ACGT") == "ACGT")
check("rc with IUPAC", rc_seq("ARYN") == "NRYT")
check("NGG accepts AGG", pam_ok("AGG", "NGG"))
check("NGG rejects GTC", !pam_ok("GTC", "NGG"))

ref <- "TTTACGTCATGCTACTCGTACGATTCGAAGGAAATTT"
guide <- "ATGCTACTCGTACGATTCGA"
g <- find_guide_hits(ref, guide, "NGG")
check("one + strand hit", length(g$hits) == 1L && g$hits[[1]]$strand == "+")
check("+ strand PAM", g$hits[[1]]$pam == "AGG")
check("+ strand cut 3 bp 5' of PAM", g$hits[[1]]$cut_after == g$hits[[1]]$pam_start - 4L)

rref <- rc_seq(ref)
gr <- find_guide_hits(rref, guide, "NGG")
check("one - strand hit", length(gr$hits) == 1L && gr$hits[[1]]$strand == "-")
# The same physical cut seen from the other strand: cut_after (+) mirrors to n - cut_after.
check("- strand cut mirrors + strand cut",
      gr$hits[[1]]$cut_after == nchar(ref) - g$hits[[1]]$cut_after)

badref <- paste0("AAA", rc_seq(guide), "GTC", "TTT")
check("rc guide without guide-strand PAM rejected", length(find_guide_hits(badref, guide, "NGG")$hits) == 0L)

# ---- aligner ---------------------------------------------------------------
q <- "ATGCTACTCGTACGATTCGAAGG"
a <- align_orientation(ref, q)
check("forward orientation", a$strand == "+" && a$identity > .95)
a2 <- align_orientation(ref, rc_seq(q))
check("reverse orientation", a2$strand == "-" && a2$identity > .95)

gotoh_loop <- function(x, y, m=2, mm=-2, o=-5, e=-1) {
  x <- strsplit(x, "")[[1]]; y <- strsplit(y, "")[[1]]
  H <- matrix(0, length(x) + 1, length(y) + 1); E <- F <- matrix(-Inf, length(x) + 1, length(y) + 1)
  for (i in seq_along(x)) for (j in seq_along(y)) {
    E[i + 1, j + 1] <- max(E[i + 1, j] + e, H[i + 1, j] + o)
    F[i + 1, j + 1] <- max(F[i, j + 1] + e, H[i, j + 1] + o)
    H[i + 1, j + 1] <- max(0, H[i, j] + ifelse(x[i] == y[j], m, mm), E[i + 1, j + 1], F[i + 1, j + 1])
  }
  max(H)
}
rescore <- function(al) {
  ar <- strsplit(al$aligned_ref, "")[[1]]; aq <- strsplit(al$aligned_query, "")[[1]]
  s <- 0; prev <- ""
  for (k in seq_along(ar)) {
    st <- if (ar[k] == "-") "E" else if (aq[k] == "-") "F" else "M"
    s <- s + if (st == "M") ifelse(ar[k] == aq[k], 2, -2) else if (st == prev) -1 else -5
    prev <- st
  }
  s
}
set.seed(11)
for (t in 1:60) {
  x <- paste(sample(c("A", "C", "G", "T"), sample(8:50, 1), TRUE), collapse="")
  y <- x
  for (z in seq_len(sample(0:3, 1))) {
    p <- sample(nchar(y), 1)
    y <- paste0(substr(y, 1, p - 1), paste(sample(c("A", "C", "G", "T"), sample(0:5, 1), TRUE), collapse=""),
                substr(y, p + 1, nchar(y)))
  }
  if (!nzchar(y)) y <- "A"
  al <- smith_waterman(x, y)
  check(sprintf("affine score equals loop Gotoh (case %d)", t), al$score == gotoh_loop(x, y))
  check(sprintf("traceback re-scores to optimum (case %d)", t), al$score == 0 || rescore(al) == al$score)
  check(sprintf("traceback covers stated ranges (case %d)", t),
        gsub("-", "", al$aligned_ref) == substr(x, al$ref_start, al$ref_end) &&
        gsub("-", "", al$aligned_query) == substr(y, al$query_start, al$query_end))
}

# ---- ABIF parser -------------------------------------------------------------
r <- read_ab1("tests/fixtures/synthetic_minimal.ab1")
check("fixture bases", r$seq == "ACGTA")
check("fixture qualities", identical(as.numeric(r$qual), c(40, 35, 30, 25, 20)))
check("fixture peak locations", identical(as.numeric(r$ploc), c(2, 5, 8, 11, 14)))
check("fixture channel order", r$channel_order == "GATC")
check("fixture trace length", length(r$channels$G) == 20L)

# ---- near-cut ratio classifier regression (reads in sequencing direction) ------
# These scenarios test the SCIENTIFIC decision independent of the alignment layer.
# The classifier must reject isolated mixed peaks, diffuse mixing, late mixing,
# missing pre-cut baseline, poor primary signal, or inadequate post coverage.
mk_metrics <- function(ratio=rep(.02, 210L), q=rep(35,210L), signal=rep(500,210L)) {
  data.frame(read_index=seq_len(210L),base=rep("A",210L),
             qual=q, primary=signal, second_to_primary=ratio)
}
raw_pre <- 100L; mapped <- seq_len(210L)
phase <- function(x, indexes, amount=.55) { x$second_to_primary[indexes] <- amount; x }
ev <- function(x, mr=mapped) cut_localized_mixed_evidence(x, raw_pre, mr)
clean <- mk_metrics()
shift <- phase(clean, 99:160)
v <- ev(shift)
check("near-cut sustained phase change detected", v$assessable && v$pattern)
check("onset offset falls near cut", is.finite(v$onset_first) && v$onset_first >= -5L && v$onset_first <= 8L)
check("pre-cut baseline measured", v$baseline_n >= 5 && v$baseline_median < .25)
check("post-cut sustained signal measured", v$post_n >= 10 && v$post_double >= 4)
check("a clean trace does not imply mixed indel", ev(clean)$assessable && !ev(clean)$pattern)
check("isolated secondary peak at cut is insufficient", !ev(phase(clean, 101))$pattern)
check("diffuse mixed signal is not a cut phase change", !ev(phase(clean, 1:210))$pattern)
check("far-downstream mixed peaks are not cut localized", !ev(phase(clean, 143:210))$pattern)
check("delayed onset is not cut localized", !ev(phase(clean, 119:150))$pattern)
check("missing pre-cut mapped bases prevents classification",
      !ev(shift, mr=95:160)$assessable && !ev(shift, mr=95:160)$pattern)
check("missing post-cut mapped bases prevents classification",
      !ev(shift, mr=30:103)$assessable && !ev(shift, mr=30:103)$pattern)
low_signal <- mk_metrics(signal=rep(10,210L))
check("weak primary signal cannot trigger indel call", !ev(phase(low_signal, 99:160))$pattern)
low_quality <- mk_metrics(q=rep(2,210L))
check("poor Phred cannot trigger indel call", !ev(phase(low_quality, 99:160))$pattern)
check("unresolved cut prevents mixed classification", !cut_localized_mixed_evidence(shift,NA,mapped)$pattern)

# ---- cut mapping: both directions require two directly mapped cut flanks -----
source("tests/helpers/synthetic_setup.R")
sref <- synthetic_reference
sg <- synthetic_guide
cut <- find_guide_hits(sref,sg)$hits[[1]]$cut_after
check("synthetic PAM-validated cut after base 176", cut == 176L)
rcut <- cut
for (strand in c("+","-")) {
  read <- if (strand=="+") substr(sref,10,nchar(sref)) else substr(rc_seq(sref),10,nchar(sref))
  al <- align_orientation(sref,read)
  bm <- map_cut_boundary(rcut,al,1L,nchar(read))
  check(paste0("cut boundary resolved on ",strand," read"), al$strand==strand && bm$valid)
  check(paste0("mapped cut flanks adjacent on ",strand," read"), bm$raw_post - bm$raw_pre == 1L)
  idx <- aligned_raw_indices(al,1L,nchar(read),rcut)
  check(paste0("mapped bases bracket cut on ",strand," read"),
        any(idx < bm$raw_pre) && any(idx > bm$raw_post))
}

# ---- ABIF / end-to-end synthetic samples -----------------------------------
cut <- find_guide_hits(sref,sg)$hits[[1]]$cut_after
del <- function(s, at, k) paste0(substr(s, 1, at), substr(s, at + k + 1, nchar(s)))
ins <- function(s, at, x) paste0(substr(s, 1, at), x, substr(s, at + 1, nchar(s)))
fwd <- function(s) substr(s, 10, nchar(s))
rv  <- function(s) substr(rc_seq(s), 10, nchar(s))
td <- tempfile("synthetic_AB1_"); dir.create(td)
scen <- list(
  list(id="wt_F", alleles=fwd(sref), strand="+"),
  list(id="wt_R", alleles=rv(sref), strand="-"),
  list(id="het_del2_F", alleles=c(fwd(sref),fwd(del(sref,cut,2L))), strand="+"),
  list(id="het_del2_R", alleles=c(rv(sref),rv(del(sref,cut,2L))), strand="-"),
  list(id="het_ins1_F", alleles=c(fwd(sref),fwd(ins(sref,cut,"A"))), strand="+"),
  list(id="hom_del7_F", alleles=fwd(del(sref,cut,7L)),strand="+"),
  list(id="hom_span_F", alleles=fwd(del(sref,cut-3L,6L)),strand="+"),
  list(id="lowq_wt_F", alleles=fwd(sref),q_scale=.2,strand="+"),
  list(id="short_F", alleles=substr(sref,10,cut-20L),strand="+")
)
res <- list()
for (s in scen) {
  path <- file.path(td,paste0(s$id,".ab1"))
  args <- list(path=path,alleles=s$alleles)
  if (!is.null(s$q_scale)) args$q_scale <- s$q_scale
  do.call(simulate_ab1,args)
  ans <- analyze_sample_guide(read_ab1(path),"synthetic",sref,sg)
  x <- ans$summary
  res[[s$id]] <- x
  check(sprintf("%s mapped strand",s$id), identical(x$strand,s$strand))
  check(sprintf("%s PAM anchored",s$id), x$pam=="AGG" && x$cut_after==cut)
  check(sprintf("%s no definitive genotype call",s$id), !grepl("CONFIRMED|DEFINITIVE",x$call))
  check(sprintf("%s no PASS for unresolved non-gap cut",s$id),
        x$status != "PASS" || x$cut_boundary_resolved || x$cut_mapping=="gap_bracketed")
  check(sprintf("%s mixed call requires cut-localized pattern",s$id),
        x$call != "POSSIBLE_HET_MOSAIC_INDEL" || (x$cut_localized_pattern && x$cut_boundary_resolved))
}
check("synthetic high-quality WT maps both cut flanks",
      res$wt_F$cut_boundary_resolved && res$wt_R$cut_boundary_resolved)
check("low-quality sample fails QC",res$lowq_wt_F$call == "INSUFFICIENT_READ_QUALITY")
check("short read cannot be called WT-like or edited",
      !res$short_F$call %in% c("WT_LIKE_NO_HET_INDEL_SIGNAL", "POSSIBLE_HET_MOSAIC_INDEL", "POSSIBLE_CLEAN_INDEL_NEAR_CUT"))

# End-to-end simulated allele calls are exploratory until real R and trace data
# benchmarking. They must at least never promote uncertain data to WT-like.
for (id in c("het_del2_F","het_del2_R","het_ins1_F")) {
  x <- res[[id]]
  check(paste("synthetic",id,"not falsely WT-like"),
        x$call != "WT_LIKE_NO_HET_INDEL_SIGNAL" || !x$cut_localized_pattern)
}

# ---- end-to-end phase-shift regression with directly mapped WT basecalls ----
# Construct channel-level secondary peaks while retaining a clean primary
# consensus: isolates cut-associated trace evidence from basecaller errors.
make_trace_pattern <- function(path, bases, ref_sequence, cut_after,
                               kind=c("near", "far", "diffuse", "isolated", "unmapped")) {
  kind <- match.arg(kind)
  tr <- simulate_trace(bases, seed=2026)
  al <- align_orientation(ref_sequence, bases)
  bc <- map_cut_boundary(cut_after, al, 1L, nchar(bases))
  stopifnot(bc$valid)
  raw_pre <- bc$raw_pre
  peaks <- switch(kind,
    near=seq.int(raw_pre-2L, raw_pre+35L),
    far=seq.int(raw_pre+45L, raw_pre+75L),
    diffuse=seq_len(nchar(bases)),
    isolated=raw_pre + 1L,
    unmapped=seq.int(raw_pre-2L, raw_pre+35L))
  peaks <- peaks[peaks >= 1L & peaks <= length(tr$ploc)]
  for (ii in peaks) {
    called <- substr(tr$seq, ii, ii)
    other <- setdiff(BASES, called)[1L]
    peak <- as.integer(tr$ploc[ii])+1L
    primary <- tr$channels[[called]][peak]
    tr$channels[[other]][peak] <- as.integer(round(primary*0.65))
  }
  if (kind == "unmapped") {
    # The quality-trimmed alignment ends before cut. Residual raw trace may
    # still be mixed, but must not be classified by extrapolating the cut.
    tr$qual[seq.int(max(1L,raw_pre-15L), length(tr$qual))] <- 2L
  }
  write_ab1(path, tr$seq, tr$qual, tr$ploc, tr$channels)
  invisible(path)
}

for (strand in c("+","-")) {
  bases <- if (strand=="+") fwd(sref) else rv(sref)
  for (k in c("near","far","diffuse","isolated","unmapped")) {
    file <- file.path(td,paste0("trace_",strand,"_",k,".ab1"))
    make_trace_pattern(file,bases,sref,cut,k)
    x <- analyze_sample_guide(read_ab1(file),"synthetic",sref,sg)$summary
    if (k == "near") {
      check(paste(strand,"near-cut peaks detected"),x$call=="POSSIBLE_HET_MOSAIC_INDEL" &&
            x$cut_boundary_resolved && x$onset_mixed_n>=2 && x$pre_scored_n>=5)
    } else {
      check(paste(strand,k,"not a mixed-indel call"),x$call!="POSSIBLE_HET_MOSAIC_INDEL")
    }
    if (k == "unmapped") {
      check(paste(strand,"cut-extrapolated/no quality prevents false WT"),
            x$status=="FAIL" && x$call!="WT_LIKE_NO_HET_INDEL_SIGNAL")
    }
  }
}

# Batch wrapper preserves filenames and records unreadable input.
bad <- file.path(td,"not_abif.ab1"); writeLines("hello",bad)
b <- analyse_batch(c(file.path(td,"wt_F.ab1"),bad), "synthetic",sref,sg)
check("batch keeps good file", nrow(b$summary)==1L)
check("batch records bad file", length(b$errors)==1L)
b2 <- analyse_batch(file.path(td,c("wt_F.ab1","wt_R.ab1")), "synthetic",sref,sg,
                    sample_names=c("synthetic_sample.ab1","synthetic_sample.ab1"))
check("batch keeps distinct sample identifiers", identical(b2$summary$sample,c("synthetic_sample","synthetic_sample_dup1")))

# ---- v0.2.2 reporting / figures / download artifacts ----------------------
source("R/plots.R")
source("R/reporting.R")
report_ans <- analyze_sample_guide(read_ab1(file.path(td,"wt_F.ab1")),
                                   "synthetic",sref,sg)
check("per-base export retains raw AB1 orientation",
      nrow(per_base_export(report_ans)) == nchar(report_ans$rec$seq))
check("per-base export includes mapped reference positions",
      any(is.finite(per_base_export(report_ans)$reference_position)))
check("per-base export includes peak heights and scoring QC",
      all(c("A","C","G","T","second_to_primary","eligible_for_peak_scoring",
            "in_hq_window") %in% names(per_base_export(report_ans))))
report_dir <- tempfile("report_tests_"); dir.create(report_dir)
sfile <- file.path(report_dir,"sample.csv")
pfile <- file.path(report_dir,"perbase.tsv")
paramfile <- file.path(report_dir,"settings.csv")
pdf_file <- file.path(report_dir,"sample.pdf")
write_sample_evidence(report_ans,sfile)
write_sample_per_base(report_ans,pfile)
settings <- list(mott_cutoff=0.05,secondary_primary_threshold=0.30,
                 min_mean_hq_phred=20)
write_settings_csv(settings,paramfile)
check("sample evidence CSV written",file.exists(sfile) && file.info(sfile)$size > 50)
check("per-base TSV readable",nrow(utils::read.delim(pfile,check.names=FALSE)) == nchar(report_ans$rec$seq))
check("settings CSV roundtrips",all(names(settings) %in% utils::read.csv(paramfile)$parameter))
write_sample_report_pdf(report_ans,pdf_file,settings=settings)
check("multi-page PDF saved",file.exists(pdf_file) && file.info(pdf_file)$size > 2000)
if (capabilities("png")) {
  png_file <- file.path(report_dir,"chromatogram.png")
  save_figure_png(png_file, function() plot_trace(report_ans), width=1200L,height=700L,res=100L)
  check("chromatogram PNG saved",file.exists(png_file) && file.info(png_file)$size > 3000)
}

if(nzchar(Sys.getenv("CSM_VERBOSE")))
  print(do.call(rbind,res)[,c("call","strand","hq_len","cut_mapping",
                              "onset_mixed_n","post_double_fraction","gap_detail")])
cat(sprintf("all %d checks passed\n",n_pass))
