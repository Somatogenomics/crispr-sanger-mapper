#!/usr/bin/env Rscript
# Command-line batch run of the same analysis the app performs.
# Intended for scripted runs and for the parity benchmark against the Python
# reference implementation. Uses base R only.
#
# Usage (from the repository root):
#   Rscript tools/run_batch.R --ref ref.fa [--record NAME] --guides G1[,G2] \
#       [--pam NGG] [--outdir results] [--figures] [--param name=value ...] reads/*.ab1
#
# Any analyze_sample_guide() argument can be set with --param, for example
#   --param gap_open=-5 --param gap_extend=-1 --param het_threshold=0.3
#
# Output:
#   <outdir>/04_results/summary.tsv     one row per read x guide
#   <outdir>/04_results/mismatches.tsv  aligned reference mismatches per read x guide
#   <outdir>/04_results/errors.tsv      files that could not be parsed or analysed
#   <outdir>/05_figures/*.png           chromatogram around the cut (with --figures)
#   <outdir>/logs/parameters.tsv        every analysis parameter actually used
#   <outdir>/logs/sessionInfo.txt       R version and platform

args <- commandArgs(trailingOnly=TRUE)
usage <- function() {
  cat("Usage: Rscript tools/run_batch.R --ref ref.fa --guides G1[,G2] [--record NAME] [--pam NGG]",
      "[--outdir results] [--figures] [--param name=value ...] reads/*.ab1\n")
  quit(status=2)
}
opt <- list(pam="NGG", outdir="results", figures=FALSE, params=character(), record=NULL)
files <- character()
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  take <- function() { if (i + 1L > length(args)) usage(); args[i + 1L] }
  if (a == "--ref") { opt$ref <- take(); i <- i + 2L }
  else if (a == "--record") { opt$record <- take(); i <- i + 2L }
  else if (a == "--guides") { opt$guides <- take(); i <- i + 2L }
  else if (a == "--pam") { opt$pam <- take(); i <- i + 2L }
  else if (a == "--outdir") { opt$outdir <- take(); i <- i + 2L }
  else if (a == "--param") { opt$params <- c(opt$params, take()); i <- i + 2L }
  else if (a == "--figures") { opt$figures <- TRUE; i <- i + 1L }
  else if (a %in% c("-h", "--help")) usage()
  else if (startsWith(a, "--")) { cat("Unknown option:", a, "\n"); usage() }
  else { files <- c(files, a); i <- i + 1L }
}
if (is.null(opt$ref) || is.null(opt$guides) || !length(files)) usage()

# Repository location: CSM_REPO if set, otherwise the current directory.
# (Run from the repository root, or set CSM_REPO in a SLURM script.)
repo <- Sys.getenv("CSM_REPO", unset=getwd())
for (f in c("utils.R", "abif.R", "qc.R", "alignment.R", "analysis.R", "plots.R")) {
  p <- file.path(repo, "R", f)
  if (!file.exists(p)) stop("Cannot find ", p, "; run from the repository root or set CSM_REPO")
  source(p)
}

recs <- read_fasta_base(opt$ref)
rec <- if (is.null(opt$record)) recs[[1]] else {
  hit <- Filter(function(z) z$name == opt$record || strsplit(z$name, "\\s+")[[1]][1] == opt$record, recs)
  if (!length(hit)) stop("FASTA record not found: ", opt$record)
  hit[[1]]
}
guides <- parse_guides(opt$guides)
if (!length(guides)) stop("No valid guide sequences supplied")

defaults <- formals(analyze_sample_guide)
defaults <- defaults[setdiff(names(defaults), c("rec", "ref_name", "ref", "guide", "pam"))]
extra <- list()
for (kv in opt$params) {
  parts <- strsplit(kv, "=", fixed=TRUE)[[1]]
  if (length(parts) != 2L || !parts[1] %in% names(defaults)) stop("Bad --param: ", kv)
  extra[[parts[1]]] <- as.numeric(parts[2])
}
used <- utils::modifyList(lapply(defaults, eval), extra)

res <- do.call(analyse_batch, c(list(files=files, ref_name=rec$name, ref=rec$seq,
                                     guides=guides, pam=toupper(opt$pam)), extra))

d_res <- file.path(opt$outdir, "04_results"); d_fig <- file.path(opt$outdir, "05_figures")
d_log <- file.path(opt$outdir, "logs")
for (d in c(d_res, d_log, if (opt$figures) d_fig)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

write.table(res$summary, file.path(d_res, "summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
mm <- do.call(rbind, lapply(names(res$results), function(k) {
  m <- res$results[[k]]$mismatches
  if (!nrow(m)) return(NULL)
  cbind(sample=res$results[[k]]$summary$sample, guide=res$results[[k]]$summary$guide, m)
}))
if (is.null(mm)) mm <- data.frame(sample=character(), guide=character(), ref_pos=integer(),
                                  ref_base=character(), read_base=character(), raw_index=integer(), qual=numeric())
write.table(mm, file.path(d_res, "mismatches.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
write.table(data.frame(item=names(res$errors), error=unlist(res$errors, use.names=FALSE)),
            file.path(d_res, "errors.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

params <- data.frame(parameter=c("reference_file", "reference_record", "reference_length", "pam", "guides",
                                 names(used)),
                     value=c(normalizePath(opt$ref), rec$name, nchar(rec$seq), toupper(opt$pam),
                             paste(guides, collapse=","), vapply(used, function(v) format(v), character(1))))
write.table(params, file.path(d_log, "parameters.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
writeLines(c(paste("run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
             paste("files:", length(files)), capture.output(sessionInfo())),
           file.path(d_log, "sessionInfo.txt"))

if (opt$figures) {
  for (k in names(res$results)) {
    x <- res$results[[k]]
    fn <- gsub("[^A-Za-z0-9_.-]+", "_", paste0(x$summary$sample, "__", x$summary$guide, ".png"))
    png(file.path(d_fig, fn), width=2000, height=800, res=200)
    plot_trace(x, het_threshold=used$het_threshold)
    dev.off()
  }
}

cat(sprintf("%d read x guide rows, %d errors -> %s\n", nrow(res$summary), length(res$errors), opt$outdir))
print(table(res$summary$call), quote=FALSE)
