# Browser-local, base-R reporting helpers. No uploads or results leave the browser
# when served by Shinylive. Download handlers write into the browser's webR FS.

export_stem <- function(x, max_chars=70L) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", as.character(x)[1L])
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) x <- "sample"
  substr(x, 1L, max_chars)
}

sample_export_name <- function(ans, extension, prefix="") {
  s <- export_stem(ans$summary$sample[1L])
  guide <- export_stem(ans$summary$guide[1L], max_chars=24L)
  stem <- if (nzchar(prefix)) paste(prefix, s, guide, sep="_") else paste(s, guide, sep="_")
  paste0(stem, ".", extension)
}

# Per-base results are preserved in ORIGINAL AB1 sequencing direction.
# The relative cut coordinate is descriptive, not a substitute for cut QC.
per_base_export <- function(ans, min_scoring_q=10L, min_primary=50) {
  m <- ans$metrics
  if (is.null(m) || !nrow(m)) return(data.frame())
  out <- m
  out$in_hq_window <- m$read_index >= ans$hq$start & m$read_index <= ans$hq$end
  out$reference_position <- rep(NA_integer_, nrow(m))
  if (!is.null(ans$aln) && nrow(ans$aln$map)) {
    mm <- ans$aln$map
    paired <- is.finite(mm[,"ref"]) & is.finite(mm[,"query"])
    if (any(paired)) {
      rawi <- as.integer(query_to_raw(as.integer(mm[paired,"query"]),
                                     ans$aln, ans$hq$start, ans$hq$len))
      good <- rawi >= 1L & rawi <= nrow(out)
      out$reference_position[rawi[good]] <- as.integer(mm[paired,"ref"])[good]
    }
  }
  out$directly_mapped <- is.finite(out$reference_position)
  resolved <- !is.null(ans$cut_boundary) && isTRUE(ans$cut_boundary$valid)
  out$read_offset_from_cut <- if (resolved) m$read_index - ans$cut_boundary$raw_pre else
    rep(NA_integer_, nrow(m))
  out$eligible_for_peak_scoring <- out$directly_mapped &
    is.finite(m$qual) & m$qual >= min_scoring_q &
    is.finite(m$primary) & m$primary >= min_primary &
    is.finite(m$second_to_primary) & m$base %in% BASES
  out
}

export_mismatch_table <- function(ans) {
  if (is.null(ans$mismatches) || !nrow(ans$mismatches)) {
    return(data.frame(ref_pos=integer(), ref_base=character(), read_base=character(),
                      raw_index=integer(), qual=numeric()))
  }
  ans$mismatches
}

export_gap_table <- function(ans) {
  if (is.null(ans$gaps) || !nrow(ans$gaps)) {
    return(data.frame(kind=character(), length=integer(), anchor=integer()))
  }
  ans$gaps
}

write_sample_evidence <- function(ans, path) {
  utils::write.csv(ans$summary, file=path, row.names=FALSE, na="")
  invisible(path)
}

write_sample_per_base <- function(ans, path) {
  utils::write.table(per_base_export(ans), file=path, row.names=FALSE,
                     sep="\t", quote=TRUE, na="")
  invisible(path)
}

# Parameter manifest makes a downloaded table reproducible against its figures.
write_settings_csv <- function(settings, path) {
  if (is.null(settings) || !length(settings)) {
    utils::write.csv(data.frame(parameter=character(),value=character()),
                     path,row.names=FALSE)
  } else {
    out <- data.frame(parameter=names(settings),
                      value=vapply(settings,as.character,character(1)),
                      stringsAsFactors=FALSE)
    utils::write.csv(out,path,row.names=FALSE)
  }
  invisible(path)
}

# All figures are derived from the analyzed sample, without silently changing
# its genotype decision. This report is a *screening record*, not confirmation.
write_sample_report_pdf <- function(ans, path, het_threshold=0.30, settings=NULL) {
  grDevices::pdf(file=path, width=11, height=7.5, onefile=TRUE, paper="special")
  on.exit(grDevices::dev.off(), add=TRUE)
  x <- ans$summary[1L, ]
  fmt <- function(v) if (length(v)==0L || is.na(v[1L])) "NA" else as.character(v[1L])
  eq <- function(nm, val) paste0(nm, ": ", fmt(val))
  lines <- c(
    "CRISPR Sanger Mapper | Sample review (v0.2.2)",
    "Research screening / QC only - genotype not independently confirmed",
    "",
    eq("Sample", x$sample), eq("Reference", x$reference), eq("Guide", x$guide),
    eq("Status", x$status), eq("Automated call", x$call),
    paste(eq("Read strand", x$strand), eq("Guide strand", x$guide_strand), sep="     "),
    paste(eq("PAM", x$pam), eq("PAM positions", paste0(fmt(x$pam_start),"-",fmt(x$pam_end))), sep="     "),
    paste(eq("Cut after reference base", x$cut_after), eq("Mapping", x$cut_mapping), sep="     "),
    paste(eq("HQ bases", x$hq_len), eq("HQ mean Phred", x$mean_q), sep="     "),
    paste(eq("Aligned bases", x$aligned_bases), eq("Identity", x$identity), sep="     "),
    paste(eq("Baseline n", x$pre_scored_n), eq("Onset mixed n", x$onset_mixed_n),
          eq("Post-cut scored n", x$post_scored_n), sep="     "),
    paste(eq("Baseline median peak ratio", x$pre_median_ratio),
          eq("Post-cut median peak ratio", x$post_median_ratio), sep="     "),
    eq("Cut-localized pattern", x$cut_localized_pattern),
    "", "Manual-review notes:"
  )
  manual <- fmt(x$manual_review)
  if (!nzchar(manual)) manual <- "None recorded; compare against controls before biological interpretation."
  lines <- c(lines, strwrap(manual, width=100L),
             "", "Settings snapshot:")
  if (!is.null(settings) && length(settings)) {
    lines <- c(lines, vapply(names(settings), function(k) paste(k, "=", as.character(settings[[k]])), ""))
  } else {
    lines <- c(lines, "No settings snapshot supplied")
  }
  lines <- c(lines,
             "", "Interpretation cautions:",
             "- A valid PAM and an aligned cut are prerequisites for a cut-associated inference.",
             "- Unresolved/near-cut mixture is not an exact allele or allele fraction.",
             "- WT-like indicates no qualifying evidence under current QC thresholds, not proof of WT.",
             "- Compare strain-matched WT and validate interesting samples independently.")
  graphics::plot.new()
  graphics::par(usr=c(0, 1, 0, 1))
  y <- seq(0.95, 0.06, length.out=length(lines))
  graphics::text(0.04, y, labels=lines, adj=c(0, .5),
                 cex=min(0.85, max(0.48, 32/length(lines))),
                 font=ifelse(seq_along(lines)==1, 2, 1),
                 col=ifelse(seq_along(lines)==2, "#aa5c18", "#222222"))
  plot_trace(ans, het_threshold=het_threshold)
  plot_peak_ratio(ans, het_threshold=het_threshold)
  plot_quality_profile(ans)
  invisible(path)
}


save_figure_png <- function(path, draw, width=1700L, height=860L, res=140L) {
  grDevices::png(filename=path, width=width, height=height, res=res, bg="white")
  on.exit(grDevices::dev.off(), add=TRUE)
  draw()
  invisible(path)
}
