base_col <- c(A="#2ca02c", C="#1f77b4", G="#111111", T="#d62728")

plot_trace <- function(ans, flank=35L, het_threshold=0.30) {
  rec <- ans$rec; m <- ans$metrics
  if (is.null(rec) || !nrow(m)) { plot.new(); text(.5,.5,"No trace"); return(invisible()) }
  if (is.finite(ans$cut_raw)) {
    lo <- max(1L, ans$cut_raw-flank); hi <- min(nrow(m), ans$cut_raw+flank)
  } else { lo <- ans$hq$start; hi <- min(ans$hq$end, lo+2L*flank) }
  if (hi < lo) { lo <- 1L; hi <- min(nrow(m),70L) }
  # PLOC is zero based.
  x0 <- max(1L, as.integer(rec$ploc[lo])+1L)
  x1 <- min(rec$trace_len, as.integer(rec$ploc[hi])+1L)
  xs <- x0:x1
  ymax <- max(vapply(BASES, function(b) max(rec$channels[[b]][xs],na.rm=TRUE), numeric(1)), na.rm=TRUE)
  par(mar=c(4.4,2.3,2.5,1), mgp=c(2.2,0.7,0), xaxs="i", yaxs="i")
  plot(range(xs), c(-0.13*ymax,1.10*ymax), type="n", xlab="Trace position", ylab="", yaxt="n",
       main=paste(ans$rec$name,"— chromatogram around target"))
  for (b in BASES) lines(xs, rec$channels[[b]][xs], col=base_col[b], lwd=1.2)
  for (i in lo:hi) {
    p <- as.integer(rec$ploc[i])+1L
    text(p, -0.055*ymax, m$base[i], col=base_col[m$base[i]] %||% "#666666", cex=.72, font=2)
    if (is.finite(m$second_to_primary[i]) && m$second_to_primary[i]>=het_threshold) rect(p-3,0,p+3,ymax, col=adjustcolor("#f3c742",alpha.f=.14), border=NA)
  }
  # Shade validated protospacer and PAM, only where directly mapped into this
  # chromatogram; the raw coordinates are always in original sequencing order.
  if (!is.null(ans$guide_hit) && !is.null(ans$aln)) {
    gh <- ans$guide_hit
    for (item in list(list(a=gh$guide_start, b=gh$guide_end, col="#548ad4"),
                      list(a=gh$pam_start, b=gh$pam_end, col="#a576bd"))) {
      anchors <- vapply(seq.int(item$a,item$b), function(rp) {
        map_ref_to_raw(rp, ans$aln, ans$hq$start, ans$hq$len)
      }, integer(1))
      anchors <- anchors[is.finite(anchors) & anchors>=lo & anchors<=hi]
      if (length(anchors)) {
        aa <- as.integer(rec$ploc[min(anchors)]) + 1L
        bb <- as.integer(rec$ploc[max(anchors)]) + 1L
        rect(aa-4, 0, bb+4, ymax, col=adjustcolor(item$col, alpha.f=.055), border=NA)
      }
    }
  }
  if (is.finite(ans$cut_raw)) {
    # cut_raw is the final base before cut in sequencing direction, or an
    # approximate location if the two immediately flanking bases do not map.
    nxt <- min(length(rec$ploc), ans$cut_raw + 1L)
    cp <- (as.numeric(rec$ploc[ans$cut_raw]) + as.numeric(rec$ploc[nxt])) / 2 + 1
    abline(v=cp, col="#a16a00", lwd=2)
    ctype <- as.character(ans$summary$cut_mapping[1])
    label <- if (ctype == "aligned_two_flanks") "PAM-validated cut"
             else paste0("cut (", ctype, "; review)")
    text(cp, 1.04*ymax, label, col="#a16a00", cex=.72, font=2)
  }
  legend("topright", legend=BASES, col=base_col, lty=1, lwd=2, bty="n", horiz=TRUE, cex=.8)
  mtext("Guide: blue background | PAM: purple | cut: gold line | mixed peaks: yellow", 1, line=3, cex=.65)
}

alignment_strings <- function(ans, width=90L) {
  if (is.null(ans$aln)) return("No alignment available")
  ar <- ans$aln$aligned_ref; aq <- ans$aln$aligned_query
  if (!nzchar(ar)) return("No alignment available")
  mid <- paste(ifelse(strsplit(ar,"")[[1]]==strsplit(aq,"")[[1]],"|"," "),collapse="")
  starts <- seq(1,nchar(ar),by=width)
  chunks <- vapply(starts,function(s){
    e <- min(nchar(ar),s+width-1L)
    paste0("REF  ",substr(ar,s,e),"\n     ",substr(mid,s,e),"\nREAD ",substr(aq,s,e))
  },character(1))
  paste(chunks,collapse="\n\n")
}
