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
       main=paste0(ans$rec$name," — cut-site chromatogram [",ans$summary$status[1L],"]"))
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
  mtext("Guide blue | PAM purple | cut gold | elevated secondary peaks yellow | QC FAIL is NOT an edit call",
        1, line=3, cex=.62)
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

# Complementary evidence figures; each occupies its own plot/page.
# Cut offsets always follow the original sequencing direction.
plot_peak_ratio <- function(ans, het_threshold=0.30, flank=48L) {
  m <- ans$metrics
  if (is.null(m) || !nrow(m)) {
    graphics::plot.new(); graphics::text(.5,.5,"No peak metrics available")
    return(invisible(NULL))
  }
  cut <- ans$cut_raw
  if (length(cut)==1L && is.finite(cut) && cut >= 1L && cut <= nrow(m)) {
    lo <- max(1L, cut-flank); hi <- min(nrow(m),cut+flank)
  } else {
    lo <- max(1L,ans$hq$start)
    hi <- min(nrow(m), max(lo, min(ans$hq$end,lo+2L*flank)))
  }
  idx <- seq.int(lo,hi)
  rat <- m$second_to_primary[idx]
  finite <- rat[is.finite(rat)]
  ylim <- c(0, min(2, max(c(0.6, het_threshold*1.3, finite),na.rm=TRUE)*1.12))
  graphics::par(mar=c(5,4.7,3.2,1.5))
  graphics::plot(idx, rat, type="n", xlim=c(lo,hi), ylim=ylim,
                 xlab="AB1 base index (original sequencing direction)",
                 ylab="Secondary / primary peak height",
                 main=paste0("Cut-anchored mixed-peak evidence [",ans$summary$status[1L],"]"))
  resolved <- !is.null(ans$cut_boundary) && isTRUE(ans$cut_boundary$valid)
  if (resolved) {
    p <- ans$cut_boundary$raw_pre
    # Render separate scientifically defined windows; overlaps are expected.
    graphics::rect(p-35,0,p-9,ylim[2],col=grDevices::adjustcolor("#74a9cf",alpha.f=.09),border=NA)
    graphics::rect(p-5,0,p+8,ylim[2],col=grDevices::adjustcolor("#e5ae4c",alpha.f=.15),border=NA)
    graphics::rect(p+3,0,p+28,ylim[2],col=grDevices::adjustcolor("#6eb0a6",alpha.f=.10),border=NA)
  }
  graphics::abline(h=het_threshold, col="#be6a15", lty=2, lwd=1.5)
  if (length(idx)>1L) graphics::lines(idx,rat,col="#606d7b",lwd=1)
  export <- per_base_export(ans)
  eligible <- export$eligible_for_peak_scoring[idx]
  ok <- eligible & is.finite(rat)
  if (any(ok)) graphics::points(idx[ok],rat[ok],pch=16,cex=.66,col="#1f598c")
  if (any(!ok & is.finite(rat))) graphics::points(idx[!ok & is.finite(rat)],
                                                rat[!ok & is.finite(rat)],pch=1,cex=.55,col="#a4a4a4")
  hit <- ok & rat >= het_threshold
  if (any(hit)) graphics::points(idx[hit],rat[hit],pch=16,cex=.86,col="#cb5b19")
  if (length(cut)==1L && is.finite(cut)) {
    graphics::abline(v=cut+0.5,col="#976400",lwd=2)
    graphics::mtext(if (resolved) "Cut directly mapped" else "Cut approximate - manual review only",
                    side=3,line=.3,cex=.73,col="#815e0b")
  } else {
    graphics::mtext("Cut not mapped; no cut-localized genotype inference",
                    side=3,line=.3,cex=.74,col="#aa5c18")
  }
  graphics::legend("topright",inset=.01,
                   legend=c("Scorable mapped base", "Excluded / unmapped", "Above peak threshold", "Threshold"),
                   col=c("#1f598c","#a4a4a4","#cb5b19","#be6a15"),
                   pch=c(16,1,16,NA),lty=c(NA,NA,NA,2),bty="n",cex=.74)
  graphics::mtext("Sequencing direction: baseline -35..-9 | onset -5..+8 | post +3..+28 | QC failure overrides apparent pattern",
                  side=1,line=3.6,cex=.62)
  invisible(NULL)
}

plot_quality_profile <- function(ans) {
  m <- ans$metrics
  if (is.null(m) || !nrow(m)) {
    graphics::plot.new(); graphics::text(.5,.5,"No quality scores available")
    return(invisible(NULL))
  }
  idx <- m$read_index
  yy <- m$qual
  upper <- max(45,yy[is.finite(yy)],na.rm=TRUE) + 5
  graphics::par(mar=c(4.4,4.6,3.0,1.5))
  graphics::plot(range(idx),c(0,upper),type="n",xlab="AB1 base index (sequencing direction)",
                 ylab="Base-call Phred quality",
                 main=paste0("Read quality and usable alignment [",ans$summary$status[1L],"]"))
  if (ans$hq$len>0L) graphics::rect(ans$hq$start,0,ans$hq$end,upper,
                                  col=grDevices::adjustcolor("#8ac5aa",alpha.f=.17),border=NA)
  graphics::abline(h=20,col="#ad720a",lty=3)
  graphics::lines(idx,yy,col="#5e6c7b",lwd=1)
  ex <- per_base_export(ans)
  mapped <- ex$directly_mapped & is.finite(yy)
  if (any(mapped)) graphics::points(idx[mapped],yy[mapped],pch=16,cex=.55,col="#236da0")
  resolved <- !is.null(ans$cut_boundary) && isTRUE(ans$cut_boundary$valid)
  if (length(ans$cut_raw)==1L && is.finite(ans$cut_raw)) {
    graphics::abline(v=ans$cut_raw+0.5,col="#976400",lwd=2)
  }
  graphics::legend("topright",legend=c("Raw Phred", "Mapped HQ bases", "Q20", "HQ window", "Expected cut"),
                   col=c("#5e6c7b","#236da0","#ad720a","#78ad97","#976400"),
                   pch=c(NA,16,NA,15,NA),lty=c(1,NA,3,NA,1),bty="n",cex=.80)
  graphics::mtext(if (resolved) "Cut flanks directly resolved" else
                    "Cut unresolved or approximate - do not infer genotype from this marker",
                  side=3,line=.3,cex=.73,col=if (resolved) "#42634d" else "#a56616")
  invisible(NULL)
}

plot_batch_quality <- function(summary) {
  if (is.null(summary) || !nrow(summary)) {
    graphics::plot.new(); graphics::text(.5,.5,"No analyzed samples")
    return(invisible(NULL))
  }
  x <- as.numeric(summary$hq_len)
  y <- as.numeric(summary$mean_q)
  valid <- is.finite(x) & is.finite(y)
  if (!any(valid)) {
    graphics::plot.new(); graphics::text(.5,.5,"No usable HQ quality values")
    return(invisible(NULL))
  }
  graphics::par(mar=c(4.8,4.7,3.0,1.5))
  xr <- range(c(0,x[valid])); yr <- range(c(0,20,y[valid]))
  graphics::plot(x[valid],y[valid],type="n",xlim=c(xr[1],xr[2]+max(5,diff(xr)*.08)),
                 ylim=c(yr[1],yr[2]+max(3,diff(yr)*.12)),
                 xlab="High-quality window length (bp)", ylab="Mean HQ Phred",
                 main="Batch QC: quality versus usable sequence")
  graphics::abline(h=20,col="#a86d18",lty=3)
  cols <- ifelse(summary$status[valid]=="PASS","#25806d","#bc7732")
  graphics::points(x[valid],y[valid],pch=19,col=cols,cex=1.15)
  if (sum(valid)<=16L) {
    lab <- substr(as.character(summary$sample[valid]),1L,17L)
    graphics::text(x[valid],y[valid],labels=lab,pos=3,cex=.65,offset=.5,col="#40505c")
  }
  graphics::legend("bottomright",legend=c("QC/screen PASS", "QC/screen FAIL", "Q20 guide"),
                   col=c("#25806d","#bc7732","#a86d18"),pch=c(19,19,NA),
                   lty=c(NA,NA,3),bty="n",cex=.8)
  graphics::mtext("PASS does not establish a confirmed edit or WT genotype",side=1,line=3.7,cex=.72)
  invisible(NULL)
}
