base_col <- c(A="#2ca02c", C="#1f77b4", G="#111111", T="#d62728")

plot_trace <- function(ans, flank=38L) {
  rec <- ans$rec; m <- ans$metrics
  if (is.null(rec) || !nrow(m)) { plot.new(); text(.5,.5,"No trace"); return(invisible()) }
  if (is.finite(ans$cut_raw)) {
    lo <- max(1L, ans$cut_raw-flank)
    hi <- min(nrow(m), ans$cut_raw+flank)
  } else {
    lo <- if (!is.null(ans$hq)) ans$hq$start else 1L
    hi <- min(nrow(m), lo+2L*flank)
  }
  if (hi < lo) { lo <- 1L; hi <- min(nrow(m),70L) }
  # ABIF PLOC uses 0-based trace positions; R vectors are 1-based.
  x0 <- max(1L, as.integer(rec$ploc[lo])+1L)
  x1 <- min(rec$trace_len, as.integer(rec$ploc[hi])+1L)
  if (x1 < x0) { plot.new(); text(.5,.5,"Trace peak indices invalid"); return(invisible()) }
  xs <- seq.int(x0,x1)
  ymax <- max(vapply(BASES, function(b) max(rec$channels[[b]][xs],na.rm=TRUE), numeric(1)), na.rm=TRUE)
  if (!is.finite(ymax) || ymax<=0) ymax <- 1
  par(mar=c(4,2.3,3.1,1), xaxs="i", yaxs="i")
  plot(range(xs), c(-0.15*ymax,1.13*ymax), type="n", xlab="Raw sequencing trace position", ylab="", yaxt="n",
       main=paste(ans$rec$name,"— sgRNA/PAM/cut-aligned trace"))

  # Shade precisely mapped protospacer/PAM sites, not an estimated fixed offset.
  gh <- ans$guide_hit
  aln <- ans$aln
  hq <- ans$hq
  if (!is.null(gh) && !is.null(aln) && !is.null(hq)) {
    shade_ref <- function(start,end,fill) {
      b1 <- map_ref_to_raw(start,aln,hq$start,hq$len)
      b2 <- map_ref_to_raw(end,aln,hq$start,hq$len)
      if (is.finite(b1) && is.finite(b2)) {
        xlo <- as.integer(rec$ploc[min(b1,b2)])+1L
        xhi <- as.integer(rec$ploc[max(b1,b2)])+1L
        if (xhi >= x0 && xlo <= x1)
          rect(max(x0,xlo),0,min(x1,xhi),ymax,
               col=adjustcolor(fill, alpha.f=.14),border=NA)
      }
    }
    shade_ref(gh$guide_start, gh$guide_end,"#4a86bd")
    shade_ref(gh$pam_start, gh$pam_end,"#e28f46")
  }
  if (is.finite(ans$cut_raw)) {
    for (i in seq.int(lo,hi)) {
      p <- as.integer(rec$ploc[i])+1L
      if (is.finite(m$second_to_primary[i]) && m$second_to_primary[i] >= 0.30) {
        rect(p-2,0,p+2,ymax,
             col=adjustcolor("#f3c742",alpha.f=.12),border=NA)
      }
    }
  }
  for (b in BASES) lines(xs,rec$channels[[b]][xs],col=base_col[b],lwd=1.1)
  for (i in seq.int(lo,hi)) {
    p <- as.integer(rec$ploc[i])+1L
    cc <- unname(base_col[m$base[i]])
    text(p,-0.06*ymax,m$base[i],col=if (length(cc) && !is.na(cc)) cc else "#777777",cex=.68,font=2)
  }
  if (!is.null(ans$cut_boundary) && isTRUE(ans$cut_boundary$valid)) {
    cb <- ans$cut_boundary
    cut_x <- mean(as.integer(rec$ploc[c(cb$raw_pre,cb$raw_post)])+1L)
    abline(v=cut_x,col="#af253a",lwd=2.5)
    text(cut_x,1.07*ymax,paste0("cut ",cb$ref_left,"|",cb$ref_right),
         col="#af253a",cex=.8,font=2,adj=c(.5,.5))
  }
  legend("topright",legend=c(BASES,"sgRNA","PAM","predicted cut"),
         col=c(base_col,"#4a86bd","#e28f46","#af253a"),
         lty=1,lwd=c(rep(2,4),7,7,2),bty="n",horiz=FALSE,cex=.73,ncol=2)
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
