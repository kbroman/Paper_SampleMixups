# supp figure with e vs e similarities (median correlations): 2nd highest vs maximum

source("colors.R")
attach("../Analysis/R/tissue_text.RData")
attach("../Analysis/R/Rcache/expr_mixup_summaries.RData")

postscript("../SuppFigs/figS6.eps", height=8.1, width=6.5, pointsize=14,
           onefile=FALSE, horizontal=FALSE)
par(mfcol=c(3,2), mar=c(4.1,4.1,3.1,2.1))
for(i in 1:6) {
  plot(mx[[i]], sec[[i]], xlim=c(0.55, 1), ylim=c(0.42, 1), xaxs="i", yaxs="i",
       xlab="", ylab="", xaxt="n", yaxt="n", type="n")
  mtext(side=3, tissues[i], col=maincolor, font=2, line=0.5)
  plim <- par("usr")
  rect(plim[1], plim[3], plim[2], plim[4], col=gray)

  abline(0,1, col=darkgray, lwd=3)

  xpos <- seq(0.6, 1, by=0.1)
  ypos <- seq(0.5, 1, by=0.1)
  abline(v=xpos, col="white")
  axis(side=1, at=xpos, mgp=c(3, 0.2, 0), tick=FALSE)
  abline(h=ypos, col="white")
  axis(side=2, at=ypos, mgp=c(3, 0.35, 0), tick=FALSE, las=1)
  title(xlab="Maximum similarity", mgp=c(1.9, 0, 0))
  title(ylab="Second highest similarity", mgp=c(2.5, 0, 0))


  # good, to fix, odd point
  grp <- match(as.numeric(cut(mx[[i]] - self[[i]], c(-Inf, 0, 0.1, Inf))), c(1,3,2))
  col <- distcolor[grp]

  points(mx[[i]][grp==1], sec[[i]][grp==1], bg=col[grp==1], pch=21)
  points(mx[[i]][grp!=1], sec[[i]][grp!=1], bg=col[grp!=1], pch=21)

  if(any(grp==3)) {
    wh <- which(grp==3)
    mult <- ifelse(i==5, -1, 1)
    xadj <- ifelse(i==5, 1, 0)
    text(mx[[i]][wh]+0.01*mult, self[[i]][wh], sub("Mouse", "", names(mx[[i]])[wh]),
         cex=1, adj=c(xadj, 0.5))
  }


  rect(plim[1], plim[3], plim[2], plim[4], col=NA)
}
dev.off()
