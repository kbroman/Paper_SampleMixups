# figure with g vs e similarity (overall proportion matches), for each tissue

rm(list=ls())
source("colors.R")
attach("../Analysis/R/tissue_text.RData")
attach("../Analysis/R/Rcache/dgve_min_and_self.RData")

# good, fixable, not found
grp <- rep(3, length(mn[[7]]))
grp[!is.na(self[[7]]) & mn[[7]] == self[[7]]] <- 1
grp[(is.na(self[[7]]) | mn[[7]] < self[[7]]) & mn[[7]] < 0.3] <- 2

postscript("../SuppFigs/figS16.eps", height=8.1, width=6.5, pointsize=14,
           onefile=FALSE, horizontal=FALSE)

hi1 <- 1.02
lo1 <- 1-0.87-0.015
hi2 <- lo1 - (hi1-lo1)*0.04
lo2 <- hi2 - (hi1-lo1)*0.10

par(mfcol=c(3,2), mar=c(3.3, 4.1, 2.6, 1.1))
for(i in 1:6) {
  plot(1-mn[[i]], 1-self[[i]], xlim=c(0.54, 1.01), ylim=c(lo2, hi1), xaxs="i", yaxs="i",
     xlab="", ylab="", xaxt="n", yaxt="n", type="n", bty="n")
  plim <- par("usr")
  rect(plim[1], lo1, plim[2], plim[4], col=gray)
  rect(plim[1], lo2, plim[2], hi2, col=gray)
  mtext(side=3, tissues[i], col=maincolor, line=0.5, font=2)

  abline(0,1, col=darkgray, lwd=3)

  xpos <- seq(0.6, 1, by=0.1)
  ypos <- seq(0.2, 1, by=0.2)
  abline(v=xpos, col="white")
  axis(side=1, at=xpos, mgp=c(3, 0.2, 0), tick=FALSE)
  abline(h=ypos, col="white")
  axis(side=2, at=ypos, mgp=c(3, 0.35, 0), tick=FALSE, las=1)
  axis(side=2, at=(lo2+hi2)/2, "N/A", mgp=c(3, 0.35, 0), tick=FALSE, las=1)
  title(xlab="Maximum similarity", mgp=c(1.9, 0, 0))
  title(ylab="Self similarity", mgp=c(2.5, 0, 0))

  library(beeswarm)
  x <- beeswarm(mn[[i]][is.na(self[[i]])], do.plot=FALSE, method="center")$x-1
  self[[i]][is.na(self[[i]])] <- x/diff(range(x))*(hi2-lo2)*0.6+(1-c(hi2+lo2)/2)

  col <- distcolor[grp]

  points(1-mn[[i]], 1-self[[i]], bg=col, pch=21)

  rect(plim[1], lo1, plim[2], plim[4], col=NA)
  rect(plim[1], lo2, plim[2], hi2, col=NA)
}

dev.off()
for(i in 1:3) detach(2)

