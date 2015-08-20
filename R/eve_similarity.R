# figure with e vs e similarities (median correlations)

rm(list=ls())
source("colors.R")
attach("../Analysis/R/tissue_text.RData")
attach("../Analysis/R/Rcache/expr_mixup_summaries.RData")

postscript("../Figs/fig3.eps", height=8.1, width=6.5, pointsize=14,
           onefile=FALSE, horizontal=FALSE)
par(mfcol=c(3,2), mar=c(4.1,4.1,3.1,2.1))

yhi1 <- 1
ylo1 <- -0.29
yhi2 <- ylo1 - (yhi1-ylo1)*0.04
ylo2 <- yhi2 - (yhi1-ylo1)*0.10

for(i in 1:6) {
  plot(mx[[i]], self[[i]], xlim=c(0.5, 1), ylim=c(ylo2, yhi1), xaxs="i", yaxs="i",
       xlab="", ylab="", xaxt="n", yaxt="n", type="n", bty="n")
  mtext(side=3, tissues[i], col=maincolor, line=0.5, font=2)
  plim <- par("usr")
  rect(plim[1], ylo1, plim[2], plim[4], col=gray)
  rect(plim[1], ylo2, plim[2], yhi2, col=gray)

  abline(0,1, col=darkgray, lwd=3)

  xpos <- seq(0.5, 1, by=0.1)
  ypos <- seq(-0.2, 1, by=0.2)
  abline(v=xpos, col="white")
  axis(side=1, at=xpos, mgp=c(3, 0.2, 0), tick=FALSE)
  abline(h=ypos, col="white")
  axis(side=2, at=ypos, mgp=c(3, 0.35, 0), tick=FALSE, las=1)
  axis(side=2, at=(ylo2+yhi2)/2, "N/A", mgp=c(3, 0.35, 0), tick=FALSE, las=1)
  title(xlab="Maximum similarity", mgp=c(1.9, 0, 0))
  title(ylab="Self similarity", mgp=c(2.5, 0, 0))


  # good, to fix, odd point
  grp <- match(as.numeric(cut(mx[[i]] - self[[i]], c(-Inf, 0, 0.1, Inf))), c(1,3,2))
  grp[is.na(self[[i]])] <- 1
  col <- distcolor[grp]

  if(i==5) {  # for kidney, include the points with missing self values
    library(beeswarm)
    x <- beeswarm(mx[[5]][is.na(self[[5]])], do.plot=FALSE, method="center", scale=0.8)$x-1
    z <- mx[[5]][is.na(self[[5]])]
    self[[5]][is.na(self[[5]])] <- x/diff(range(x))*(yhi2-ylo2)*0.5+(yhi2+ylo2)/2
  }

  points(mx[[i]], self[[i]], bg=col, pch=21)

  if(any(grp==3)) {
    wh <- which(grp==3)
    text(mx[[i]][wh]+0.01, self[[i]][wh], sub("Mouse", "", names(mx[[i]])[wh]),
         cex=1, adj=c(0, 0.5))
  }

  wh <- which(grp==2)
  if(i==1) {
    wh <- wh[order(mx[[i]][wh])]
    xalign <- c(1, 1, 0.5, 0.5, 0)
    yalign <- c(0.5, 0.5, 1, 0, 0.5)
    xadj <- c(-0.01, -0.01, 0, 0, 0.01)
    yadj <- c(0, 0, -0.045, 0.045, 0)
    for(j in seq(along=wh))
      text(mx[[i]][wh[j]]+xadj[j], self[[i]][wh[j]]+yadj[j], sub("Mouse", "", names(mx[[i]])[wh[j]]),
           cex=1, adj=c(xalign[j], yalign[j]))
  }    

  if(i==2) {
    text(mx[[i]][wh], self[[i]][wh]+0.045, sub("Mouse", "", names(mx[[i]])[wh]),
         cex=1, adj=c(0.5, 0))
  }

  if(i==4) {
    wh <- wh[order(mx[[i]][wh])]
    xalign <- c(1, 1, 0)
    yalign <- rep(0.5, 3)
    xadj <- c(-0.01, -0.01, 0.01)
    yadj <- rep(0,3)
    for(j in seq(along=wh))
      text(mx[[i]][wh[j]]+xadj[j], self[[i]][wh[j]]+yadj[j], sub("Mouse", "", names(mx[[i]])[wh[j]]),
           cex=1, adj=c(xalign[j], yalign[j]))
  }    

  if(i==5) {
    wh <- wh[order(mx[[i]][wh])]
    xalign <- c(1, 0)
    yalign <- rep(0.5, 2)
    xadj <- c(-0.01, 0.01)
    yadj <- rep(0, 2)
    for(j in seq(along=wh))
      text(mx[[i]][wh[j]]+xadj[j], self[[i]][wh[j]]+yadj[j], sub("Mouse", "", names(mx[[i]])[wh[j]]),
           cex=1, adj=c(xalign[j], yalign[j]))
  }    

  if(i==6) {
    wh <- wh[order(mx[[i]][wh])]
    xalign <- c(1, 0)
    yalign <- rep(0.5, 2)
    xadj <- c(-0.01, 0.01)
    yadj <- rep(0, 2)
    for(j in seq(along=wh))
      text(mx[[i]][wh[j]]+xadj[j], self[[i]][wh[j]]+yadj[j], sub("Mouse", "", names(mx[[i]])[wh[j]]),
           cex=1, adj=c(xalign[j], yalign[j]))
  }    


  rect(plim[1], ylo1, plim[2], plim[4], col=NA)
  rect(plim[1], ylo2, plim[2], yhi2, col=NA)
}
dev.off()
for(i in 1:3) detach(2)
