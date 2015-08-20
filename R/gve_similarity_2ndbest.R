# figure with g vs e similarity: Maximum vs 2nd best

rm(list=ls())

source("colors.R")
attach("../Analysis/R/tissue_text.RData")
attach("../Analysis/R/Rcache/dgve_min_and_self.RData")

postscript("../SuppFigs/figS17.eps", height=6.5, width=6.5, pointsize=14,
           onefile=FALSE, horizontal=FALSE)

par(mar=c(3.3, 4.1, 1.1, 1.1))
plot(1-mn[[7]], 1-sec[[7]], xlim=c(0.51, 1.01), ylim=c(0.51, 1.01), xaxs="i", yaxs="i",
     xlab="", ylab="", xaxt="n", yaxt="n", type="n")
plim <- par("usr")
rect(plim[1], plim[3], plim[2], plim[4], col=gray)

abline(0,1, col=darkgray, lwd=3)

xpos <- seq(0.6, 1, by=0.1)
ypos <- seq(0.6, 1, by=0.1)
abline(v=xpos, col="white")
axis(side=1, at=xpos, mgp=c(3, 0.2, 0), tick=FALSE)
abline(h=ypos, col="white")
axis(side=2, at=ypos, mgp=c(3, 0.35, 0), tick=FALSE, las=1)
title(xlab="Maximum similarity", mgp=c(1.9, 0, 0))
title(ylab="Second highest similarity", mgp=c(2.5, 0, 0))


# good, fixable, not found
grp <- rep(3, length(mn[[7]]))
grp[!is.na(self[[7]]) & mn[[7]] == self[[7]]] <- 1
grp[(is.na(self[[7]]) | mn[[7]] < self[[7]]) & mn[[7]] < 0.2] <- 2

col <- distcolor[grp]

points(1-mn[[7]], 1-sec[[7]], bg=col, pch=21)

distcolor.text <- distcolor
for(i in 2:3) {
  tmp <- col2rgb(distcolor[i])
  distcolor.text[i] <- rgb(tmp[1]*0.7, tmp[2]*0.7, tmp[3]*0.7, maxColorValue=255)
}
distcolor.text[1] <- "gray30"

text(1-0.15, 1-0.44, "Good/", col=distcolor.text[1], font=2)
text(1-0.15, 1-0.46, "Fixable", col=distcolor.text[2], font=2)
text(1-0.35, 1-0.468, "Not found", col=distcolor.text[3], font=2)


rect(plim[1], plim[3], plim[2], plim[4], col=NA)

dev.off()
for(i in 1:2) detach(2)
