# example between-tissue expression scatterplots: with all probes and with a subset of correlated ones

library(qtl)
source("colors.R")
attach("../Analysis/R/Rcache/mlratios.RData")
attach("../Analysis/R/calls.RData")
attach("../Analysis/R/Rcache/expr_corr.RData")

#set.seed(50527372)
#ind <- sample(rownames(calls)[rownames(calls) == calls[,"kidney"] & rownames(calls)==calls[,"liver"]], 1)
ind <- "Mouse3567"

jpeg("../SuppFigs/figS3.jpg", height=650*2, width=650*2, pointsize=36, quality=100)
par(mar=c(3.1, 3.6, 0.6, 0.6))
plot(kidney.mlratio[,ind], liver.mlratio[,ind],
     xlab="", ylab="", xlim=c(-2.05, 2.05), ylim=c(-2.05,2.05),
     xaxt="n", yaxt="n", xaxs="i", yaxs="i",
     type="n")
plim <- par("usr")
rect(plim[1], plim[3], plim[2], plim[4], col=gray)

abline(0,1, col=darkgray, lwd=2)

pos <- seq(-2, 2, by=0.5)
axis(side=1, at=pos, mgp=c(3, 0.2, 0), tick=FALSE)
axis(side=2, at=pos, mgp=c(3, 0.35, 0), tick=FALSE, las=1)
abline(h=pos, v=pos, col="white", lwd=2)
title(xlab="Gene expression in kidney", mgp=c(1.8, 0, 0))
title(ylab="Gene expression in liver", mgp=c(2.5, 0, 0))

points(kidney.mlratio[,ind], liver.mlratio[,ind], cex=0.7, pch=16, col="gray50")
wh <- which(expr.corr[,"kl"] > 0.75)
points(kidney.mlratio[wh,ind], liver.mlratio[wh,ind], bg=sexcolor[1], pch=21, cex=0.9)

rect(plim[1], plim[3], plim[2], plim[4], col=NA)
dev.off()
