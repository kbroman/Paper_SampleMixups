######################################################################
# eQTL genotype vs expression: main figure
######################################################################

source("colors.R")
library(lineup)


attach("../Analysis/R/Rcache/dgve.RData")
y <- attr(dgve[["islet"]], "y")
obsg <- attr(dgve[["islet"]], "obsg")
infg <- attr(dgve[["islet"]], "infg")

attach("../Analysis/R/Rcache/pmark.RData")

gene <- "499541"
wh <- which(sapply(y, function(a,b) any(colnames(a) == b),gene))
e <- y[[wh]]
og <- obsg[,wh]
ig <- infg[,wh]

id <- findCommonID(rownames(e), names(og))
e <- e[id$first,]
ig <- ig[id$first]
ig[is.na(ig)] <- 0
og <- og[id$second]

set.seed(35410287)
u <- rnorm(length(og), 0, 0.08)

pmar <- pmark[gene, 2]

gnames <- attr(dgve[["islet"]], "genonames")

postscript("../Figs/fig5.eps", height=4, width=4, onefile=FALSE, horizontal=FALSE, pointsize=10)
par(mar=c(3.1,4.1,1.1,1.1))
plot((og+u), e, type="n", ylim=c(-1.5, 0.5), xlim=c(0.5,3.5), xaxs="i",
     xaxt="n", xlab="", yaxt="n", ylab="")
title(xlab=paste("Genotype at", pmar), mgp=c(1.9,0.5,0))
title(ylab=paste("expression of probe", gene), mgp=c(3, 0.5, 0))
plim <- par("usr")
rect(plim[1], plim[3], plim[2], plim[4], col=gray)
abline(v=1:3, col=darkgray, lwd=2)
abline(h=seq(-1.5, 0.5, by=0.5), col="white")
points((og+u)[ig != 0], e[ig != 0], bg=f2color[ig[ig != 0]], pch=21)
points((og+u)[ig == 0], e[ig == 0], bg=CCcolor[3], pch=21)
axis(side=1, at=1:3, labels=gnames, tick=FALSE, mgp=c(3, 0.1, 0))
axis(side=2, at=seq(-1.5, 0.5, by=0.5), mgp=c(3, 0.3, 0), las=1, tick=FALSE)
rect(plim[1], plim[3], plim[2], plim[4], col=NA, border="black", lend=1, ljoin=1)
dev.off()
