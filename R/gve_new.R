######################################################################
# eQTL genotype vs expression with the revised data
######################################################################

source("colors.R")
library(lineup)

attach("../Analysis/FinalData/islet_mlratio_final.RData")
attach("../Analysis/FinalData/aligned_geno_with_pmap.RData")

library(qtl)
f2g$pheno$id <- f2g$pheno$MouseNum
id <- findCommonID(f2g, islet.mlratio)

postscript("../SuppFigs/figS18.eps", height=9.75, width=6.5,
           onefile=FALSE, horizontal=FALSE, pointsize=12)

layout(rbind(c(2,3), c(4,5), c(1,6)))

gene <- "499541"
y <- islet.mlratio[id$second, gene]
f2g <- fill.geno(f2g, err=0.002, map.function="c-f")
pmar <- "rs13476177"
g <- pull.geno(f2g[,id$first])[,pmar]
set.seed(35410287)
u <- rnorm(length(g), 0, 0.08)

gnames <- c("BB", "BR", "RR")

par(mar=c(4.1,4.1,2.1,1.1))
plot((g+u), y, type="n", ylim=c(-1.5, 0.5), xlim=c(0.5,3.5), xaxs="i",
     xaxt="n", xlab="", yaxt="n", ylab="")
title(xlab=paste("Genotype at", pmar), mgp=c(1.9,0.5,0))
title(ylab=paste("expression of probe", gene), mgp=c(3, 0.5, 0))
plim <- par("usr")
rect(plim[1], plim[3], plim[2], plim[4], col=gray)
abline(v=1:3, col=darkgray, lwd=2)
abline(h=seq(-1.5, 0.5, by=0.5), col="white")
points((g+u), y, bg=f2color[g], pch=21)
axis(side=1, at=1:3, labels=gnames, tick=FALSE, mgp=c(3, 0.1, 0))
axis(side=2, at=seq(-1.5, 0.5, by=0.5), mgp=c(3, 0.3, 0), las=1, tick=FALSE)
rect(plim[1], plim[3], plim[2], plim[4], col=NA, border="black", lend=1, ljoin=1)

text(plim[1]-diff(plim[1:2])*0.15, plim[4]+diff(plim[3:4])*0.08, LETTERS[5], font=2, xpd=TRUE, cex=1.3)


par(mar=c(4.1,4.1,2.1,1.1))
probes <- list(c("502129", "517583"),
               c("10002898467", "510446"),
               c("592144", "10003116152"),
               c("10004035256", "513432"))
markers <- c("rs13478396", "rs13478978",
             "rs13480320", "rs3023193")


selected <- c(62, 3, 11, 22)

for(i in 1:4) {
  cat("****", i, "****\n")

  e <- islet.mlratio[id$second, probes[[i]]]
  g <- pull.geno(f2g[,id$first])[,markers[i]]

  pmar <- markers[i]

  plot(e[,1], e[,2], type="n", xlab="", ylab="", xaxt="n", yaxt="n")
  title(ylab=paste("expression of probe", colnames(e)[2]), mgp=c(2.2, 0.5, 0))
  title(xlab=paste("expression of probe", colnames(e)[1]), mgp=c(1.5, 0.5, 0))

  plim <- par("usr")
  rect(plim[1], plim[3], plim[2], plim[4], col=gray)
  xl <- pretty(e[,1])
  if(i==3) { # trying to force it!
    xl <- c(-0.3, -0.1, 0.1, 0.3)
    abline(v=xl, col="white")
    axis(side=1, at=xl, mgp=c(3, 0.1, 0), tick=FALSE)
    xl <- c(-0.2, 0, 0.2)
  }
  abline(v=xl, col="white")
  axis(side=1, at=xl, mgp=c(3, 0.1, 0), tick=FALSE)

  yl <- pretty(e[,2])
  abline(h=yl, col="white")
  axis(side=2, at=yl, mgp=c(3, 0.3, 0), las=1, tick=FALSE)

  text(plim[1]-diff(plim[1:2])*0.15, plim[4]+diff(plim[3:4])*0.08, LETTERS[i], font=2, xpd=TRUE, cex=1.3)


  points(e[,1], e[,2], bg=f2color[g], pch=21)

  location <- "bottomright"
  if(i==1) location <- "topleft"
  if(i==4) location <- "bottomleft"

  if(i==1) {
    legend(location, pch=21, pt.bg=f2color, title=pmar,
           legend=gnames, bg="white")
  } else {
    pmar.w <- strwidth(pmar)*1.1
    pmar.h <- strheight(pmar)*2
    if(location=="bottomright") {
      rect(plim[2]-pmar.w, plim[3], plim[2], plim[3]+pmar.h, col="white")
      text(plim[2]-pmar.w/2, plim[3]+pmar.h/2, pmar)
    } else {
      rect(plim[1]+pmar.w, plim[3], plim[1], plim[3]+pmar.h, col="white")
      text(plim[1]+pmar.w/2, plim[3]+pmar.h/2, pmar)
    }
  }

  rect(plim[1], plim[3], plim[2], plim[4], col=NA, border="black", lend=1, ljoin=1)
}
for(i in 1:2) detach(2)
dev.off()
