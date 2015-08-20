# histograms of gve similarities

source("colors.R")
attach("../Analysis/R/tissue_text.RData")
attach("../Analysis/R/Rcache/dgve.RData")
postscript("../SuppFigs/figS15.eps", height=7.9, width=6.5, pointsize=14,
           onefile=FALSE, horizontal=FALSE)
library(lineup)
par(mfrow=c(6,2), mar=c(4.1,1.1,1.1,1.1), oma=c(0, 1.5, 1.5, 0))

for(i in tissues) {
  y <- 1-pulldiag(dgve[[i]])
  hist(y, breaks=seq(0, 1, by=0.01),
       ylab="", yaxt="n", main="", xlab="Similarity",
       xaxt="n", mgp=c(2.1, 0, 0), bty="o", xaxs="i", axes=FALSE)
  mtext(side=2, i, line=1, col=maincolor, font=2)
  if(i=="adipose")
    mtext(side=3, "Self-self similarity", line=1, col=maincolor, font=2)
  axis(side=1, seq(0, 1, by=0.25), mgp=c(1.7, 0.7, 0.5), line=0.5)
  rug(y[y < 0.8], col=sexcolor[1], lwd=1, ticksize=-0.12,
      lend=1, ljoin=1)
  
  rn <- rownames(dgve[[i]])
  y <- 1-dgve[[i]]
  for(j in rn) 
    if(j %in% colnames(y)) y[j,j] <- NA
  hist(y, breaks=seq(0, 1, by=0.02),
       ylab="", yaxt="n", main="", xlab="Similarity",
       xaxt="n", mgp=c(2.1, 0, 0), bty="o", xaxs="i", axes=FALSE)
  axis(side=1, seq(0, 1, by=0.25), mgp=c(1.7, 0.7, 0.5), line=0.5)
  if(i=="adipose")
    mtext(side=3, "Self-nonself similarity", line=1, col=maincolor, font=2)
  rug(y[y > 0.8], col=sexcolor[1], lwd=1, ticksize=-0.12,
      lend=1, ljoin=1)
}
       
dev.off()
