# histograms of eve similarities

source("colors.R")
attach("../Analysis/R/tissue_text.RData")
attach("../Analysis/R/Rcache/expr_corr_betw_tissues.RData")
postscript("../SuppFigs/figS5.eps", height=7.9, width=6.5, pointsize=14,
           onefile=FALSE, horizontal=FALSE)
library(lineup)
par(mfrow=c(6,2), mar=c(4.1,1.1,1.1,1.1), oma=c(0, 1.5, 1.5, 0))

rng <- range(sapply(dcs, range, na.rm=TRUE))
for(i in tissues) {
  y <- pulldiag(dcs[[i]])
  hist(y, breaks=seq(-0.5, 1, by=0.02),
       ylab="", main="", xlab="Similarity",
       mgp=c(2.1, 0, 0), axes=FALSE)
  mtext(side=2, i, line=1, col=maincolor, font=2)
  if(i=="adipose")
    mtext(side=3, "Self-self similarity", line=1, col=maincolor, font=2)
  axis(side=1, seq(-0.5, 1, by=0.5), mgp=c(1.7, 0.7, 0.5), line=0.5)
  axis(side=1, seq(-0.25, 0.75, by=0.5), labels=rep("", 3), line=0.5)
  rug(y[y < 0.8], col=sexcolor[1], lwd=1, ticksize=-0.12,
      lend=1, ljoin=1)
  
  rn <- rownames(dcs[[i]])
  temp <- dcs[[i]]
  for(j in rn) 
    if(j %in% colnames(temp)) temp[j,j] <- NA
  hist(temp[temp >= -0.5], breaks=seq(-0.5, 1, by=0.02),
       ylab="", main="", xlab="Similarity",
       mgp=c(2.1, 0, 0), axes=FALSE)
  axis(side=1, seq(-0.5, 1, by=0.5), mgp=c(1.7, 0.7, 0.5), line=0.5)
  axis(side=1, seq(-0.25, 0.75, by=0.5), labels=rep("", 3), line=0.5)
  if(i=="adipose")
    mtext(side=3, "Self-nonself similarity", line=1, col=maincolor, font=2)
  rug(temp[temp > 0.8], col=sexcolor[1], lwd=1, ticksize=-0.12,
      lend=1, ljoin=1)
}

       
dev.off()
