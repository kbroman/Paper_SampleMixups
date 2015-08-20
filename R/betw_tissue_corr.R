# figure with densities of between-tissue correlations

library(broman) # for grayplot
source("colors.R")
attach("../Analysis/R/tissue_text.RData")
attach("../Analysis/R/Rcache/expr_corr.RData") # expr.corr, 40572 x 15

corr.dens <- vector("list", 15)
names(corr.dens) <- colnames(expr.corr)
for(i in 1:ncol(expr.corr))
  corr.dens[[i]] <- density(expr.corr[,i], from=-1, to=+1, n=2001)

postscript("../SuppFigs/figS4.eps", height=8.4, width=6.5, pointsize=12,
           onefile=FALSE, horizontal=FALSE)
par(mfcol=c(3,2), mar=c(4.1,2.1,3.1,2.1))
ymx <- max(sapply(corr.dens, function(a) max(a$y)))*1.05
set.seed(64163398)
subCCcolor <- CCcolor[-c(2,4)]
for(i in 1:6) {
  wh <- sample(which(tissuepairs[,1]==tissues[i] | tissuepairs[,2]==tissues[i]))
  short <- tissuepairs[wh,"short"]
  other <- short
  for(j in seq(along=wh)) {
    tmp <- unlist(tissuepairs[wh[j],1:2])
    other[j] <- tmp[tmp != tissues[i]]
  }
  thecolors <- subCCcolor[match(other, tissues)]

  grayplot(corr.dens[[wh[1]]]$x, corr.dens[[wh[1]]]$y, type="l", col=thecolors[1],
           vlines=seq(-1, 1, by=0.25), bgcolor=gray, ylim=c(0, ymx),
           yaxs="i", xaxs="i", hlines=NULL, yat=NA, xlab="Inter-tissue correlation",
           lwd=2, mgp=c(2.1, 0.5, 0))
  mtext(side=3, tissues[i], col=maincolor, line=0.5, font=2)
  for(j in seq(along=wh)[-1])
    lines(corr.dens[[wh[j]]]$x, corr.dens[[wh[j]]]$y, lwd=2, col=thecolors[j])

  if(i==1) legend("topleft", lwd=2, col=subCCcolor, tissues[1:6], bg=gray)
}
dev.off()
