# scatterplots for example of swap of expression arrays

source("colors.R")
attach("../Analysis/R/Rcache/mlratios.RData")
attach("../Analysis/R/Rcache/expr_corr.RData")
attach("../Analysis/R/tissue_text.RData")

tissue <- "gastroc"
tissue.index <- match(tissue, tissues)
ind <- c("Mouse3655", "Mouse3659")

postscript("../SuppFigs/figS7.eps", height=8.1, width=6.5, pointsize=12,
           onefile=FALSE, horizontal=FALSE)
par(mfcol=c(5,4), mar=c(3.5,3.6,0.6,1.1), oma=c(0, 0, 1.6, 0))
for(i in ind) {
  for(j in ind) {
    for(k in seq(along=tissues)[-tissue.index]) {
      inum <- sub("Mouse", "", i)
      jnum <- sub("Mouse", "", j)
      plot(0,0,
           xlab="", ylab="", xlim=c(-2.1, 2.1), ylim=c(-2.1,2.1),
           xaxt="n", yaxt="n", xaxs="i", yaxs="i",
           type="n")
      plim <- par("usr")
      rect(plim[1], plim[3], plim[2], plim[4], col=gray)

      abline(0,1, col=darkgray, lwd=2)

      pos <- seq(-2, 2, by=1)
      axis(side=1, at=pos, mgp=c(3, 0.2, 0), tick=FALSE)
      axis(side=2, at=pos, mgp=c(3, 0.35, 0), tick=FALSE, las=1)
      abline(h=pos, v=pos, col="white")
      title(xlab=paste(tissue, inum), mgp=c(1.6, 0, 0))
      title(ylab=paste(tissues[k], jnum), mgp=c(1.8, 0, 0))

      if(k == 1)
        mtext(side=3, at=mean(plim[1:2]), line=0.8, cex=0.9, col=maincolor,
              paste(jnum, "vs.", inum), font=2)

      tp <- paste0(substr(tissues[k], 1, 1), substr(tissue, 1, 1))
      if(!(tp %in% colnames(expr.corr)))
        tp <- paste0(substr(tissue, 1, 1), substr(tissues[k], 1, 1))

      wh <- which(expr.corr[,tp] > 0.75)
      points(get(arr[tissue.index])[wh,i], get(arr[k])[wh,j], bg=sexcolor[2], pch=21)

      rect(plim[1], plim[3], plim[2], plim[4], col=NA)
    }
  }
}
dev.off()
