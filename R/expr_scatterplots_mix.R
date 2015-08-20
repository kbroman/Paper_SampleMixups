# scatterplots for apparent mixture in kidney 

source("colors.R")
attach("../Analysis/R/Rcache/mlratios.RData")
attach("../Analysis/R/Rcache/expr_corr.RData")
attach("../Analysis/R/tissue_text.RData")

tissue <- "kidney"
tissue.index <- match(tissue, tissues)
ind <- c("Mouse3484", "Mouse3503")

pairs <- tissuepairs[(tissuepairs[,1]==tissue | tissuepairs[,2]==tissue),"short"]
probes.sametissue <- which(apply(expr.corr[,pairs], 1, median) > 0.75)


postscript("../SuppFigs/figS10.eps", height=8.1, width=6.5, pointsize=12,
           onefile=FALSE, horizontal=FALSE)
par(mfcol=c(6,4), mar=c(3.5,3.6,0.6,1.1), oma=c(0, 0, 1.6, 0))
for(i in ind) {
  for(j in ind) {
    for(k in seq(along=tissues)) {
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
      
      if(k==tissue.index)
        wh <- probes.sametissue
      else
        wh <- which(expr.corr[,tp] > 0.75)

      if(k==tissue.index) {
        if(i==j) color <- darkgray
        else color <- sexcolor[1]
      }
      else color <- sexcolor[2]
      
      points(get(arr[tissue.index])[i,wh], get(arr[k])[j,wh], bg=color, pch=21)

      rect(plim[1], plim[3], plim[2], plim[4], col=NA)
    }
  }
}
dev.off()
