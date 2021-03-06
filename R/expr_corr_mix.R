# plots of inter-tissue correlations for apparent sample mixture in kidney

source("colors.R")
attach("../Analysis/R/tissue_text.RData")
tissue <- "kidney"
tissue.index <- match(tissue, tissues)
hilitpairs <- (tissuepairs[,1]==tissue | tissuepairs[,2]==tissue)

file <- "Rcache/expr_corr_mix.RData"
if(file.exists(file)) {
  load(file)
} else {
  attach("../Analysis/R/Rcache/mlratios.RData")
  attach("../Analysis/R/Rcache/expr_corr.RData")

  ind <- c("Mouse3484", "Mouse3503")

  # calculate correlations
  corr.mix <- matrix(nrow=4, ncol=nrow(tissuepairs))
  dimnames(corr.mix) <- list(as.character(1:4), tissuepairs[,"short"])
  kk <- 1
  for(i in ind) {
    for(j in ind) {
      rownames(corr.mix)[kk] <- paste(i, "vs.", j)
      for(k in 1:nrow(tissuepairs)) {
        wh <- which(expr.corr[,k] > 0.75)
        corr.mix[kk,k] <- cor(get(arr[tissuepairs[k,1]])[wh,i],
                          get(arr[tissuepairs[k,2]])[wh,j],
                          use="complete")
      }
      kk <- kk + 1
    }
  }

  corr.mix.onetissue <- rep(0, 6)
  for(k in 1:6) {
    wh <- which(tissuepairs[,1]==tissues[k] | tissuepairs[,2]==tissues[k])
    wh <- which(apply(expr.corr[,wh], 1, max) > 0.75)
    corr.mix.onetissue[k] <- cor(get(arr[k])[wh,ind[1]],
                                 get(arr[k])[wh,ind[2]],
                                 use="complete")
  }

  # this vs that as one vector
  corr.mix.ivj <- c(corr.mix[2,], corr.mix[3,], corr.mix.onetissue)

  # assign it names
  cn <- colnames(corr.mix)
  rcn <- sapply(strsplit(cn, ""), function(a) paste(rev(a), collapse=""))
  names(corr.mix.ivj) <- c(cn, rcn, paste0(tissues.short, tissues.short))
  corr.mix.ivj <- corr.mix.ivj[order(names(corr.mix.ivj))]

  save(corr.mix, corr.mix.ivj, file=file)

  for(i in 1:2) detach(2)
}

colors <- c(rev(sexcolor), distcolor[2], CCcolor[8])

postscript("../SuppFigs/figS11.eps", height=6.5, width=6.5, pointsize=12,
           onefile=FALSE, horizontal=FALSE)
layout(cbind(c(1,3),c(2,3)))
par(mar=c(4.0,3.6,2.1,2.1))
for(i in c(1,4)) {
  plot(0,0,
       xlab="", ylab="", xlim=c(0.25,15.75), ylim=c(0.31,1.01),
       xaxt="n", yaxt="n", xaxs="i", yaxs="i",
       type="n")

  mtext(side=3, rownames(corr.mix)[i], col=maincolor, line=0.5, font=2)

  plim <- par("usr")
  rect(plim[1], plim[3], plim[2], plim[4], col=gray)

  thesecolors <- colors[c(1, ifelse(i==1, 2, 3))]

  thesecolors.text <- thesecolors
  for(j in 1:2) {
    tmp <- col2rgb(thesecolors[j])
    thesecolors.text[j] <- rgb(tmp[1]*0.8, tmp[2]*0.8, tmp[3]*0.8, maxColorValue=255)
  }

  axis(side=1, at=(1:15)[!hilitpairs], tissuepairs[!hilitpairs,"short"], tick=FALSE,
       mgp=c(3, 0.35, 0), las=2, col.axis=thesecolors.text[1])
  axis(side=1, at=(1:15)[hilitpairs], tissuepairs[hilitpairs,"short"], tick=FALSE,
       mgp=c(3, 0.35, 0), las=2, col.axis=thesecolors.text[2])
  abline(v=1:15, col=darkgray)
  ypos <- seq(0.4, 1, by=0.1)
  axis(side=2, at=ypos, mgp=c(3, 0.35, 0), tick=FALSE, las=1)
  abline(h=ypos, col="white")
  title(xlab="Tissue pair", mgp=c(2.1, 0, 0))
  title(ylab="Correlation", mgp=c(2.1, 0, 0))

  points((1:15)[!hilitpairs], corr.mix[i,!hilitpairs], pch=21, bg=thesecolors[1])
  points((1:15)[hilitpairs], corr.mix[i,hilitpairs], pch=21, bg=thesecolors[2])

  rect(plim[1], plim[3], plim[2], plim[4], col=NA)
}

par(mar=c(3.5,3.6,2.6,2.1))

nam <- names(corr.mix.ivj)
grp <- rep(1, length(corr.mix.ivj))
first <- sapply(strsplit(nam, ""), "[", 1)
second <- sapply(strsplit(nam, ""), "[", 2)
tissue.short <- tissues.short[tissue.index]
grp[first==tissue.short] <- 2
grp[second==tissue.short] <- 3
grp[first==tissue.short & second==tissue.short] <- 4

plot(0,0,
     xlab="", ylab="", xlim=c(0.25,36.75), ylim=c(0.31,1.01),
     xaxt="n", yaxt="n", xaxs="i", yaxs="i",
     type="n")

mtext(side=3, rownames(corr.mix)[2], col=maincolor, line=0.5, font=2)

plim <- par("usr")
rect(plim[1], plim[3], plim[2], plim[4], col=gray)

thesecolors <- colors
thesecolors.text <- thesecolors
for(j in 1:4) {
  tmp <- col2rgb(thesecolors[j])
  thesecolors.text[j] <- rgb(tmp[1]*0.8, tmp[2]*0.8, tmp[3]*0.8, maxColorValue=255)
}

for(i in 1:4)
  axis(side=1, at=(1:36)[grp==i], nam[grp==i], tick=FALSE,
       mgp=c(3, 0.35, 0), las=2, col.axis=thesecolors.text[i])

abline(v=1:36, col=darkgray)
ypos <- seq(0.4, 1, by=0.1)
axis(side=2, at=ypos, mgp=c(3, 0.35, 0), tick=FALSE, las=1)
abline(h=ypos, col="white")
title(xlab="Tissue pair", mgp=c(2.1, 0, 0))
title(ylab="Correlation", mgp=c(2.1, 0, 0))

for(i in 1:4)
  points((1:36)[grp==i], corr.mix.ivj[grp==i], pch=21, bg=thesecolors[i])

rect(plim[1], plim[3], plim[2], plim[4], col=NA)


dev.off()
