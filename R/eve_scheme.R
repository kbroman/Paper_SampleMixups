######################################################################
# expression vs. expression
######################################################################

source("colors.R")
attach("../Analysis/OrigData/F2.mlratio.kidney.RData")
attach("../Analysis/OrigData/F2.mlratio.liver.RData")
attach("../Analysis/R/Rcache/expr_corr.RData")
attach("../Analysis/R/Rcache/expr_corr_betw_tissues.RData")
dkl <- d$liver$kidney
library(lineup)
id <- findCommonID(colnames(kidney.mlratio), colnames(liver.mlratio))
set.seed(73952596)
transcript <- sample(which(expr.corr[,"kl"] > 0.75), 1)

filetypes <- c("pdf", "eps")
for(ft in filetypes) {
  if(ft=="eps") postscript("../Figs/fig1.eps", width=6.5, height=4.7, pointsize=10, onefile=FALSE, horizontal=FALSE)
  else pdf("../Figs/fig1.pdf", width=6.5, height=4.7, pointsize=10)

  par(mar=rep(0.1, 4), bty="n")

  plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,100), ylim=c(107,0),
       xaxs="i", yaxs="i")

  xw <- 22
  xgap <- c(7, 8)
  x1 <- c(xgap[1], xgap[1]+xw)
  x2 <- c(x1[2]+xgap[2], x1[2]+xgap[2]+xw)
  ygap <- 10
  yh <- c(28, 25)
  y1 <- c(ygap, ygap+yh[1])
  y2 <- c(ygap, ygap+yh[2])

  ytgap <- 3
  xtgap <- 2

  rect(x1[1], y1[1], x1[2], y2[2])
  rect(x1[1], y2[2]+0.5, x1[2], y1[2]+0.5)
  text(mean(x1), y1[1]-ytgap, "expression in kidney", cex=1.2, col=maincolor, font=2)
  text(x1[1]-xtgap, mean(y1), "kidney samples", srt=90)
  text(mean(x1), y1[2]+0.5+ytgap, "array probes")

  rect(x2[1], y1[1], x2[2], y2[2])
  rect(x2[1], y2[2]+0.5, x2[2], y1[2]+0.5)
  text(mean(x2), y1[1]-ytgap, "expression in liver", cex=1.2, col=maincolor, font=2)
  text(x2[1]-xtgap, mean(y1), "liver samples", srt=90)
  text(mean(x2), y1[2]+ytgap, "array probes")

  recw <- 1.2
  rect(x1[1]+diff(x1)*0.2, y1[1], x1[1]+diff(x1)*0.2+recw, y2[2],
       col=CCcolor[1])
  rect(x2[1]+diff(x2)*0.2, y1[1], x2[1]+diff(x2)*0.2+recw, y2[2],
       col=CCcolor[1])

  arrowgap <- 4.5
  arroww <- 5
  arrows(x2[2]+arrowgap, mean(y1), x2[2]+arrowgap+arroww, mean(y1), len=0.1, angle=15, col=maincolor, lwd=2)

  text(x1[1]-5, y1[1]-5, "A", font=2, cex=1.3)


  x3 <- c(99-xw, 99)
  rect(x3[1], y1[1], x3[2], y1[2], col=gray)
  text(x3[1]-5, y1[1]-5, "B", font=2, cex=1.3)
  text(x3[1]-xtgap, mean(y1), "kidney expression", srt=90)
  text(mean(x3), y1[2]+ytgap, "liver expression")

  x <- liver.mlratio[transcript,id$second]
  y <- kidney.mlratio[transcript,id$first]
  x <- x/diff(range(x))/1.3*diff(x3)+mean(x3)
  y <- mean(y1)-(y/diff(range(y))/1.3*diff(y1))
  points(x, y, bg=sexcolor[1], pch=21, cex=0.8)
  text(mean(x3), y1[1]-ytgap, paste("corr =", sprintf("%.2f", cor(x, -y, use="complete"))),
       cex=1.2, col=maincolor, font=2)


  y3 <- y1+y1[2]+ygap+10
  y4 <- y3
  rect(x1[1], y3[1], x1[2], y3[2])
  text(mean(x1), y3[1]-ytgap, "expression in kidney", cex=1.2, col=maincolor, font=2)
  text(x1[1]-xtgap, mean(y3), "kidney samples", srt=90)
  text(mean(x1), y3[2]+0.5+ytgap, "array probes")


  rect(x2[1], y4[1], x2[2], y4[2])
  rect(x2[1], y4[1], x2[2], y4[2])
  text(mean(x2), y4[1]-ytgap, "expression in liver", cex=1.2, col=maincolor, font=2)
  text(x2[1]-xtgap, mean(y4), "liver samples", srt=90)
  text(mean(x2), y4[2]+ytgap, "array probes")

  text(x1[1]-5, y3[1]-5, "C", font=2, cex=1.3)

  arrows(x2[2]+arrowgap, mean(y3), x2[2]+arrowgap+arroww, mean(y3), len=0.1, angle=15, col=maincolor, lwd=2)

  grayw <- diff(x1)*0.25

  rect(x1[1], y3[1], x1[1]+grayw, y3[2], col=gray)
  rect(x2[1], y4[1], x2[1]+grayw, y4[2], col=gray)


  rech <- 1.3
  rect(x1[1], y3[1]+diff(y3)*0.3, x1[1]+grayw, y3[1]+diff(y3)*0.3+rech, col=CCcolor[1])
  rect(x2[1], y4[1]+diff(y4)*0.6, x2[1]+grayw, y4[1]+diff(y4)*0.6+rech, col=sexcolor[2])


  x3 <- c(99-xw, 99)
  y5 <- c(y3[1], y3[1]+diff(y3)*1.2)
  rect(x3[1], y5[1], x3[2], y5[2])
  text(x3[1]-5, y5[1]-5, "D", font=2, cex=1.3)
  text(x3[1]-xtgap, mean(y5), "kidney samples", srt=90)
  text(mean(x3), y5[2]+ytgap, "liver samples")
  text(mean(x3), y5[1]-ytgap, "similarity matrix", cex=1.2, col=maincolor, font=2)

  lnames <- rownames(dkl)
  knames <- colnames(dkl)
  mice <- c(paste0("Mouse", c("3179", "3208", "3210", "3347", "3510", "3523", "3136", "3141")),
            sample(knames[is.na(match(knames, lnames))], 1),
            sample(lnames[is.na(match(lnames, knames))], 1))
  z <- matrix(ncol=length(mice), nrow=length(mice))
  dimnames(z) <- list(mice, mice)
  m <- match(mice, lnames)
  for(i in mice) {
    if(i %in% knames)
      z[i,!is.na(m)] <- dkl[m[!is.na(m)],i]
  }
                                        # modify similarity matrix for illustration purposes
  z[!is.na(z)] <- (z[!is.na(z)]+1)/2
  z[!is.na(z) & z < 0.82] <- z[!is.na(z) & z < 0.82]*0.4

                                        # colors
  blues<-colorRampPalette(c("white","blue"))(256)
  zcol <- matrix("", ncol=ncol(z), nrow=nrow(z))
  zcol[is.na(z)] <- "orange"
  zcol[!is.na(z)] <- blues[z[!is.na(z)]*255+1]

  xgap <- 0.3
  x <- x3[1]+xgap + (0:9)*(diff(x3)-xgap*2)/10
  xn <- x + diff(x[1:2])
  ygap <- 0.5
  y <- y5[1]+ygap + (0:9)*(diff(y5)-2*ygap)/10
  yn <- y + diff(y[1:2])
  for(i in seq(along=y))
    for(j in seq(along=x))
      rect(x[j], y[i], xn[j], yn[i], col=zcol[11-i,j], border="white")

  rect(x3[1], y5[1], x3[2], y5[2])


  dev.off()

}
