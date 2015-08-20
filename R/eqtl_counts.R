# figure with counts of eQTL, before and after correction

library(broman) # for grayplot
source("colors.R")
attach("../Analysis/R/tissue_text.RData")
attach("../Analysis/R/Rcache/neqtl_new.RData")
attach("../Analysis/R/Rcache/neqtl_old.RData")
neqtl.old <- t(sapply(neqtl.old, colSums))
neqtl.new <- t(sapply(neqtl.new, colSums))


postscript("../Figs/fig9.eps", height=4.5, width=6.5,
           pointsize=12, onefile=FALSE, horizontal=FALSE)
par(mfrow=c(1,2), las=1, mar=c(5.1, 4.1, 4.1, 2.1))
x <- c(1:6, 1:6)
y <- c(neqtl.old[,1], neqtl.new[,1])
grayplot(x, y, xlab="", ylab="No. local-eQTL", type="n",
         ylim=c(0, 3600), yat = seq(0, 4000, by=1000), hlines=seq(0, 4000, by=1000),
         vlines=1:6, vlines.col=darkgray, vlines.lwd=2,
         xat=NA, yaxs="i", xlim=c(0.7, 6.3), 
         mgp.x=c(3, 0.3, 0), mgp.y=c(3, 0.2, 0))
segments(x-0.2, y, x+0.2, y, col=rep(sexcolor, rep(6,2)), lwd=2)
plim <- par("usr")
mtext(side=3, "local-eQTL", col=maincolor, font=2, line=0.5, cex=1.2)
axis(side=1, 1:6, tissues, mgp=c(3, 0.35, 0), tick=FALSE, las=2)
mtext(side=1, line=3.6, "Tissue")
#text(plim[1]-diff(plim[1:2])*0.18, plim[4]+diff(plim[3:4])*0.08, LETTERS[1], font=2, xpd=TRUE, cex=1.2)


y <- c(neqtl.old[,2], neqtl.new[,2])
grayplot(x, y, xlab="", type="n",
         ylab=expression(paste(plain("No. "), italic("trans"), plain("-eQTL"))), 
         ylim=c(0, 27500), hlines=seq(0, 25000, by=5000),
         vlines=1:6, xat=NA, yaxs="i", vlines.col=darkgray, vlines.lwd=2,
         xlim=c(0.7, 6.3), mgp=c(3, 0.3, 0))
segments(x-0.2, y, x+0.2, y, col=rep(sexcolor, rep(6,2)), lwd=2)
plim <- par("usr")
mtext(side=3, expression(paste(bolditalic("trans"), bold("-eQTL"))),
      col=maincolor, font=2, line=0.5, cex=1.2)
axis(side=1, 1:6, tissues, mgp=c(3, 0.35, 0), tick=FALSE, las=2)
mtext(side=1, line=3.6, "Tissue")

legend("topleft", lwd=2, col=rev(sexcolor), c("Corrected",
       "Original"), bg="white", cex=0.7, seg.len=0.9)

#text(plim[1]-diff(plim[1:2])*0.18, plim[4]+diff(plim[3:4])*0.08, LETTERS[2], font=2, xpd=TRUE, cex=1.2)

dev.off()
