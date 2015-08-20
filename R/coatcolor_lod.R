# lod curves for coat color traits, before and after corrections
source("colors.R")
attach("../Analysis/R/Rcache/agouti_scan.RData")
attach("../Analysis/R/Rcache/tufted_scan.RData")
source("my_plot_scanone.R")

postscript("../SuppFigs/figS19.eps", height=6.5, width=6.5, onefile=FALSE, horizontal=FALSE, pointsize=12)
par(mar=c(2.8,3.1, 3.1, 0.6), mfrow=c(2,1))
my.plot.scanone(out.agouti.n, out.agouti.o, bgrect=gray, bandcol=darkgray, col=rev(sexcolor),
                yaxt="n", incl.markers=FALSE, ylim=c(0,115), yaxs="i", mgp=c(1.5,0,0),
                hlines=seq(0, 100, by=20))
mtext(side=3, "Agouti coat", font=2, line=0.5, col=maincolor, cex=1.2)
axis(side=2, at=seq(0, 100, by=20), mgp=c(3, 0.35, 0), las=1, tick=FALSE, cex.axis=0.9)
legend("topright", col=rev(sexcolor), lwd=2, c("Corrected", "Original"), bg="white")
plim <- par("usr")
text(plim[1]-diff(plim[1:2])*0.05, plim[4]+diff(plim[3:4])*0.08, "A", font=2, xpd=TRUE, cex=1.3)

my.plot.scanone(out.tufted.n, out.tufted.o, bgrect=gray, bandcol=darkgray, col=rev(sexcolor),
                yaxt="n", incl.markers=FALSE, ylim=c(0,115), yaxs="i", mgp=c(1.5,0,0),
                hlines=seq(0, 100, by=20))
mtext(side=3, "Tufted coat", font=2, line=0.5, col=maincolor, cex=1.2)
axis(side=2, at=seq(0, 100, by=20), mgp=c(3, 0.35, 0), las=1, tick=FALSE, cex.axis=0.9)
plim <- par("usr")
text(plim[1]-diff(plim[1:2])*0.05, plim[4]+diff(plim[3:4])*0.08, "B", font=2, xpd=TRUE, cex=1.3)
dev.off()
