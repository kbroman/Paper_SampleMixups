# lod curves for insulin, before and after corrections
source("colors.R")
attach("../Analysis/R/Rcache/insulin_scan.RData")
source("my_plot_scanone.R")

postscript("../Figs/fig8.eps", height=4.5, width=6.5, onefile=FALSE, horizontal=FALSE, pointsize=12)
par(mar=c(2.8,3.1, 0.6, 0.6))
my.plot.scanone(out.ins.n, out.ins.o, bgrect=gray, bandcol=darkgray, col=rev(sexcolor),
                yaxt="n", incl.markers=FALSE, ylim=c(0,9.5), yaxs="i", mgp=c(1.5,0,0))
axis(side=2, at=seq(0, 8, by=2), mgp=c(3, 0.35, 0), las=1, tick=FALSE, cex.axis=0.9)
legend("topright", col=rev(sexcolor), lwd=2, c("Corrected", "Original"), bg="white")
dev.off()

detach(2)
