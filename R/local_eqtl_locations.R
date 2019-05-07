# plots of local eQTL locations
library(qtl)
tissues <- c("adipose", "gastroc", "hypo", "islet", "kidney", "liver")
source("colors.R")

attach("../Analysis/FinalData/aligned_geno_with_pmap.RData")

attach("../Analysis/R/Rcache/dgve.RData")

minloc <- sapply(pull.map(f2g), min)

source("my_plot_map.R")

postscript("../SuppFigs/figS14.eps", height=8.4, width=6.5, pointsize=12,
           onefile=FALSE, horizontal=FALSE)
par(mar=c(3.0,3.3,0.6,0.6))
my.plot.map(f2g, chr=1:19, tichwidth=0.1, xlim=c(0.5, 19.9), main="", darkgray=darkgray)

for(i in seq(along=tissues)) {
  eqtl <- matrix(as.numeric(unlist(strsplit(colnames(attr(dgve[[i]], "obsg")), "@"))), ncol=2, byrow=TRUE)
  eqtl <- data.frame(chr=factor(eqtl[,1], 1:19), pos=eqtl[,2])
  rownames(eqtl) <- as.character(1:nrow(eqtl))
  eqtl <- qtl::interpPositions(eqtl, pmap, pull.map(f2g))
  eqtl[,1] <- as.numeric(as.character(eqtl[,1]))
  points(eqtl[,1]+0.12*i, eqtl[,3]-minloc[eqtl[,1]], bg=CCcolor[c(1:3,5,6,8)[i]], pch=21, cex=0.6)
}

legend("bottomright", pch=21, pt.bg=CCcolor[c(1:3,5,6,8)], tissues, bg="white")
dev.off()
