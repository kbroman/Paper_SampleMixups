# plots of Xist and Y chr expression, colored by sex, before and after correction of sample mixups
tissues <- c("adipose", "gastroc", "hypo", "islet", "kidney", "liver")
source("colors.R")
library(lineup)
library(broman)


file <- "Rcache/xist_and_y.RData"
if(file.exists(file)) {
  load(file)
} else {
  attach("../Analysis/OrigData/annot.amit_rev.RData")
  xist.probe <- annot$a_gene_id[!is.na(annot$officialgenesymbol) & annot$officialgenesymbol=="Xist"]
  ychr.probes <- annot$a_gene_id[annot$chr=="Y"]
  goodyprobes <- c("10002897002", "512831", "10002915709")
  necropsy <- read.csv("../Analysis/FinalData/necropsy.csv", as.is=TRUE)
  necropsy$Sex <- sub("F", "Female", sub("M", "Male",necropsy$Sex))

  sex <- necropsy$Sex
  names(sex) <- necropsy$MouseNum

  mlratio.old <- mlratio.new <- vector("list", length(tissues))
  names(mlratio.old) <- names(mlratio.new) <- tissues

  cat("old\n")
  attach("../Analysis/R/Rcache/mlratios.RData")
  for(i in tissues) {
    cat("\t", i, "\n")

    if(i=="hypo") # save indicators of bad arrays
      hypo.bad <- colnames(hypo.mlratio)[apply(hypo.mlratio, 2, function(a) median(a[!is.na(a) & a > -2 & a < 2])) > 0.016]

    mlratio.old[[i]] <- get(paste.(i, "mlratio"))[,c(xist.probe, goodyprobes)]
  }
  detach(2)

  cat("new\n")
  attach("../Analysis/R/Rcache/mlratios_revised.RData")
  for(i in tissues) {
    cat("\t", i, "\n")

    mlratio.new[[i]] <- get(paste.(i, "mlratio"))[,c(xist.probe, goodyprobes)]
  }
  detach(2)

  detach(2)
  save(sex, mlratio.old, mlratio.new, hypo.bad, file=file)
}

sexnum <- match(sex, c("Female", "Male"))

postscript("../SuppFigs/figS12.eps", height=8.4, width=6.5, pointsize=12,
           onefile=FALSE, horizontal=FALSE)
par(mfrow=c(length(tissues),2), oma=c(0,2.1,1.6,0), mar=c(3.1, 3.6, 0.6, 0.6))
for(i in tissues) {
  id.old <- findCommonID(rownames(mlratio.old[[i]]), names(sex))
  x <- mlratio.old[[i]][id.old$first,1]
  y <- rowMeans(mlratio.old[[i]][,-1], na.rm=TRUE)
  plot2nd <- (y > 0.5*x & sex[id.old$second]=="Female") | (y < 0.5*x & sex[id.old$second]=="Male")

  grayplot(x[!plot2nd], y[!plot2nd], xlab="", ylab="", bgcolor=gray, 
           vlines=seq(-2, 1, by=0.5), hlines=seq(-2, 1, by=0.5), ylim=c(-2.1, 1.1), yaxs="i",
           xlim=c(-2.05, 1.05), xaxs="i", mgp.y=c(3, 0.4, 0), mgp.x=c(3, 0.2, 0),
           bg=sexcolor[sexnum[id.old$second]][!plot2nd], pch=21)
  title(xlab=expression(paste(italic("Xist "), plain("expression"))), mgp=c(1.7, 0, 0))
  title(ylab="Y chr expression", mgp=c(2.5, 0, 0))
  plim <- par("usr")
  mtext(side=2, i, at=mean(plim[3:4]), cex=1.1, line=4.1, col="darkslateblue", font=2)
  if(i=="adipose") mtext(side=3, at=mean(plim[1:2]), "Original", line=0.8, cex=1.1, col=maincolor, font=2)
  if(sum(plot2nd) > 0)
    points(x[plot2nd], y[plot2nd], bg=sexcolor[sexnum[id.old$second]][plot2nd], pch=21)
  rect(plim[1], plim[3], plim[2], plim[4], col=NA, border="black", lend=1, ljoin=1)

  id.new <- findCommonID(rownames(mlratio.new[[i]]), names(sex))
  x <- mlratio.new[[i]][id.new$first,1]
  y <- rowMeans(mlratio.new[[i]][,-1], na.rm=TRUE)
  plot2nd <- (y > 0.5*x & sex[id.new$second]=="Female") | (y < 0.5*x & sex[id.new$second]=="Male")

  grayplot(x[!plot2nd], y[!plot2nd], xlab="", ylab="", bgcolor=gray, 
           vlines=seq(-2, 1, by=0.5), hlines=seq(-2, 1, by=0.5), ylim=c(-2.1, 1.1), yaxs="i",
           xlim=c(-2.05, 1.05), xaxs="i", mgp.y=c(3, 0.4, 0), mgp.x=c(3, 0.2, 0),
           bg=sexcolor[sexnum[id.new$second]][!plot2nd], pch=21)
  title(xlab=expression(paste(italic("Xist "), plain("expression"))), mgp=c(1.7, 0, 0))
  title(ylab="Y chr expression", mgp=c(2.5, 0, 0))
  if(i=="adipose") mtext(side=3, at=mean(plim[1:2]), "Corrected", line=0.8, cex=1.1, col=maincolor, font=2)
  if(sum(plot2nd) > 0)
    points(x[plot2nd], y[plot2nd], bg=sexcolor[sexnum[id.new$second]][plot2nd], pch=21)
  rect(plim[1], plim[3], plim[2], plim[4], col=NA, border="black", lend=1, ljoin=1)

}
dev.off()
