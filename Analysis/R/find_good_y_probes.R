load("../OrigData/annot.amit_rev.RData")
xist.probe <- annot$a_gene_id[!is.na(annot$officialgenesymbol) & annot$officialgenesymbol=="Xist"]
ychr.probes <- annot$a_gene_id[annot$chr=="Y"]

load("../FinalData/aligned_geno_with_pmap.RData")

tissues <- c("adipose", "gastroc", "hypo", "islet", "kidney", "liver")
for(i in tissues)
  load(paste0("../FinalData/", i, "_mlratio_final.RData"))

par(mfrow=c(6,1))
for(i in tissues) {
  tmp <- paste.(i, "mlratio")
  id <- findCommonID(f2g$pheno$MouseNum, rownames(get(tmp)))
  sex <- f2g$pheno$Sex[id$first]
  tmp <- get(tmp)[id$second,ychr.probes]
  bs <- beeswarm(tmp ~ col(tmp), method="center", pch=1, spacing=0.2, cex=0.5, do.plot=FALSE)

  plot(0, 0, type="n", xlim=c(0.5, ncol(tmp)+0.5), ylim=c(-2, 2),
       main=i, xlab="Probe", ylab="ML Ratio")
  for(k in 1:ncol(tmp))
    points(k + runif(nrow(tmp), -0.3, 0.3), tmp[,k], cex=0.5,
           col=c("blue", "red")[match(sex, c("Male", "Female"))])
}  
