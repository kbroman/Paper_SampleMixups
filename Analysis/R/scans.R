load("../FinalData/aligned_geno_with_pmap.RData")
load("../FinalData/lipomics_final_rev2.RData")

lipomics <- lipomics[match(f2g$pheno$MouseNum, lipomics$MouseNum),]
tuft <- match(lipomics$TUFT, c("N", "Y"))-1
agouti <- match(lipomics$AGOUTI, c("T","B"))-1

f2g <- calc.genoprob(f2g, step=0.5, err=0.002, map.function="c-f", stepwidth="max")
out.agouti <- scanone(f2g, phe=agouti, model="binary")
out.tuft <- scanone(f2g, phe=tuft, model="binary")

p.agouti <- f2g$geno[[2]]$prob[,"loc75",]
p.tuft <- f2g$geno[[17]]$prob[,"rs3700924",]


tab.agouti <- table(apply(p.agouti, 1, function(a) if(max(a)> 0.99) return(which(a==max(a))) else return(NA) ), agouti) # 7/480 mismatches
tab.tuft <- table(apply(p.tuft, 1, function(a) if(max(a)> 0.99) return(which(a==max(a))) else return(NA) ), tuft)       # 4/519 mismatches

######################################################################

out <- scanone(f2g, pheno.col=log(lipomics$"insulin 10 wk"), intcovar=as.numeric(f2g$pheno$Sex))
out2 <- scanone(f2g, pheno.col=log(lipomics$"INSULIN (ng/ml) 10 wk"), intcovar=as.numeric(f2g$pheno$Sex))
out3 <- scanone(f2g, pheno.col=log(lipomics$"insulin 10 wk") + log(lipomics$"INSULIN (ng/ml) 10 wk"), intcovar=as.numeric(f2g$pheno$Sex))

######################################################################

