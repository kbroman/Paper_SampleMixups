load("Rcache/f2gc.RData")
load("Rcache/tufted_scan.RData")
#tmp <- lodint(out.tufted.n, chr=17)[,-3]
#tufted.qtl <- interpPositions(tmp, newmap, pmap)
tufted.qtl <- c(26376132, 28359193)
load("../OrigData/annot.final.RData")
tuftqtl.probes <- annot[annot$chromosome=="chr17" & ((annot$end > tufted.qtl[1] & annot$end < tufted.qtl[2]) |
                          (annot$start > tufted.qtl[1] & annot$start < tufted.qtl[2])),"a_gene_id"]

load("Rcache/mlratios_revised.RData")
load("tissue_text.RData")
mlr <- vector("list", length(tissues))
names(mlr) <- tissues
for(i in tissues) {
  id <- findCommonID(f2gc$pheno$MouseNum, rownames(get(arr[i])))
  mlr[[i]] <- list(tuft=tufted.n[id$first],
                   mlratio=get(arr[i])[id$second,tuftqtl.probes])
}

logp <- lapply(mlr, function(a) apply(a[[2]], 2, function(a,b) -log10(t.test(a ~ b)$p.value), a[[1]]))
logp.m <- matrix(unlist(logp), ncol=length(tissues))
dimnames(logp.m) <- list(names(logp[[1]]), names(logp))
logp <- logp.m
rm(logp.m)

plot(mlr[[3]][[2]][,"10002912429"], mlr[[3]][[2]][,"10002916945"], col=c("blue", "red")[mlr[[3]][[1]]+1])
