load("Rcache/mlratios_revised.RData")
load("tissue_text.RData")

library(lineup)

load("Rcache/expr_corr_rev.RData") # pairwise expr corr with realigned data

id <- findCommonID(rownames(adipose.mlratio), rownames(islet.mlratio))

pval.ai <- apply((adipose.mlratio[id$first,] - islet.mlratio[id$second,]), 2, function(a) log10(t.test(a)$p.value))

o <- order(pval.ai)[1:150]
x <- rbind(adipose.mlratio[,o], islet.mlratio[,o])

x[x == -2 | x == 2] <- NA
corr <- cor(t(x), use="pairwise")
plot(colMeans(corr[-(1:493), (1:493)]))
plot(colMeans(corr[(1:493), -(1:493)]))
plot(colMeans(corr[(1:493), (1:493)]))
cm <- colMeans(corr[(1:493), (1:493)])

######################################################################


load("Rcache/expr_corr_rev.RData") # pairwise expr corr with realigned data
expr.corr.med <- apply(expr.corr, 1, median)
expr.corr.max <- apply(expr.corr, 1, max)

for(i in tissues) {
  load(paste0("../OrigData/F2.int2.", i, ".RData"))
  x <- paste0(i, ".int2")
  assign(x, log2(get(x)))
}

pair.corr <- vector("list", nrow(tissuepairs))
names(pair.corr) <- tissuepairs$short
for(i in 1:nrow(tissuepairs)) {
  cat(i, "\n")
  x <- get(paste0(tissuepairs[i,1], ".int2"))[expr.corr.max > 0.5,]
  y <- get(paste0(tissuepairs[i,2], ".int2"))[expr.corr.max > 0.5,]

  colnames(x) <- sub("^Mouse", tissuepairs[i,1], colnames(x))
  colnames(y) <- sub("^Mouse", tissuepairs[i,2], colnames(y))

  z <- cbind(x,y)
  pair.corr[[i]] <- cor(z, use="pair")
}


par(mfrow=c(3,5), las=1)
for(i in 1:nrow(tissuepairs)) {
  tis <- substr(colnames(pair.corr[[i]]), 1, 1)
  utis <- unique(tis)
  x <- rowMeans(pair.corr[[i]][,tis==utis[1]])
  y <- rowMeans(pair.corr[[i]][,tis==utis[2]])
  plot(x, y,
       xlab=paste("corr vs", utis[1]),
       ylab=paste("corr vs", utis[2]), xlim=c(0.3, 1), ylim=c(0.3, 1),
       col=c("blue","red")[(tis==utis[1])+1], lwd=2, cex=0.8)
}
