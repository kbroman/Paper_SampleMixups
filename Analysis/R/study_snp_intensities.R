x <- read.csv("snp_intensities.csv")
avesigsum <- tapply(log2(x$SignalSum), x$MouseNum, function(a) mean(a[is.finite(a)]))
avesigsum <- avesigsum[rownames(calls)]
hilit <- rownames(calls)[grep("omit", calls[,"DNA"])]

hist(avesigsum, breaks=100)
rug(avesigsum)
rug(avesigsum[hilit], col="green", lwd=2)

