# address question from Reviewer 2 about mismatches in coat color

attach("../Analysis/R/Rcache/agouti_tab.RData")
attach("../Analysis/R/Rcache/tufted_tab.RData")



attach("../Analysis/R/Rcache/dgve_min_and_self.RData")
library(qtlcharts)
mice <- names(mn[[7]])
grp <- rep(1, length(mice))
names(grp) <- mice
grp[agouti.mismatch] <- 2
grp[tufted.mismatch] <- 3

x <- 1-mn[[7]]
y <- 1-self[[7]]
z <- 1-sec[[7]]
o <- order(grp)

iplot(x[o], y[o], indID=mice[o], group=grp[o])

iplot(x[grp==2], y[grp==2], indID=mice[grp==2],
      chartOpts=list(xlim=c(0.5, 1), ylim=c(0.15, 1)))
iplot(x[grp==3], y[grp==3], indID=mice[grp==3],
      chartOpts=list(xlim=c(0.5, 1), ylim=c(0.15, 1)))


iplot(x[o], z[o], indID=mice[o], group=grp[o],
      chartOpts=list(xlim=c(0, 1), ylim=c(0, 1)))

iplot(x[grp==2], y[grp==2], indID=mice[grp==2],
      chartOpts=list(xlim=c(0.5, 1), ylim=c(0.15, 1)))
iplot(x[grp==3], y[grp==3], indID=mice[grp==3],
      chartOpts=list(xlim=c(0.5, 1), ylim=c(0.15, 1)))
