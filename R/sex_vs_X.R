# figure with X chromosome genotypes inconsistent with sex

rm(list=ls())
source("colors.R")

attach("../Analysis/OrigData/rawg.RData") # contains notfemale and notmale
attach("../Analysis/OrigData/clean_cross.RData") # contains notfemale and notmale
necropsy <- read.csv("../Analysis/FinalData/necropsy.csv", as.is=TRUE)
sex <- sub("F", "Female", sub("M", "Male",necropsy$Sex))
names(sex) <- necropsy$MouseNum

rawg <- rawg["X", rawg$pheno$Strain == "F2"]
rawg$pheno$Sex <- sex[as.character(rawg$pheno$MouseNum)] # fix sex

mnr <- markernames(rawg)
mnc <- markernames(f2g)
rawg <- pull.markers(rawg, mnr[!is.na(match(mnr, mnc))])

female <- (rawg$pheno$Sex=="Female")
female.prob <- apply(rawg[,female]$geno$X$data, 1, function(a) sum(a==1, na.rm=TRUE))
names(female.prob) <- as.character(rawg$pheno$MouseNum)[female]
female.prob <- sort(female.prob[female.prob > 0])

male <- (rawg$pheno$Sex=="Male")
male.prob <- apply(rawg[,male]$geno$X$data, 1, function(a) sum(a==2, na.rm=TRUE))
names(male.prob) <- as.character(rawg$pheno$MouseNum)[male]
male.prob <- sort(male.prob[male.prob > 0])

nf <- length(female.prob)
nm <- length(male.prob)
map <- rawg$geno$X$map
map <- map - min(map)
jit <- c(0, -0.15, 0.15, 0, -0.4, -0.1, +0.2, 0.35, 0, -0.25, -0.05, +0.3, 0.15, 0, 0, 0, -0.35, 0.05, 0.15, 0)*2
omap <- map
map <- map + jit


postscript("../SuppFigs/figS2.eps", height=7, width=8.5, horizontal=TRUE, onefile=FALSE, pointsize=12)

mappad <- 1
mappad2 <- 7
par(mar=c(3.1,2.6,3.1,5.1))
plot(0, 0, type="n", xlab="", xaxt="n", ylab="", yaxt="n", bty="n",
     xlim=c(-mappad, max(map)+mappad), ylim=c(nf+nm+1.75, 0.25), xaxs="i", yaxs="i")
plim <- par("usr")
rect(plim[1], plim[3], plim[2], nf+1.25, col=gray)
rect(plim[1], plim[4], plim[2], nf+0.75, col=gray)
axis(side=1, at=seq(0, 70, by=10), mgp=c(3, 0.2, 0), tick=FALSE, las=1)
title(xlab="Location (cM)", mgp=c(1.8, 0, 0))
abline(v=seq(0, 70, by=10), col="white")
yf <- 1:nf
segments(rep(0,nf), yf, rep(max(map), nf), yf)

ym <- (1:nm)+nf+1
segments(rep(0,nf), ym, rep(max(map), nf), ym)
abline(h=nf+1.25)
fx <- rawg$geno$X$data[match(names(female.prob), rawg$pheno$MouseNum),]
for(i in 1:nrow(fx)) {
  notna <- !is.na(fx[i,])
  points(map[notna], rep(yf[i], length(map))[notna], pch=21, bg=f2color[fx[i,notna]])
  text(max(map)+mappad2, yf[i], female.prob[i], xpd=TRUE, adj=c(1, 0.5))
}

mx <- rawg$geno$X$data[match(names(male.prob), rawg$pheno$MouseNum),]
for(i in 1:nrow(mx)) {
  notna <- !is.na(mx[i,])
  points(map[notna], rep(ym[i], length(map))[notna], pch=21, bg=f2color[mx[i,notna]])
  text(max(map)+mappad2, ym[i], male.prob[i], xpd=TRUE, adj=c(1, 0.5))
}

mtext(side=2, at=mean(yf), line=1, "Nominally female", font=2)
mtext(side=2, at=mean(ym), line=1, "Nominally male", font=2)
for(i in 1:2) detach(2)
mtext(side=3, at=max(map)+mappad2*0.95, "number\nincompatible", line=0.5, font=2)

xpad <- 1
genolabels <- c("BB or BY", "BR", "RR or RY")
x <- c(26.2, 37.7, 43.8)
y <- rep(plim[4]+diff(plim[3:4])*0.05, 3)
points(x, y, bg=f2color, pch=21, xpd=TRUE, cex=1.3)
text(x+xpad, y, genolabels, adj=c(0, 0.5), xpd=TRUE)

rect(plim[1], plim[3], plim[2], nf+1.25, col=NA, lend=1, ljoin=1)
rect(plim[1], plim[4], plim[2], nf+0.75, col=NA, lend=1, ljoin=1)

dev.off()
