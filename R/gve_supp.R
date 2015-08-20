######################################################################
# eQTL genotype vs expression: supplemental figures
######################################################################

source("colors.R")
library(lineup)


attach("../Analysis/R/Rcache/dgve.RData")
y <- attr(dgve[["islet"]], "y")
obsg <- attr(dgve[["islet"]], "obsg")
infg <- attr(dgve[["islet"]], "infg")

attach("../Analysis/R/Rcache/pmark.RData")

gnames <- attr(dgve[["islet"]], "genonames")

postscript("../SuppFigs/figS13.eps", height=6.5, width=6.5, onefile=FALSE, horizontal=FALSE, pointsize=12)
par(mfrow=c(2,2), mar=c(4.1,4.1,2.1,1.1))
selected <- c(62, 3, 11, 22)
for(i in selected) {
  cat("****", i, "****\n")


  e <- y[[i]]
  pmar <- pmark[colnames(e)[1],2]

  og <- obsg[,i]
  ig <- infg[,i]

  id <- findCommonID(rownames(e), names(og))
  e <- e[id$first,]
  ig <- ig[id$first]
  og <- og[id$second]
  plotsecond <- !is.na(ig) & !is.na(og) & og != ig

  plot(e[,1], e[,2], type="n", xlab="", ylab="", xaxt="n", yaxt="n")
  title(ylab=paste("expression of probe", colnames(e)[2]), mgp=c(2.2, 0.5, 0))
  title(xlab=paste("expression of probe", colnames(e)[1]), mgp=c(1.5, 0.5, 0))

  plim <- par("usr")
  rect(plim[1], plim[3], plim[2], plim[4], col=gray)
  xl <- pretty(e[,1])
  if(i==3) { # trying to force it!
    xl <- c(-0.3, -0.1, 0.1, 0.3)
    abline(v=xl, col="white")
    axis(side=1, at=xl, mgp=c(3, 0.1, 0), tick=FALSE)
    xl <- c(-0.2, 0, 0.2)
  }
  abline(v=xl, col="white")
  axis(side=1, at=xl, mgp=c(3, 0.1, 0), tick=FALSE)

  yl <- pretty(e[,2])
  abline(h=yl, col="white")
  axis(side=2, at=yl, mgp=c(3, 0.3, 0), las=1, tick=FALSE)

  text(plim[1]-diff(plim[1:2])*0.15, plim[4]+diff(plim[3:4])*0.08, LETTERS[match(i,selected)], font=2, xpd=TRUE, cex=1.3)


  points(e[!plotsecond,1], e[!plotsecond,2], bg=f2color[og[!plotsecond]], pch=21, cex=0.75)
  points(e[plotsecond,1], e[plotsecond,2], bg=f2color[og[plotsecond]], pch=21, cex=0.75)

  location <- "bottomright"
  if(i==62) location <- "topleft"
  if(i==22) location <- "bottomleft"

  if(i==62) {
    legend(location, pch=21, pt.bg=f2color, title=pmar,
           legend=gnames, bg="white")
  } else {
    pmar.w <- strwidth(pmar)*1.1
    pmar.h <- strheight(pmar)*2
    if(location=="bottomright") {
      rect(plim[2]-pmar.w, plim[3], plim[2], plim[3]+pmar.h, col="white")
      text(plim[2]-pmar.w/2, plim[3]+pmar.h/2, pmar)
    } else {
      rect(plim[1]+pmar.w, plim[3], plim[1], plim[3]+pmar.h, col="white")
      text(plim[1]+pmar.w/2, plim[3]+pmar.h/2, pmar)
    }
  }

  rect(plim[1], plim[3], plim[2], plim[4], col=NA, border="black", lend=1, ljoin=1)
}
for(i in 1:2) detach(2)
dev.off()
