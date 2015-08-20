# illustration of the swaps of expression arrays

source("colors.R")
attach("../Analysis/R/tissue_text.RData")

arraySwaps <- vector("list", length(tissues))
names(arraySwaps) <- tissues

arraySwaps$adipose <- list(c("Mouse3583", "Mouse3584"), # all swaps
                           c("Mouse3187", "Mouse3188", "Mouse3200"))
                               
arraySwaps$gastroc <- list(c("Mouse3655", "Mouse3659")) # swap

arraySwaps$hypo <- list(c("Mouse3179", "Mouse3188"), # all swaps
                        c("Mouse3208", "Mouse3210"),
                        c("Mouse3347", "Mouse3348"),
                        c("Mouse3367", "Mouse3369"),
                        c("Mouse3381", "Mouse3382"),
                        c("Mouse3449", "Mouse3451"),
                        c("Mouse3452", "Mouse3454"),
                        c("Mouse3589", "Mouse3590"),
                        c("Mouse3592", "Mouse3594"))

arraySwaps$islet <- list(c("Mouse3295", "Mouse3296"), # right is dup
                         c("Mouse3598", "Mouse3599")) # swap
                         

arraySwaps$kidney <- list(c("Mouse3484", "Mouse3503"), # mixture?
                          c("Mouse3510", "Mouse3523")) # swap

arraySwaps$liver <- list(c("Mouse3136", "Mouse3141"), # right is dup
                         c("Mouse3142", "Mouse3143")) # swap

postscript("../Figs/fig4.eps", height=6.5, width=6.5, pointsize=14,
           onefile=FALSE, horizontal=FALSE)
par(mfcol=c(3,2), mar=c(1.1, 1.1, 3.1, 1.1))
for(i in 1:6) {
  plot(0,0,type="n", xlim=c(0, 100), ylim=c(0, 100),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i", bty="n")
  mtext(side=3, tissues[i], col=maincolor, line=0.5, font=2)
  plim <- par("usr")
  rect(plim[1], plim[3], plim[2], plim[4], col=gray)
      
  if(i==1) {
    x <- c(35, 65)
    xd <- 5
    y <- c(80, 80)
    points(x, y, cex=1.7, lwd=2)
    arrows(x[1]+xd, y[1], x[2]-xd, y[1], len=0.07, code=3, lwd=2)

    text(x[1]-xd, y[1], sub("Mouse", "", arraySwaps$adipose[[1]][1]),
         adj=c(1, 0.5))
    text(x[2]+xd, y[2], sub("Mouse", "", arraySwaps$adipose[[1]][2]),
         adj=c(0, 0.5))


    x <- c(50, 65, 35)
    xd <- 5
    y <- c(50, 20, 20)
    yd <- 3
    points(x, y, cex=1.7, lwd=2)
    m <- diff(y[1:2])/diff(x[1:2])
    # following arrows were pointed the wrong way; code=1 fixes it by putting head at the beginning
    arrows(x[1]+xd*0.6, y[1]+m*xd*0.6, x[2]-xd*0.6, y[2]-m*xd*0.6, len=0.07, lwd=2, code=1)
    arrows(x[2]-xd, y[2], x[3]+xd, y[3], len=0.07, lwd=2, code=1)
    arrows(x[3]+xd*0.6, y[3]-m*xd*0.6, x[1]-xd*0.6, y[1]+m*xd*0.6, len=0.07, lwd=2, code=1)

    text(x[1]+xd, y[1], sub("Mouse", "", arraySwaps$adipose[[2]][1]),
         adj=c(0, 0.5))
    text(x[2]+xd, y[2], sub("Mouse", "", arraySwaps$adipose[[2]][2]),
         adj=c(0, 0.5))
    text(x[3]-xd, y[3], sub("Mouse", "", arraySwaps$adipose[[2]][3]),
         adj=c(1, 0.5))
  }

  if(i==2) {
    x <- c(35, 65)
    xd <- 5
    y <- c(50, 50)
    points(x, y, cex=1.7, lwd=2)
    arrows(x[1]+xd, y[1], x[2]-xd, y[1], len=0.07, code=3, lwd=2)

    text(x[1]-xd, y[1], sub("Mouse", "", arraySwaps$gastroc[[1]][1]),
         adj=c(1, 0.5))
    text(x[2]+xd, y[2], sub("Mouse", "", arraySwaps$gastroc[[1]][2]),
         adj=c(0, 0.5))

  }

  if(i==4 || i==6) {
    x <- c(35, 65)
    xd <- 5
    y <- c(65, 65)
    points(x, y, cex=1.7, lwd=2, col=c(distcolor[3], "black"))
    points(x[1], y[1], cex=0.8, pch=16)
    arrows(x[1]+xd, y[1], x[2]-xd, y[1], len=0.07, lwd=2)

    text(x[1]-xd, y[1], sub("Mouse", "", arraySwaps[[i]][[1]][1]),
         adj=c(1, 0.5))
    text(x[2]+xd, y[2], sub("Mouse", "", arraySwaps[[i]][[1]][2]),
         adj=c(0, 0.5))


    y <- c(35, 35)
    points(x, y, cex=1.7, lwd=2)
    arrows(x[1]+xd, y[1], x[2]-xd, y[1], len=0.07, lwd=2, code=3)

    text(x[1]-xd, y[1], sub("Mouse", "", arraySwaps[[i]][[2]][1]),
         adj=c(1, 0.5))
    text(x[2]+xd, y[2], sub("Mouse", "", arraySwaps[[i]][[2]][2]),
         adj=c(0, 0.5))
    
  }
  if(i==5) {
    x <- c(35, 65)
    xd <- 5
    y <- c(65, 65)
    points(x, y, cex=1.7, lwd=2, col=distcolor[3])
    points(x, y, cex=0.8, pch=16)
    arrows(x[1]+xd, y[1], x[2]-xd, y[1], len=0.07, code=3, lwd=2)

    text(x[1]-xd, y[1], sub("Mouse", "", arraySwaps[[i]][[1]][1]),
         adj=c(1, 0.5))
    text(x[2]+xd, y[2], sub("Mouse", "", arraySwaps[[i]][[1]][2]),
         adj=c(0, 0.5))

    text(50, 65+7, "?", col=distcolor[3], font=2)


    y <- c(35, 35)
    points(x, y, cex=1.7, lwd=2)
    arrows(x[1]+xd, y[1], x[2]-xd, y[1], len=0.07, lwd=2, code=3)

    text(x[1]-xd, y[1], sub("Mouse", "", arraySwaps[[i]][[2]][1]),
         adj=c(1, 0.5))
    text(x[2]+xd, y[2], sub("Mouse", "", arraySwaps[[i]][[2]][2]),
         adj=c(0, 0.5))
    
  }

  if(i==3) { # hypo!

    xx <- 16
    xd <- 3

    for(j in 1:9) {
      if(j <=5) {
        y <- rep(100-100/6*j, 2)
        x <- c(xx, 50-xx)
      }
      else {
        y <- rep(100-100/6*(j-5)-100/12, 2)
        x <- c(50+xx, 100-xx)
      }
      points(x, y, cex=1.7, lwd=2)
      arrows(x[1]+xd, y[1], x[2]-xd, y[2], code=3, lwd=2, len=0.07)

      text(x[1]-xd, y[1], sub("Mouse", "", arraySwaps[[i]][[j]][1]),
           adj=c(1, 0.5))
      text(x[2]+xd, y[2], sub("Mouse", "", arraySwaps[[i]][[j]][2]),
           adj=c(0, 0.5))
    

    }

  }


}
dev.off()

detach(2)
