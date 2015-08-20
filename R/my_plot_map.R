my.plot.map <- 
function(x, map2, chr, horizontal=FALSE, shift=TRUE,
         show.marker.names=FALSE, alternate.chrid=FALSE, 
         tickwidth=0.15, ...) 
{
  dots <- list(...)
  if("main" %in% names(dots)) {
    themain <- dots$main
    usemaindefault <- FALSE
  }
  else usemaindefault <- TRUE

  if("xlim" %in% names(dots)) {
    xlim <- dots$xlim
    usexlimdefault <- FALSE
  }
  else usexlimdefault <- TRUE
  
  if("ylim" %in% names(dots)) {
    ylim <- dots$ylim
    useylimdefault <- FALSE
  }
  else useylimdefault <- TRUE

  if("xlab" %in% names(dots)) 
    xlab <- dots$xlab
  else {
    if(horizontal)
      xlab <- "Location (cM)"
    else
      xlab <- "Chromosome"
  }

  if("ylab" %in% names(dots)) 
    ylab <- dots$ylab
  else {
    if(horizontal)
      ylab <- "Chromosome"
    else
      ylab <- "Location (cM)"
  }

  map <- x
  # figure out if the input is a cross (containing a map)
  #    or is the map itself
  if(any(class(map) == "cross"))
    map <- pull.map(map)
  if(!missing(map2) && any(class(map2) == "cross"))
    map2 <- pull.map(map2)
  
  if(!any(class(map)=="map")  || (!missing(map2) && !any(class(map2) == "map")))
    warning("Input should have class \"cross\" or \"map\".")

  if(!missing(map2) && is.matrix(map[[1]]) != is.matrix(map2[[1]]))
      stop("Maps must be both sex-specific or neither sex-specific.")

  if(!missing(chr)) {
    map <- map[qtl:::matchchr(chr, names(map))]
    if(!missing(map2)) map2 <- map2[qtl:::matchchr(chr, names(map2))]
  }

  sex.sp <- FALSE

  if(is.matrix(map[[1]])) { # sex-specific map
    one.map <- FALSE
    sex.sp <- TRUE
    if(!missing(map2)) {
      if(is.logical(map2)) {
        horizontal <- map2
        map2 <- lapply(map,function(a) a[2,])
        map <- lapply(map,function(a) a[1,])
      }
      else {
        Map1 <- lapply(map,function(a) a[1,,drop=TRUE])
        Map2 <- lapply(map,function(a) a[2,,drop=TRUE])
        Map3 <- lapply(map2,function(a) a[1,,drop=TRUE])
        Map4 <- lapply(map2,function(a) a[2,,drop=TRUE])
        old.mfrow <- par("mfrow")
        on.exit(par(mfrow=old.mfrow))

        par(mfrow=c(2,1))
        class(Map1) <- class(Map2) <- class(Map3) <- class(Map4) <- "map"
        plotMap(Map1,Map3,horizontal=horizontal,shift=shift,
                 show.marker.names=show.marker.names,alternate.chrid=alternate.chrid)
        plotMap(Map2,Map4,horizontal=horizontal,shift=shift,
                 show.marker.names=show.marker.names,alternate.chrid=alternate.chrid)
        return(invisible(NULL))
      }
    }
    else {
      map2 <- lapply(map,function(a) a[2,])
      map <- lapply(map,function(a) a[1,])
    }
  }
  else { # single map
    # determine whether a second map was given
    if(!missing(map2)) 
      one.map <- FALSE
    else one.map <- TRUE
  }
       
  if(one.map) {

    n.chr <- length(map)
    if(!show.marker.names) {  # locations of chromosomes
      chrpos <- 1:n.chr
      thelim <- range(chrpos)+c(-0.5, 0.5)
    }
    else {
      chrpos <- seq(1, n.chr*2, by=2)
      thelim <- range(chrpos)+c(-0.35, 2.35)
    }

    if(shift) map <- lapply(map, function(a) a-a[1])
    maxlen <- max(unlist(lapply(map,max)))

    if(horizontal) {
      old.xpd <- par("xpd")
      old.las <- par("las")
      par(xpd=TRUE, las=1)
      on.exit(par(xpd=old.xpd,las=old.las))
      
      if(usexlimdefault) xlim <- c(0,maxlen)
      if(useylimdefault) ylim <- rev(thelim)
      plot(0,0,type="n",xlim=xlim, ylim=ylim,yaxs="i",
           xlab="", ylab="", yaxt="n", xaxt="n")
      title(xlab=xlab, mgp=c(1.5, 0, 0))
      title(ylab=ylab, mgp=c(2.1, 0, 0))
      a <- par("usr")
      rect(a[1], a[3], a[2], a[4], col="gray80")
      abline(v=1:19, col=darkgray, xpd=FALSE)
      abline(h=seq(0, 100, by=20), col="white", xpd=FALSE)
      rect(a[1], a[3], a[2], a[4], col=NA)
      axis(side=2, tick=FALSE, mgp=c(3, 0.2, 0))
      
      for(i in 1:n.chr) {
	segments(min(map[[i]]), chrpos[i], max(map[[i]]), chrpos[i])
        segments(map[[i]], chrpos[i]-tickwidth, map[[i]], chrpos[i]+tickwidth)

        if(show.marker.names)
          text(map[[i]], chrpos[i]+0.35, names(map[[i]]), srt=90, adj=c(1,0.5))
      }

      # add chromosome labels
      if(!alternate.chrid || length(chrpos) < 2) {
        for(i in seq(along=chrpos))
          axis(side=2, at=chrpos[i], labels=names(map)[i], tick=FALSE, mgp=c(3, 0.1, 0))
      }
      else {
        odd <- seq(1, length(chrpos), by=2)
        even <- seq(2, length(chrpos), by=2)
        for(i in odd) {
          axis(side=2, at=chrpos[i], labels="")
          axis(side=2, at=chrpos[i], labels=names(map)[i], line=-0.4, tick=FALSE, mgp=c(3, 0.1, 0))
        }
        for(i in even) {
          axis(side=2, at=chrpos[i], labels="")
          axis(side=2, at=chrpos[i], labels=names(map)[i], line=+0.4, tick=FALSE, mgp=c(3, 0.1, 0))
        }
      }
    }
    else {
      old.xpd <- par("xpd")
      old.las <- par("las")
      par(xpd=TRUE,las=1)
      on.exit(par(xpd=old.xpd,las=old.las))
      
      if(usexlimdefault) xlim <- thelim
      if(useylimdefault) ylim <- c(maxlen+1, -1)
      plot(0,0,type="n",ylim=ylim,xlim=xlim,xaxs="i",
           xlab="", ylab="", yaxt="n", xaxt="n", yaxs="i")
      title(xlab=xlab, mgp=c(1.5, 0, 0))
      title(ylab=ylab, mgp=c(2.1, 0, 0))
      
      a <- par("usr")
      rect(a[1], a[3], a[2], a[4], col="gray80")
      abline(v=1:19, col=darkgray, xpd=FALSE)
      abline(h=seq(0, 100, by=20), col="white", xpd=FALSE)
      rect(a[1], a[3], a[2], a[4], col=NA)
      axis(side=2, tick=FALSE, mgp=c(3, 0.2, 0))
      
      for(i in 1:n.chr) {
        segments(chrpos[i], min(map[[i]]), chrpos[i], max(map[[i]]))
        segments(chrpos[i]-tickwidth, map[[i]], chrpos[i]+tickwidth, map[[i]])

        if(show.marker.names)
          text(chrpos[i]+0.35, map[[i]], names(map[[i]]), adj=c(0,0.5))

      }
      # add chromosome labels
      if(!alternate.chrid || length(chrpos) < 2) {
        for(i in seq(along=chrpos))
          axis(side=1, at=chrpos[i], labels=names(map)[i], tick=FALSE, mgp=c(3, 0.1, 0))
      }
      else {
        odd <- seq(1, length(chrpos), by=2)
        even <- seq(2, length(chrpos), by=2)
        for(i in odd) {
          axis(side=1, at=chrpos[i], labels="")
          axis(side=1, at=chrpos[i], labels=names(map)[i], line=-0.4, tick=FALSE, mgp=c(3, 0.1, 0))
        }
        for(i in even) {
          axis(side=1, at=chrpos[i], labels="")
          axis(side=1, at=chrpos[i], labels=names(map)[i], line=+0.4, tick=FALSE, mgp=c(3, 0.1, 0))
        }
      }
    }
    if(usemaindefault)
      title(main="Genetic map")
    else if(themain != "")
      title(main=themain)
      
  }
  else {
    map1 <- map

    # check that maps conform
    if(is.matrix(map2[[1]]))
      stop("Second map appears to be a sex-specific map.")
    if(length(map1) != length(map2))
      stop("Maps have different numbers of chromosomes.")
    if(any(names(map1) != names(map2))) {
      cat("Map1: ", names(map1), "\n")
      cat("Map2: ", names(map2), "\n")
      stop("Maps have different chromosome names.")
    }

    if(shift) {
      map1 <- lapply(map1,function(a) a-a[1])
      map2 <- lapply(map2,function(a) a-a[1])
    }

    n.mar1 <- sapply(map1, length)
    n.mar2 <- sapply(map2, length)
    markernames1 <- lapply(map1, names)
    markernames2 <- lapply(map2, names)
    if(any(n.mar1 != n.mar2)) {
      if(show.marker.names) {
        warning("Can't show marker names because of different numbers of markers.")
        show.marker.names <- FALSE
      }
    }
    else if(any(unlist(markernames1) != unlist(markernames2))) {
      if(show.marker.names) {
        warning("Can't show marker names because markers in different orders.")
        show.marker.names <- FALSE
      }
    }

    n.chr <- length(map1)
    maxloc <- max(c(unlist(lapply(map1,max)),unlist(lapply(map2,max))))

    if(!show.marker.names) {  # locations of chromosomes
      chrpos <- 1:n.chr
      thelim <- range(chrpos)+c(-0.5, 0.5)
    }
    else {
      chrpos <- seq(1, n.chr*2, by=2)
      thelim <- range(chrpos)+c(-0.4, 2.4)
    }

    if(!horizontal) {
      old.xpd <- par("xpd")
      old.las <- par("las")
      par(xpd=TRUE,las=1)
      on.exit(par(xpd=old.xpd,las=old.las))

      if(usexlimdefault) xlim <- thelim
      if(useylimdefault) ylim <- c(maxloc, 0)
      plot(0,0,type="n",ylim=ylim,xlim=xlim, xaxs="i",
           xlab="", ylab="", yaxt="n", xaxt="n")
      title(xlab=xlab, mgp=c(1.5, 0, 0))
      title(ylab=ylab, mgp=c(2.1, 0, 0))

      a <- par("usr")
      rect(a[1], a[3], a[2], a[4], col="gray80")
      abline(v=1:19, col=darkgray, xpd=FALSE)
      abline(h=seq(0, 100, by=20), col="white", xpd=FALSE)
      rect(a[1], a[3], a[2], a[4], col=NA)
      axis(side=2, tick=FALSE, mgp=c(3, 0.2, 0))
    
      for(i in 1:n.chr) {
        if(max(map2[[i]]) < max(map1[[i]])) 
          map2[[i]] <- map2[[i]] + (max(map1[[i]])-max(map2[[i]]))/2

        else 
          map1[[i]] <- map1[[i]] + (max(map2[[i]])-max(map1[[i]]))/2
        
        segments(chrpos[i]-0.3, min(map1[[i]]), chrpos[i]-0.3, max(map1[[i]]))
        segments(chrpos[i]+0.3, min(map2[[i]]), chrpos[i]+0.3, max(map2[[i]]))
        
        # lines between markers
        wh <- match(markernames1[[i]], markernames2[[i]])
        for(j in which(!is.na(wh)))
          segments(chrpos[i]-0.3, map1[[i]][j], chrpos[i]+0.3, map2[[i]][wh[j]])
        if(any(is.na(wh)))
          segments(chrpos[i]-0.4, map1[[i]][is.na(wh)], chrpos[i]-0.2, map1[[i]][is.na(wh)])
        wh <- match(markernames2[[i]], markernames1[[i]])
        if(any(is.na(wh)))
          segments(chrpos[i]+0.4, map2[[i]][is.na(wh)], chrpos[i]+0.2, map2[[i]][is.na(wh)])

        if(show.marker.names)
          text(chrpos[i]+0.35, map2[[i]], names(map2[[i]]), adj=c(0,0.5))

      }
      # add chromosome labels
      if(!alternate.chrid || length(chrpos) < 2) {
        for(i in seq(along=chrpos))
          axis(side=1, at=chrpos[i], labels=names(map1)[i], tick=FALSE, mgp=c(3, 0.1, 0))
      }
      else {
        odd <- seq(1, length(chrpos), by=2)
        even <- seq(2, length(chrpos), by=2)
        for(i in odd) {
          axis(side=1, at=chrpos[i], labels="")
          axis(side=1, at=chrpos[i], labels=names(map1)[i], line=-0.4, tick=FALSE, mgp=c(3, 0.1, 0))
        }
        for(i in even) {
          axis(side=1, at=chrpos[i], labels="")
          axis(side=1, at=chrpos[i], labels=names(map1)[i], line=+0.4, tick=FALSE, mgp=c(3, 0.1, 0))
        }
      }
    }
    else {
      old.xpd <- par("xpd")
      old.las <- par("las")
      par(xpd=TRUE,las=1)
      on.exit(par(xpd=old.xpd,las=old.las))

      if(usexlimdefault) xlim <- c(0,maxloc)
      if(useylimdefault) ylim <- rev(thelim)
      plot(0,0,type="n",xlim=xlim,ylim=ylim, xaxt="n",
           xlab="", ylab="", yaxt="n", yaxs="i")
      title(xlab=xlab, mgp=c(1.5, 0, 0))
      title(ylab=ylab, mgp=c(2.1, 0, 0))

      a <- par("usr")
      rect(a[1], a[3], a[2], a[4], col="gray80")
      abline(v=1:19, col=darkgray, xpd=FALSE)
      abline(h=seq(0, 100, by=20), col="white", xpd=FALSE)
      rect(a[1], a[3], a[2], a[4], col=NA)
      axis(side=2, tick=FALSE, mgp=c(3, 0.2, 0))
    
      for(i in 1:n.chr) {
      
        if(max(map2[[i]]) < max(map1[[i]])) 
          map2[[i]] <- map2[[i]] + (max(map1[[i]])-max(map2[[i]]))/2
        else 
          map1[[i]] <- map1[[i]] + (max(map2[[i]])-max(map1[[i]]))/2
        
        segments(min(map1[[i]]), chrpos[i]-0.3, max(map1[[i]]), chrpos[[i]]-0.3)
        segments(min(map2[[i]]), chrpos[i]+0.3, max(map2[[i]]), chrpos[[i]]+0.3)

        # lines between markers
        wh <- match(markernames1[[i]], markernames2[[i]])
        for(j in which(!is.na(wh)))
          segments(map1[[i]][j], chrpos[i]-0.3, map2[[i]][wh[j]], chrpos[i]+0.3)
        if(any(is.na(wh)))
          segments(map1[[i]][is.na(wh)], chrpos[i]-0.4, map1[[i]][is.na(wh)], chrpos[i]-0.2)
        wh <- match(markernames2[[i]], markernames1[[i]])
        if(any(is.na(wh)))
          segments(map2[[i]][is.na(wh)], chrpos[i]+0.4, map2[[i]][is.na(wh)], chrpos[i]+0.2)

        if(show.marker.names)
          text(map2[[i]], chrpos[i]+0.35, names(map2[[i]]), srt=90, adj=c(1,0.5))
          
      }
      # add chromosome labels
      if(!alternate.chrid || length(chrpos) < 2) {
        for(i in seq(along=chrpos))
          axis(side=2, at=chrpos[i], labels=names(map1)[i], tick=FALSE, mgp=c(3, 0.1, 0))
      }
      else {
        odd <- seq(1, length(chrpos), by=2)
        even <- seq(2, length(chrpos), by=2)
        for(i in odd) {
#          axis(side=2, at=chrpos[i], labels="")
          axis(side=2, at=chrpos[i], labels=names(map1)[i], line=-0.4, tick=FALSE, mgp=c(3, 0.1, 0))
        }
        for(i in even) {
#          axis(side=2, at=chrpos[i], labels="")
          axis(side=2, at=chrpos[i], labels=names(map1)[i], line=+0.4, tick=FALSE, mgp=c(3, 0.1, 0))
        }
      }

    }
    if(usemaindefault) {
      if(!sex.sp) title(main="Comparison of genetic maps")
      else title(main="Genetic map")
    }
    else if(themain != "")
      title(main=themain)
  }    
  invisible()
}
