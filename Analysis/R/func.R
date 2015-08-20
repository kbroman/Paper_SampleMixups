myknn <-
function(x, y, k=40, l=35)
{
  if(!is.matrix(y)) y <- as.matrix(y)
  noy <- apply(y, 1, function(a) any(is.na(a)))
  nox <- is.na(x)

  ytr <- y[(!nox & !noy),,drop=FALSE]
  xtr <- x[(!nox & !noy)]
  ytest <- y[!noy,,drop=FALSE]

  require(class)
  knn.out <- knn(train=ytr, test=ytest, cl=xtr, k=k, l=l)
  cat("predictions: ", sum(!is.na(knn.out)), "/", length(knn.out), "\n")

  temp <- rep(NA, length(x))
  temp[!noy] <- knn.out
  if(is.factor(knn.out)) {
    temp <- as.factor(temp)
    levels(temp) <- levels(knn.out)
  }
  data.frame(y, x, knn.out=temp)
}
  

######################################################################

count.qtl <-
function(scanone_result, theprobepos=probepos, threshold=5, distance=2.5)
{
  mx <- apply(scanone_result[,-(1:2)], 2, tapply, scanone_result[,1], max)
  neqtl <- colSums(mx >= 5)

  theprobepos <- theprobepos[names(out)[-(1:2)],]

  chr <- as.character(theprobepos[,1])
  pos <- theprobepos[,3]

  ncis <- rep(0, length(neqtl))
  names(ncis) <- names(neqtl)
  for(i in seq(along=cis)) {
    if(mx[chr[i], i] >= 5) {
      if(!is.na(distance)) {
        tmp <- scanone_result[scanone_result[,1]==chr[i], c(2,i+2)]
        tmp <- mean(tmp[tmp[,2]==max(tmp[,2]),1])
        if(abs(tmp-pos[i]) <= distance) ncis[i] <- 1
      } else {
        li <- lodint(scanone_result[,c(1,2,i+2)], chr=chr[i])[,2]
        if(li[1] <= pos[i] && pos[i] <= li[length(li)])
          ncis[i] <- 1
      }
    }
  }
  cbind(ncis=ncis, ntrans=neqtl-ncis)
}
