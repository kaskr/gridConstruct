## Create simple sparse model.matrix from a factor
modMat <- function(fac,x=rep(1,length(fac))){
  if(is(fac,"sparseMatrix"))return(fac)
  Dim <- c(length(fac),nlevels(fac))
  j <- as.integer(as.integer(fac))
  i <- as.integer(1:Dim[1])
  keep <- !is.na(j) ## Insert zero if a value is missing
  ans <- spMatrix(Dim[1],Dim[2],i[keep],j[keep],x[keep])
  as(ans,"dgCMatrix")
}

## Image plotting of abundance surface.
## Transform log-surface to contour levels of abundance suface.  x is
## mapped to [0,1] such that equidistant subdivision of [0,1]
## corresponds to D(10%), D(20%) etc
## x0 -> sum(exp(x)*(x<x0))/sum(exp(x)) = sum(p*(x<x0)) where
## p=exp(x)/sum(exp(x))
concTransform <- function(x){
  i <- order(x)
  ys <- sort(exp(x))
  p <- ys/sum(ys)
  x[i] <- cumsum(p)
  x
}
colorstrip <- function(col){
  ousr <- par("usr")
  on.exit(par(usr=ousr))
  par(usr=c(-40,1,0,1))
  levels <- seq(0,1,length=length(col)+1)
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col,border=NA)
}

## FIXME: Remove probably
## Utility to construct "distFactor"
distFactor <- function(...,distfun=match.fun("-"),sign=TRUE){
  sc <- sys.call()
  nm <- as.character(sc)[-1]
  dotargs <- list(...)
  names(dotargs) <- nm[seq(dotargs)]
  dotargs <- lapply(dotargs,as.factor)
  levs <- lapply(dotargs,levels)
  nlevs <- sapply(levs,length)
  e <- rev(expand.grid(rev(lapply(levs,seq))))
  chr2num <- function(x)sapply(strsplit(x,"[^0-9\\.-]+"),function(x)mean(as.numeric(x),na.rm=TRUE))
  fac <- interaction(dotargs,sep=":",lex.order=TRUE)
  var <- lapply(levs,chr2num)
  mat <- lapply(var,function(x)outer(x,x,distfun))
  if(!sign)mat <- lapply(mat,abs)
  dist <- Map(function(x,i)x[i,i],mat,e)
  names(dist) <- paste("d",names(dist),sep="")
  dist$I <- diag(prod(nlevs))
  attr(dist,"Dim") <- rev(nlevs)
  levelAttrib(fac)$dist <- dist
  class(fac) <- c("distFactor",class(fac))
  fac
}
print.distFactor <- function(x,...){
  cat("Factor of length",length(x),"with",nlevels(x),"levels\n\n")
  cat("Attached distance matrices\n")
  print(lapply(levelAttrib(x)$dist,dim))
  cat("Grid dimension\n")
  print(attr(levelAttrib(x)$dist,"Dim"))
}
"[.distFactor" <- function(x,...){
  lattrx <- levelAttrib(x)
  y <- NextMethod("[")
  levelAttrib(y) <- lattrx
  class(y) <- oldClass(x)
  y
}
