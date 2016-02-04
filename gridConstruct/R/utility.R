## Observation number within group:
## Assign to each observation its number within the group-level that it belongs
obsnum <- function(fac){
  ans <- numeric(n <- length(fac))
  for(i in split(1:n,fac))ans[i] <- 1:length(i)
  ans
}

## Create simple sparse model.matrix from a factor
## modMat <- function(fac,x=rep(1,length(fac))){
##   Dim <- as.integer(c(length(fac),nlevels(fac)))
##   j <- as.integer(as.integer(fac)-1)
##   i <- as.integer((1:Dim[1])-1)
##   x <- as.numeric(x)
##   as(new("dgTMatrix",i=i,j=j,x=x,Dim=Dim),"dgCMatrix")
## }
## --- New
modMat <- function(fac,x=rep(1,length(fac))){
  if(is(fac,"sparseMatrix"))return(fac)
  Dim <- c(length(fac),nlevels(fac))
  j <- as.integer(as.integer(fac))
  i <- as.integer(1:Dim[1])
  keep <- !is.na(j) ## Insert zero if a value is missing
  ans <- spMatrix(Dim[1],Dim[2],i[keep],j[keep],x[keep])
  as(ans,"dgCMatrix")
}


## ------ Shave a large fitted object ----------------
## shave <- function(x,mb=5){
##   size <- function(x){
##     ans <- try(length(serialize(x,NULL)))
##     if(inherits(ans,"try-error"))return(Inf)
##     ans
##   }
##   i <- sapply(x,size)/1024^2 < mb
##   ans <- x[i]
##   class(ans) <- class(x)
##   ans
## }
## Reduce the very large fitted object
shave <- function(obj){
  ans <- obj[c("par","hessian","data","vcov","sd","eta","Dim","hessian.full","fixefnames","call","value")]
  class(ans) <- class(obj)
  ans
}

## General function for object returned by "optim".
## Computes simulated confidence bands of any vector-function
simconfint <- function(fit,
                       fun=function(par)par,
                       n=1000,
                       reduce=function(x)quantile(x,c(.025,.5,.975)),
                       ...){
  ## par <- GMRFsample(fit$hessian,n)+fit$par
  ## nam <- names(fit$par)
  par <- GMRFsample(solve(forceSymmetric(vcov(fit))),n)+coef(fit)
  nam <- names(coef(fit))
  sim <- apply(par,2,function(x){names(x) <- nam;fun(x,...)})
  if(is.null(reduce))return(sim)
  if(is.vector(sim))sim <- t(sim)
  t(apply(sim,1,reduce))
}

##' Likelihood ratio tests of nested models.
##'
##' @param ... Nested models to test.
##' @param main Test all models against main instead of succesive testing?
##' @param latex Generate latex output?
##' @return matrix
LRtest <- function(...,main=NULL,latex=FALSE){
  suc <- is.null(main)  #Succesive testing?
  args <- list(...)
  if(is.null(names(args)))names(args) <- sapply(substitute(f(...)),deparse)[-1]
  nparms <- sapply(args,function(x)length(x$par))
  lik <- sapply(args,function(x)sum(x$value))

  ## ----- For fitlgc class
  if("fitlgc" %in% class(args[[1]])){
    nparms <- nparms + sapply(args,function(x)x$nfixed)
    lik <- -sapply(args,logLik)    
  }
  ## -----------------------
  
  casenames <- names(args)
  if(suc)minus2logQ <- c(NA,2*diff(lik)) else minus2logQ <- 2*(lik-main$value)
  if(suc)df <- c(NA, -diff(nparms)) else df <- length(main$par)-nparms
  p <- 1-pchisq(minus2logQ,df=df)
  m <- cbind(lik,minus2logQ,nparms,df,p)
  rownames(m) <- casenames
  if(latex)colnames(m)[1:2] <- c('$-\\log L$','$-2\\log Q$')
  else colnames(m)[1:2] <- c('-log L','-2log Q')
  m
}

## Image plotting of abundance surface.
## Transform log-surface to contour levels of abundance suface.
## x is mapped to [0,1] such that equidistant subdivision of [0,1] corresponds to D(10%), D(20%) etc
## x0 -> sum(exp(x)*(x<x0))/sum(exp(x))  =  sum(p*(x<x0)) where p=exp(x)/sum(exp(x))
concTransform <- function(x){
  i <- order(x)
  ys <- sort(exp(x))
  p <- ys/sum(ys)
  x[i] <- cumsum(p)
  x
}
colorstrip <- function(col){
  ##opar <- par()
  ousr <- par("usr")
  on.exit(par(usr=ousr))
  par(usr=c(-40,1,0,1))
  levels <- seq(0,1,length=length(col)+1)
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col,border=NA)
  ##par(opar)
}
 

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

## Vector containing NAs --> sparse vector with NA replaced by zero.
na2sparse <- function(v) {
  i <- which(is.finite(v))
  as(spMatrix(length(v), 1, i = i, j = rep(1, length(i)), 
              v[i]), "dgCMatrix")
}

chol.solve <- function (a, b, ...){
  u <- chol(a) ## a=t(u)%*%u
  if(!missing(b)){
    y <- backsolve(u,b,transpose=TRUE)
    x <- backsolve(u,y)
    x
  } else {
    x <- chol2inv(u)
    attr(x,"logdet") <- 2*sum(log(diag(u)))
    x
  }
}
## Special function to combine huge lists
## Overloaded by lgcopt
c.list <- .Primitive("c")

## LAPPLY and REMOVE - USE WITH CARE!
## Used some places to avoid too many copies of very large objects.
rmLapply <- function(env,x,FUN,...){
  ans <- lapply(get(x,envir=env),FUN,...)
  rm(list=x,envir=env)
  ans
}
