##' \code{gridConstruct} constructs a grid.
##'
##' These functions help to construct the object required to
##' build Gaussian Markov Random fields with the formula interface.
##' Most situations can be handled by
##' \enumerate{
##' \item Creating a grid using \code{grid <- gridConstruct(data)}.
##' \item Building a gridFactor object using \code{gridFactor(data,grid)}.
##' }
##' 
##' Construction of grids in practice generally involves 3 steps:
##' \enumerate{
##' \item \code{gridConstruct} - Construct the grid sufficiently fine and
##' sufficiently large to contain all data points.
##' \item \code{gridFilter} - Filter off unwanted grid points. For instance
##' grid points on land or grid points too far away from the region of
##' interest.
##' \item \code{gridLocate} - For each data point locate the nearest grid
##' point.
##' }
##' 
##' @title Construct a grid based on observation locations.
##' @param data data.frame.
##' @param type Neighborhood structure type.
##' @param filter Call \code{gridFilter} after grid construction?
##' @param ...
##' @param center Optional list to control origo of grid.
##' @param km Optional distance between neighboring grid points in km.
##' @return grid object
##' @rdname grid
##' @example demo/gridConstruction.R
gridConstruct <- function(data,type=
                          c("squareGrid",
                            "triangularGrid",
                            "scatterGrid"),
                          filter=TRUE,
                          wet=!wetEdges,
                          wetEdges=FALSE,
                          connected=TRUE,
                          ordertol=1,
                          ...){                              
  type <- match.arg(type)
  constructor <- match.fun(type)
  ans <- constructor(data,...)
  if(type=="scatterGrid")filter <- FALSE
  if(filter)ans <- gridFilter(grid=ans,data=data,wet=wet,
                              connected=connected,ordertol=ordertol,
                              wetEdges=wetEdges,...)
  ans
}
##' \code{gridFilter} filters off unwanted grid points.
##'
##' @title Filter off unwanted grid points.
##' @param icesSquare Remove grid points outside ICES squares in the data?
##' @param nearestObs Remove grid points with closest data point greater than \code{nearestObs}.
##' @param wet Remove grid points on land.
##' @param wetEdges Alternative: Keep edges passing through water. 
##' @param ordertol Require at least \code{ordertol} neighbors to every grid point.
##' @param connected Keep only largest connected component.
##' @return Filtered grid object
##' @rdname grid
gridFilter <- function(grid,data,
                       icesSquare=FALSE,
                       nearestObs=Inf,
                       wet=FALSE,
                       connected=FALSE,
                       ordertol=0,
                       wetEdges=FALSE,
                       ...){
  keep <- rep(TRUE,nrow(grid))
  ## icesSquare filtering
  if(icesSquare){
    require(DATRAS)
    sq <- unique(DATRAS::icesSquare(data))
    ##pol <- icesSquare2coord(sq,format="polygons")
    ##keep <- keep & !is.na(cut(grid,pol))
    sqg <- DATRAS::icesSquare(grid)
    keep <- keep & (sqg %in% sq)
  }
  ## nearest observation filtering
  if(is.finite(nearestObs)){
    ## i <- gridLocate(data,grid)
    ## d <- dist.km(grid,data[i,],outer=FALSE)
    ## keep <- keep & (d<=nearestObs)
    i <- gridLocate(data,grid[keep,])
    d <- dist.km(grid[keep,],data[i,],outer=FALSE)
    keep[keep] <- (d<=nearestObs)
  }
  ## map filtering
  if(wet){
    ##keep <- keep & is.na(cut(grid))
    keep[keep] <- is.na(cut(grid[keep,]))
  }
  if(wetEdges){
    ## Reduce pattern matrix first
    Q <- attr(grid,"pattern")
    Q[!keep,] <- 0
    Q[,!keep] <- 0
    attr(grid,"pattern") <- Q
    ## Now "wetEdgeFilter" has less work
    grid <- wetEdgeFilter(grid,...)
  }
  grid <- grid[keep,]
  ## Remove gridpoints with too few neighbors
  if(ordertol>0){
    rem <- function(i0) {
      p0 <- attr(grid, "pattern")[i0, i0]
      cs <- colSums(p0) - 1
      keep <- cs > ordertol
      i0 <- i0[keep]
      if (any(!keep)) 
        i0 <- rem(i0)
      i0
    }
    i0 <- rem(seq(length=nrow(grid)))
    grid <- grid[i0,]
  }
  ## Keep largest connected component
  if(connected){
    c <- connectedComponents(grid)
    if(length(c)>1){
      c <- c[[which.max(sapply(c,length))]]
      grid <- grid[c,]
    }
  }
  grid
}

## ==============================================================
## Binary search of closest gridpoint.
## ==============================================================
gridLocateBin <- function(grid,points){
  recursiveLocate <- function(x,y){
    ##print(nrow(x))
    ##if(nrow(x)==1)return(x)
    if(nrow(x)==0 | nrow(y)==0)return(NULL)
    if(nrow(x)==1)return(rep(x$key,nrow(y)))
    x[1:2] <- x[2:1]; y[1:2] <- y[2:1]
    med <- median(x[,1])
    i <- factor(x[,1]<med,levels=c(TRUE,FALSE))
    spx <- split(x,i)
    j <- factor(y[,1]<med,levels=c(TRUE,FALSE))
    spy <- split(y,j)
    rm(x,y)
    ans <- mapply(recursiveLocate,spx,spy,SIMPLIFY=FALSE)
    if(any(sapply(ans,length)==0))return(unlist(ans))
    else return(unsplit(ans,j))
  }
  grid <- as.data.frame(grid[c("lon","lat")])
  points <- as.data.frame(points[c("lon","lat")])
  grid$key <- 1:nrow(grid)
  recursiveLocate(grid,points)
}

## ==============================================================
## Brute force search of closest gridpoint.
## ==============================================================
##' Locate nearest grid point to a set of locations.
##'
##' \code{gridLocate} finds the closest grid point for a set of locations.
##' 
##' \code{gridLocate} performs a brute force search of closest gridpoint
##' to each data point. The index of the closest grid point is returned.
##' @param grid 
##' @param points 
##' @rdname grid
gridLocate <- function(grid,points){
  nearestkm(points$lon,points$lat,grid$lon,grid$lat)
}


##' \code{gridFactor} constructs a gridFactor object.
##' @rdname grid
gridFactor <- function(data,grid,...){
  data <- data[c("lon","lat")]
  ch <- paste(data$lon,data$lat)
  nd <- which(!duplicated(ch))
  fac <- unclass(factor(ch,levels=ch[nd]))
  i <- gridLocate(grid,data[nd,])
  fac <- factor(i[fac],1:nrow(grid))
  levelAttrib(fac)$grid <- grid
  class(fac) <- c("gridFactor", class(fac))
  return(fac)
}

## Calling "factor(gridFactor)" should not remove levels
## Workaround: Overload factor function (because "factor" is not generic)
## 
## factor <- function(x,...){
##   if(inherits(x,"gridFactor"))return(x)
##   if(inherits(x,"distFactor"))return(x)
##   base::factor(x,...)
## }

## NOT WORKING: because:
## - other packages like DATRAS ignores the overloaded method
## - when e.g. PROJ(GMRF(gf)) tries to factor out un-used levels it is ignored!
## setGeneric("factor")
## setMethod("factor","gridFactor",function (x = character(), levels, labels = levels, exclude = NA, 
##     ordered = is.ordered(x)){warning("gridFactor not modified by factor()");x})


## scatterGrid. Used to build covariance models directly
scatterGrid <- function(data,sample.id){
  if(missing(sample.id))stop("scatterGrid requires a sample.id")
  id <- data[[sample.id]]
  if(!is.factor(id))stop("sample.id must be the name of a factor in the data.")
  i <- tapply(seq(length=nrow(data)),id,function(x)x[1])
  data[i,,drop=FALSE]
}

wetEdgeFilter <- function(grid,wettol=0.1,edgeDiscretize=11,...){
  Q <- attr(grid,"pattern")
  Q@x[Q@i<=Q@j] <- 0
  Q <- drop0(Q)
  Q <- as(Q,"dgTMatrix")
  i <- Q@i+1
  j <- Q@j+1
  gr1 <- as.data.frame(grid)[i,]
  gr2 <- as.data.frame(grid)[j,]
  f <- function(t){
    print(t)
    v <- gr2-gr1
    is.na(cut.polygonGrid(gr1+t*v))
  }
  cat("Edge discretize [0,1]:\n")
  mat <- sapply(seq(0,1,length=edgeDiscretize),f)
  wet <- rowMeans(mat)
  k <- wet>wettol
  Q@x[!k] <- 0
  Q <- drop0(Q)
  Q <- Q+t(Q)
  diag(Q) <- 1
  attr(grid,"pattern") <- as(Q,"dgTMatrix")
  grid
}


## -----------------------------------------------------------
## Make extension of class "factor" to have "level atributes"
## Ideas for level attributes:
## * distance matrix "dist"
## * Pattern matrix "pattern"
levelAttrib <- function(x){
  if(!is.factor(x))stop("x must be a factor")
  attr(x,"levelAttrib")
}
"levelAttrib<-" <- function(x,value){
  if(!is.factor(x))stop("x must be a factor")  
  attr(x,"levelAttrib") <- value
  x
}

"[.gridFactor" <- function(x,...){
  lattrx <- levelAttrib(x)
  y <- NextMethod("[")
  levelAttrib(y) <- lattrx
  class(y) <- oldClass(x)
  y
}
## Special case: wrap levels og a factor around a circle
## * Call factor(...)
## * Set levelAtrrib
circFactor <- function(...){
  fac <- factor(...)
  gr <- data.frame(x=levels(fac))
  n <- nlevels(fac)
  i <- 1:n
  P <- spMatrix(n,n,i=i,j=c(i[-1],i[1]),x=rep(1,n))
  P <- P+(t(P))+.symDiagonal(n)
  attr(gr, "pattern") <- as(P,"dgCMatrix")
  levelAttrib(fac)$grid <- gr
  fac
}
