## ==============================================================
## Geographic coordinates <-> Cartesian coordinates.
## Found on wikipedia.
## ==============================================================
sph2car <- function(obj){
  G2R <- pi/180
  lat <- obj$lat*G2R; lon <- obj$lon*G2R
  data.frame(x=cos(lat)*cos(lon),
             y=cos(lat)*sin(lon),
             z=sin(lat))
}
car2sph <- function(obj){
  R2G <- 180/pi
  local( data.frame(lon=R2G*atan2(y,x),
                    lat=R2G*atan2(z,sqrt(x*x+y*y))) , obj)
}
sphMidpoint <- function(obj1,obj2){
  carmidpoint <- .5*(sph2car(obj1)+sph2car(obj2))
  car2sph(carmidpoint)
}
carDist <- function(obj1,obj2){
  sqrt(rowSums( data.frame(mapply("-",obj1,obj2,SIMPLIFY=FALSE)) ^2))
}
sphDist <- function(obj1,obj2,radius=6378.1){
  dcar <- carDist(sph2car(obj1),sph2car(obj2))
  dsph <- 2*asin(.5*dcar) ## distance on unit-sphere
  radius*dsph
}

## ==============================================================
## Rotation around x,y or z- axis
## (default y-axis  corresponding to North/South movement)
## Example: Northpole (lon=0,lat=90) is mapped to (lon=0,lat=0)
##          if angle=90.
## ==============================================================
sphRot <- function(obj,angle=0,axis=2){
  axis <- match(as.character(axis),1:3)
  G2R <- pi/180; angle <- G2R*angle
  car <- sph2car(obj)
  ## Example: axis=2  =====>   
  ## Q=
  ##   |  cos(angle) ,   0,   sin(angle) |
  ##   |           0 ,   1,           0  |
  ##   |  -sin(angle),   0,   cos(angle) |
  u <- axis == (1:3) ## Unit vector to rotate around
  Q <- tcrossprod(u)
  Q[!u,!u] <- c( cos(angle), -sin(angle), sin(angle), cos(angle) )
  ans <- Q %*% t(car)
  rownames(ans) <- c("x","y","z")
  ans <- car2sph( as.data.frame(t(ans)) )
  attributes(ans) <- attributes(obj)
  ans
}

## ==============================================================
## Move an object centered at Northpole by first rotating the
## object "angle" degrees around z-axis followed by moving the
## center to the spherical coordinate "final".
## ==============================================================
sphMoveNorthpole <- function(obj,final,angle=0){
  obj$lon <- obj$lon+angle
  roty <- 90-final$lat
  ans <- sphRot(obj,roty)
  ans$lon <- ans$lon+final$lon
  ans
}

## ==============================================================
## Given coordinates (lon,lat).
## Compute circle surrounding the points.
## ==============================================================
surroundCircle <- function(data,scale.radius=1,
                           center = data.frame(lon=mean(range(data$lon)),lat=mean(range(data$lat))),
                           ...){
  angle <- 180/pi*max(2*asin(carDist(sph2car(center),sph2car(data))/2))
  angle <- scale.radius*angle
  list(center=center,angle=angle)
}

## ==============================================================
## Given coordinates (lon,lat).
## Compute regular triangulated hexagon surrounding the points.
## ==============================================================
initialHexagon <- function(data,...){
  sc <- surroundCircle(data,...)
  ## The circle centered at the North-pole with latitude=90-maxangle surrounds
  ## the point pattern.
  angle <- sc$angle; center <- sc$center
  maxangle <- angle*2/sqrt(3) ## Polygon must surround circle
  shape <- data.frame(lon=seq(-180,180,by=60),lat=90-maxangle)
  newshape <- sphMoveNorthpole(shape,center,00)
  pol0 <- lapply(c(lapply(1:5,"+",0:1),list(c(6,1))),function(i)rbind(newshape[i,],center))
  do.call("rbind",pol0)
}

## ==============================================================
## Given coordinates (lon,lat).
## Compute square surrounding the points centered at (0,0).
## ==============================================================
initialSquare <- function(data,...){
  sc <- surroundCircle(data,...)
  v <- sc$angle*c(-1,1)
  expand.grid(lon=v,lat=v)[c(1,2,4,3),]
}

## ==============================================================
## Triangle storage:
##    First three rows of data.frame ~ first triangle corners
##    Next three rows of data.frame ~ second triangle corners
##    ...
##
## Vectorized subdivision:
##
##
##   P1____M3___P3        M1          |
##     \  / \  /          M2          | 1st triangle
##      \/_ _\/           M3          |
##    M1 \   /M2          M1             |
##        \ /             M2 <- P1       | 2nd triangle
##         V              M3             |
##         P2             M1                |
##                        M2                | 3rd triangle
##                        M3 <- P2          |
##                        M1 <- P3            |
##                        M2                  | 4th triangle
##                        M3                  |
##
## Algorithm:
##    1. From P1,P2,P3 compute M1,M2,M3
##    2. Repeat M1,M2,M3 four times
##    3. Insert P1,P2,P3 on components 5,9,10
makeTriangles <- function(data,n=5,triangle0=initialHexagon(data,...),km,...){
  data <- as.data.frame(data[c("lon","lat")])
  
  if(!missing(km)){
    m <- max(dist.km(triangle0))/2
    n <- ceiling(log2(m/km))
    scale <- 2^n*km/m
    triangle0 <- initialHexagon(triangle0,scale=sqrt(3)/2*scale)
  }
  
  ## -------------- Perform subdivision
  midpoint <- function(d1,d2){ ## Spherical midpoint
    G2R <- pi/180;
    d1 <- d1*G2R; d2 <- d2*G2R
    lon1 <- d1$lon;lon2 <- d2$lon
    lat1 <- d1$lat;lat2 <- d2$lat
    a <- sin(lat1)+sin(lat2)
    b <- cos(lat1)*cos(lon1)+cos(lat2)*cos(lon2)
    c <- cos(lat1)*sin(lon1)+cos(lat2)*sin(lon2)
    x <- atan(c/b)
    y <- atan(a/(b*cos(x) + c*sin(x)))
    data.frame(lon=x/G2R,lat=y/G2R)
  }
  
  subdiv4 <- function(x){
    extend <- function(v,m,n){outer( v, (0:(n-1))*m, "+")}
    n <- nrow(x)/3 ## Number of triangles
    y <- x[extend( c(2,3,1), 3, n ), ]
    mp <- midpoint(x,y)
    ans <- mp[ extend(rep(1:3,4), 3, n ) ,]
    ans[ extend(c(5,9,10), 12, n) ,] <- x
    ans
  }

  subdiv <- function(x,n=3){ ## Subdivide recursively
    for(i in 1:n)x <- subdiv4(x)
    x
  }

  triangles <- subdiv(triangle0,n=n)
  
  ans <- triangles
  attr(ans,"data") <- data
  
  class(ans) <- c("triangles","polygons",class(ans))
  ans
}

## ==============================================================
## Square subdivision:
##
##   P1____M4___P4        M5
##    |    |    |         M5 <- M1
##    |    |    |         M5 <- P1
##  M1|----M5---|M3       M5 <- M4
##    |    |    |         M5
##    |____|____|         M5 <- M4
##   P2    M2   P3        M5 <- P4
##                        M5 <- M3
##                        M5
##                        M5 <- M3
##                        M5 <- P3
##                        M5 <- M2
##                        M5
##                        M5 <- M2
##                        M5 <- P2
##                        M5 <- M1
## M: 2,16,12,14,8,10,4,6
## P: 3,15,11,7
makeSquares <- function(data,n=5,square0=initialSquare(data,...),
                        center=data.frame(lon=mean(range(data$lon)), lat=mean(range(data$lat))),...){
  center <- as.list(center)
  midpoint <- function(d1,d2){ 
    .5*(d1+d2)
  }

  subdiv4 <- function(x){
    extend <- function(v,m,n){outer( v, (0:(n-1))*m, "+")}
    n <- nrow(x)/4 ## Number of squares
    y <- x[extend( c(2,3,4,1), 4, n ), ]
    mp <- midpoint(x,y)
    mp1 <- mp[extend(1,4,n),]
    mp3 <- mp[extend(3,4,n),]
    mp5 <- midpoint(mp1,mp3)
    ans <- mp5[extend(rep(1,16),1 ,n),]
    ans[ extend(c(2,16,12,14,8,10,4,6), 16, n), ] <- mp[extend(c(1,1,2,2,3,3,4,4), 4, n),]
    ans[ extend(c(3,15,11,7), 16, n), ] <- x
    ans
  }

  subdiv <- function(x,n=3){ ## Subdivide recursively
    for(i in 1:n)x <- subdiv4(x)
    x
  }

  squares <- subdiv(square0,n=n)

  ans <- sphRot(squares,-center$lat)
  ans$lon <- ans$lon+center$lon

  attr(ans,"data") <- data

  class(ans) <- c("squares","polygons",class(ans))
  ans
}

## ==============================================================
## Methods for class "polygons"
## ==============================================================
nsidesPolygons <- function(x){
  ch <- class(x)[1]
  switch(ch,
         "triangles"=3,
         "squares"=4,
         stop("Unknown polygon"))
}
splitPolygons <- function(x){
  ns <- nsidesPolygons(x)
  rbind(x,NA)[rbind(matrix(1:nrow(x),ns),nrow(x)+1),]
}
summary.polygons <- function(object,...){
  extend <- function(v,m,n){outer( v, (0:(n-1))*m, "+")}
  nsides <- nsidesPolygons(object)
  x <- object
  n <- nrow(x)/nsides ## Number of polygons
  y <- x[extend( c((1:nsides)[-1],1), nsides, n ), ]
  side.length <- t(matrix(sphDist(x,y),nsides))
  s <- .5*rowSums(side.length)
  if(nsides==3){
    area <- exp(.5*rowSums(log(cbind(s-side.length,s)))) ## Herons formula
  } else area <- rep(NA,n)
  list(n=n,
       side.length=side.length,
       area=area
       )
}
plot.polygons <- function(x,add=FALSE,...){
  data <- attr(x,"data")
  if(!add)plot(as.data.frame(lapply(x,range)),col="white")
  polygon(splitPolygons(x),...)
  points(data$lon,data$lat,col="red",pch=16,cex=.5)
}

heatpal <- colorRampPalette(c("white", "yellow","orange" ,"red"))
image.polygons <- function(x,y,add=FALSE,col=heatpal(16), ##col=heat.colors(16),
                           codes=cut(y,length(col)),border=FALSE,
                           strip=colorstrip(col),
                           map=polygon(makeMap(x,map="world")),
                           breaks=NULL,...){
  if(!is.null(breaks)){
    stopifnot(length(breaks)==length(col)+1)
    codes <- cut(y,breaks)
  }
  if(!add)plot(as.data.frame(lapply(x,range)),type="n")
  polygon(splitPolygons(x),col=as.character(col[codes]),border=border,...)
  eval(map)
  eval(strip)
}
print.polygons <- function(x,...){
  s <- summary(x)
  cat("Number of polygons:",s$n,"\n")
  cat("\nDistribution of side lengths:\n")
  print(quantile(s$side.length))
}

## ==============================================================
## Methods for class "triangles"
## ==============================================================
splitTriangles <- function(x){
  rbind(x,NA)[rbind(matrix(1:nrow(x),3),nrow(x)+1),]
}
print.triangles <- function(x,...){
  s <- summary(x)
  cat("Number of triangles:",s$n,"\n")
  cat("\nDistribution of side lengths:\n")
  print(quantile(s$side.length))
  cat("\nDistribution of triangle areas:\n")
  print(quantile(s$area))  
}

## ==============================================================
## Based on "triangles" object:
##     Compute grid
##     Compute neighborhood pattern
## ==============================================================
.incidenceMatrix <- function(fac,ns,n=length(fac)/ns){
  ## Incidence matrix of n independent ns-cycles
  c <- c((1:ns)[-1],1)
  C <- as(outer(1:ns,1:ns,function(x,y)(x-y)%%ns %in% c(0,1,(ns-1))),"sparseMatrix")
  M <- Diagonal(n) %x% C
  ## Collect nodes corresponding to identical positions
  A <- modMat(fac) ## Alternatively: use as(fac,"sparseMatrix")
  P <- t(A)%*%M%*%A
  P@x[] <- 1
  as(P,"dgTMatrix")
}
polygonGrid <- function(object,...){
  ## -------------- Construct corresponding grid and pattern-matrix
  ##triangles <- triangles$triangles
  ##grid <- unique(triangles)
  ns <- nsidesPolygons(object)
  n <- nrow(object)/ns
  ch <- do.call("paste",object)
  fac <- factor(ch,levels=unique(ch))
  P <- .incidenceMatrix(fac,ns)
  grid <- object[!duplicated(fac),]
  attr(grid,"pattern") <- P
  attr(grid,"polygons") <- fac
  attr(grid,"nsides") <- ns
  class(grid) <- c("polygonGrid","data.frame")
  grid
}
triangularGrid <- function(data,triangles=makeTriangles(data,...),...){
  polygonGrid(triangles)
}
squareGridSubdiv <- function(data,squares=makeSquares(data,...),...){
  polygonGrid(squares)
}
rectGrid <- function(lon,lat){
  nx <- length(lon)
  ny <- length(lat)
  e <- expand.grid(lon = lon, lat = lat)
  m <- matrix(1:nrow(e), nx, ny)
  i <- t(matrix(c(m[-nx, -ny], m[-1, -ny], m[-1, -1], m[-nx, -1]), 
                ncol = 4))
  grid <- e
  attr(grid, "pattern") <- .incidenceMatrix(factor(i), 4)
  attr(grid, "polygons") <- factor(i)
  attr(grid, "nsides") <- 4
  class(grid) <- c("polygonGrid", "data.frame")
  grid
}
squareGrid <- function(data,n=32,square0=initialSquare(data,...),
                       center=data.frame(lon=mean(range(data$lon)), lat=mean(range(data$lat))),
                       rot=0,km,...){
  center <- as.list(center)
  if(!missing(km)){
    deg <- km/6378.1*180/pi
    lon <- seq(min(square0$lon),max(square0$lon),by=deg)
    lat <- seq(min(square0$lat),max(square0$lat),by=deg)    
  } else {
    lon <- seq(min(square0$lon),max(square0$lon),length=n)
    lat <- seq(min(square0$lat),max(square0$lat),length=n)
  }
  grid <- rectGrid(lon,lat)
  grid <- sphRot(grid,rot,1)
  grid <- sphRot(grid,-center$lat,2)
  grid$lon <- grid$lon+center$lon
  grid
}

## ==============================================================
## Methods for class "polygonGrid"
## ==============================================================
as.data.frame.polygonGrid <- function(x,...){
  attr(x,"pattern") <- NULL
  attr(x,"polygons") <- NULL
  class(x) <- "data.frame"
  x
}
plot.polygonGrid <- function(x,type="l",add=FALSE,...){
  P <- attr(x,"pattern")
  i <- P@i+1; j <- P@j+1 ## edges
  keep <- i<j
  i <- i[keep]; j <- j[keep]
  class(x) <- "data.frame"
  obj <- rbind(x,NA)[  t(cbind(i,j,length(i)+1)) , ]
  if(add)plot <- points
  plot(obj,type=type,...)
}
"[.polygonGrid" <- function(x, i,...){
  ##if(is.character(i))return(NextMethod("["))
  P <- attr(x,"pattern")
  fac <- attr(x,"polygons")
  ns <- attr(x,"nsides")

  x <- as.data.frame(x)
  x <- x[i,...]

  subgrid <- !any(duplicated(i)) & is.numeric(i) 
  subgrid <- subgrid | is.logical(i) 
  
  if(subgrid){
    attr(x,"pattern") <- P[i,i]
    newlevels <- levels(fac)[i]
    k <- rep(as.logical(ns == colSums(matrix(fac %in% newlevels,ns))),each=ns)
    g <- factor(fac[k],levels=newlevels)
    attr(x,"polygons") <- g
    class(x) <- c("polygonGrid","data.frame")
  } else {
    warning("Loosing neighborhood and polygon information")
  }
  x
}

as.polygons <- function(x,...)UseMethod("as.polygons")
as.polygons.polygonGrid <- function(x,...){
  fac <- attr(x,"polygons")
  x <- as.data.frame(x)
  ans <- x[fac,]
  ns <- attr(x,"nsides")
  cl <- switch(as.character(ns),
               "4"="squares",
               "3"="triangles",
               stop("Unknown polygon"))
  class(ans) <- c(cl,"polygons",class(ans))
  ans
}
## Points in polygon
cut.polygonGrid <- function(x,polygons=splitMap(makeMap(x,"worldHires")),...){
  if(is.list(polygons)){
    if(!is.list(polygons[[1]])) polygons <- list(polygons)
  } else stop("polygons must be list or data.frame")
  if(is.null(names(polygons))) names(polygons) <- 1:length(polygons)
  require(sp)
  mat <- sapply(polygons,function(y)point.in.polygon(x$lon,x$lat,y[[1]],y[[2]]))
  fac <- mat %*% (1:ncol(mat))
  fac[!fac] <- NA
  fac <- factor(names(polygons)[fac])
  fac
}
## Model matrix.
## Assign the average values at the corners to each polygon.
model.matrix.polygonGrid <- function(object,...){
  nsides <- attr(object,"nsides")
  fac <- attr(object,"polygons")
  n <- length(fac)/nsides
  i <- rep(1:n,each=nsides)
  j <- as.numeric(fac)
  A <- spMatrix(n,nrow(object),i=i,j=j,x=0*j+1/nsides)
  A
}
image.polygonGrid <- function(x,y,...){
  A <- model.matrix(x)
  polygons <- as.polygons(x)
  image(polygons,as.vector(A%*%y),...)
}

boundary <- function(x,...)UseMethod("boundary")
boundary.polygonGrid <- function(x){
  P <- attr(x,"pattern")
  nb <- as.numeric(rowSums(P))
  which(nb<max(nb))
}

## ==============================================================
## Map stuff
## ==============================================================
makeMap <- function(x,map="worldHires",...){
  require(maps);require(mapdata)
  xlim <- range(x$lon)
  ylim <- range(x$lat)
  map.poly  <-  map(map,xlim=xlim,ylim=ylim,
                    fill=TRUE,col="grey",
                    interior=FALSE,plot=FALSE)
  map.poly
}

splitMap <- function(map){
  i <- cumsum(is.na(map$x))
  ans <- split(as.data.frame(map[c("x","y")]),i)
  names(ans) <- map$names
  lapply(ans,na.omit)
}

## Find connected components of a grid ordered after size of components.
connectedComponents <- function (gr0) 
{
  P <- attr(gr0, "pattern")
  P <- as(P,"dgCMatrix")
  col <- function(i)P@i[ (P@p[i]+1):(P@p[i+1]) ] +1
  root <- Vectorize(function(i){while(p[i]!=i)i <- p[i];i})
  p <- 1:ncol(P)
  for(i in seq(p)){
    ## i and all of its neighbors
    n <- col(i)
    ## Merge trees by letting each root point to the lowest root
    r <- root(n)
    p[r] <- min(r)
  }
  ## Label points with tree roots
  for(i in seq(p))p[i] <- root(p[i])
  split(1:nrow(P),p)
}

## Workaround: To make imageanimations work with raw polygons (no polygonGrid).
model.matrix.polygons <- function(x,...){
  .symDiagonal(summary(x)$n)
}
as.polygons.polygons <- function(x,...)x
