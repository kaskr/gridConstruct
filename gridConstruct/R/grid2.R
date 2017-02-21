## Add transformed coordinates to DATRASraw object:
addCoords <- function(d, map=Europe) {
    ## Define coordinates of HH records:
    tmp <- d[["HH"]][c("lon", "lat")]
    names(tmp) <- c("x", "y")
    coordinates(tmp) <- ~x + y
    proj4string(tmp) <- CRS("+proj=longlat")
    ## Transform coordinates
    proj4 <- proj4string(Europe)
    tmp <- spTransform(tmp, CRS(proj4))
    ## Attach coordinates to DATRASraw
    d[["HH"]]$coords <- tmp
    d
}

## Get grid from DATRASraw
## FIXME: Improve interface !

makeGrid <- function(d,
                     km=1,
                     scale=1,
                     center=colMeans( apply(as.data.frame(d[["HH"]]$coords), 2, range) ),
                     map=Europe) {
    if(is.null(d[["HH"]]$coords)) d <- addCoords(d)
    coords <- d[["HH"]]$coords
    dx <- ceiling( scale * diff(range(d[["HH"]]$coords$x)) / 2 )
    dy <- ceiling( scale * diff(range(d[["HH"]]$coords$y)) / 2 )
    gr <- expand.grid(x = seq(-dx, dx, by=km),
                      y = seq(-dy, dy, by=km) )
    gr <- as.data.frame( t( (t(gr) + center) ) )
    coordinates(gr) <- ~x + y
    proj4string(gr) <- CRS(proj4string(map))
    wet <- over(gr, map)
    gr <- gr[is.na(wet$FID_Europe), ]
    gr <- as.data.frame(gr)
    attr(gr, "map") <- map
    class(gr) <- c("makeGrid", class(gr))
    gr
}
plot.makeGrid <- function(x,map=TRUE,...){
    plot(x$x, x$y)
    if(map)plot(attr(gr,"map"), add=TRUE, col="grey")
}


if(FALSE){
    ## Experimental setup of new grid construction methods.
    ## Read data
    library(DATRAS)
    d <- readICES("~/Projects/muslinger/data/Exchange_Data_all_years_v1.csv")
    d <- subset(d,lon<9.5 & lat>56)
    library(gridConstruct)

    gr <- makeGrid(d, km=1)
    plot(gr)
    
    plot(qw$lon,qw$lat)
    plot(map, col="grey", add=TRUE)
    plot(map, col="grey", add=TRUE, xlim=range(qw$lon), ylim=range(qw$lat))

    plot(map, col="grey", xlim=range(qw$lon), ylim=range(qw$lat))
    points(qw$lon,qw$lat)
#####
    ## devtools:::install_github("s-u/fastshp")
    ## shape2 <- read.shp(system.file("shpfiles/Europe",
    ##                                package="gridConstruct"))

    
}
