library(rgdal)
Europe <- readOGR("shpfiles/Europe", "Kystlinie")
new_proj <- sub("units=m", "units=km", proj4string(Europe))
Europe <- spTransform(Europe, CRS(new_proj))
save(Europe, file="gridConstruct/inst/Europe.RData")

map <- readOGR("land-polygons-complete-4326", "land_polygons")

map <- readOGR("TM_WORLD_BORDERS-0.3", "TM_WORLD_BORDERS-0.3")

map <- readOGR("gshhg-shp-2.3.6/GSHHS_shp/f", "GSHHS_f_L1")

map0 <- subset(map, 0 < coordinates(map)[,1] & coordinates(map)[,1] < 14)

library(raster) 
data("wrld_simpl", package="maptools")

## Crop to the desired extent, then plot
out <- crop(wrld_simpl, extent(0, 14, 50, 60))
plot(out, col="khaki", bg="azure2")


out <- crop(map, extent(0, 14, 50, 60))
plot(out, col="khaki", bg="azure2")
