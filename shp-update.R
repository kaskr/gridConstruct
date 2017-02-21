library(rgdal)
Europe <- readOGR("shpfiles/Europe", "Kystlinie")
save(Europe, file="gridConstruct/inst/Europe.RData")
