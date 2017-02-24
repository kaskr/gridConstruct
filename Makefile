PACKAGE=gridConstruct
VERSION=1.0
TARBALL=${PACKAGE}_${VERSION}.tar.gz
ZIPFILE=${PACKAGE}_${VERSION}.zip

all:
	make doc-update
	make build-package
	make install
	make pdf

doc-update:
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\",roclets = c(\"collate\", \"rd\"))" | R --slave

build-package:
	R CMD build --resave-data=no $(PACKAGE)

install:
	R CMD INSTALL --preclean $(PACKAGE)

unexport TEXINPUTS
pdf:
	rm -f $(PACKAGE).pdf
	R CMD Rd2pdf --no-preview $(PACKAGE)

shp-update:
	R --slave < shp-update.R

## Get map data from OpenStreetMapData
## ODbL 1.0 license:
## You are free:
##      To Share (copy, distribute and use the database)
##      Create (To produce works from the database)
##      Adapt (To modify, transform and build upon the database)
## As long as you:
##      Attribute, Share-Alike, Keep open.
land-polygons-complete-4326.zip:
	wget http://data.openstreetmapdata.com/land-polygons-complete-4326.zip

land-polygons-complete-4326: land-polygons-complete-4326.zip
	unzip land-polygons-complete-4326.zip

TM_WORLD_BORDERS-0.3.zip:
	wget http://www.thematicmapping.org/downloads/TM_WORLD_BORDERS-0.3.zip

TM_WORLD_BORDERS-0.3: TM_WORLD_BORDERS-0.3.zip
	mkdir TM_WORLD_BORDERS-0.3
	cd TM_WORLD_BORDERS-0.3; unzip ../TM_WORLD_BORDERS-0.3.zip

## GPL !!! :)
gshhg-shp-2.3.6.zip:
	wget https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.6.zip

gshhg-shp-2.3.6: gshhg-shp-2.3.6.zip
	mkdir gshhg-shp-2.3.6
	cd gshhg-shp-2.3.6; unzip ../gshhg-shp-2.3.6.zip
