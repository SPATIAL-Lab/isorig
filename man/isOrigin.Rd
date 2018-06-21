\name{isOrigin}
\alias{isOrigin}

\title{isOrigin}

\examples{
library(isOrigin)
library(ggplot2)
d <- getOrigData()
library(sp)
data("naMap") # load North America SpatialPolygonsDataFrame data in isOrigin package
d5 = subOrigData(taxon = "Turdus merula", mask = naMap) # No isotope data found in mask you choose!
d6 = subOrigData(taxon = c("Danaus plexippus", "Setophaga ruticilla", "Turdus migratorius"), mask = naMap)
d7 = subOrigData(group = "human", mask = naMap)
library(raster)
data(d2h_world)
d = subOrigData(taxon = c("Homo sapiens"), mask = naMap)
# exclude some outliers. This step is optional, which depends on your data quality
dd = d[d$Longitude<(-80),]
r = calRaster(known = dd, isoscape = d2h_world, mask = naMap)
# one unknown-origin example
id = c("Harry_Potter", "Superman")
d2H = c(-110, -180)
un = data.frame(id,d2H)
asn = pdRaster(r$isoscape.rescale$mean,r$isoscape.rescale$sd,unknown=un,mask=naMap)
}
