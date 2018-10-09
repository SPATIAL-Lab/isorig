\name{jointP}
\alias{jointP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
jointP
}
\description{
Joint probability for individuals of common origin (product of probabilities)
}
\usage{
jointP(pdR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdR}{
probability density map which must be raster or RasterLayer
}
}

\examples{
d <- getOrigData()
# load North America SpatialPolygonsDataFrame data in isOrigin package
data("naMap")
data(d2h_world)
d = subOrigData(taxon = c("Homo sapiens"), mask = naMap)
# exclude some outliers. This step is optional, which depends on your data quality
dd = d[d$Longitude<(-80),]
r = calRaster(known = dd, isoscape = d2h_world, mask = naMap)
# four unknown-origin example
id = c("A", "B", "C", "D")
d2H = c(-110, -180, -130, -150)
un = data.frame(id,d2H)
# assignment
asn = pdRaster(r$isoscape.rescale,unknown=un,mask=naMap)
# joint probability for individuals of common origin
jointP(asn)
}
}
