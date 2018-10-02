\name{qtlRaster}
\alias{qtlRaster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binary reclassification
}
\description{
Select top X percent of probability surface (defined by area percentage or cumulative probability)
}
\usage{
qtlRaster(pdR, threshold, thresholdType)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdR}{
probability density map which must be raster or RasterLayer
}
  \item{threshold}{
numeric from 0 to 1, top percentage to choose
}
  \item{thresholdType}{
numeric 1 or 2. 1 is to choose top percentage by area; 2 is to choose top percentage by cumulative probability
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
asn = pdRaster(r$isoscape.rescale$mean,r$isoscape.rescale$sd,unknown=un,mask=naMap)
# top 10% of probability surface (defined by % area)
qtlRaster(asn, threshold = 0.1, thresholdType = 2)
# top 10% of probability surface (defined by % cumulative probability)
qtlRaster(asn, threshold = 0.1, thresholdType = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line