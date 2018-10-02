\name{unionP}
\alias{unionP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
unionP
}
\description{
Probability that at least one individual came from the location (union of probabilities)
}
\usage{
unionP(pdR)
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
asn = pdRaster(r$isoscape.rescale$mean,r$isoscape.rescale$sd,unknown=un,mask=naMap)
# Probability that at least one individual came from the location (union of probabilities)
unionP(asn)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line