\name{knownOrig}
\alias{knownOrig}
\docType{data}
\title{
Isotope values for tissues of known origin
}
\description{
This dataset consists of hydrogen values including human hair, insect chitin and birds feathers, with location information.
}
\usage{data("knownOrig")}
\format{
  A SpatialPointsDataFrame
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\seealso{
\code{\link[isOrigin]{subOrigData}}
}
\examples{
library(maptools)
data("knownOrig")
data("wrld_simpl")
plot("wrld_simpl")
points(knownOrig$Longitude,knownOrig$Latitude,col="red",cex=0.5)
}
\keyword{datasets}
