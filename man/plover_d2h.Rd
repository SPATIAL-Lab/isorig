\name{plover_d2h}
\alias{plover_d2h}
\docType{data}
\title{
plover_d2h}
\description{
Hydrogen isotope of the feather of known origin mountain plover and associated precipitation hydrogen isotope of growing season}
\usage{data("plover_d2h")}
\format{
  A data frame with 465 observations on the following 6 variables.
  \describe{
    \item{\code{age}}{a factor with levels \code{adult} \code{adultR} \code{chick} \code{chickR}}
    \item{\code{year}}{a numeric vector}
    \item{\code{lat}}{a numeric vector}
    \item{\code{long}}{a numeric vector}
    \item{\code{d2h_feather}}{a numeric vector}
    \item{\code{d2h_precipitation}}{a numeric vector}
  }
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
\examples{
data(plover_d2h)
plot(plover_d2h)
}
\keyword{datasets}
