\name{pdRaster}
\alias{pdRaster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
assign function
}
\description{
Assign unknown tissue isotope to calibrated precipitation isoscape
}
\usage{
pdRaster(rescaled_raster, unknown, mask = NULL, genplot=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rescaled_raster}{
raster. The tissue-calibrated isoscape raster (mean and sd) obtained from rescale function, or user provided
}
  \item{unknown}{
data.frame. First column should be unique sample ID, and second column should be tissue d2H values of the individuals for which the assignments will be made.
}
\item{mask}{
SpatialPolygonsDataFrame. This mask will constrain the assignment area. If this is not provided, a default of mask rescaled_raster is used.
}
\item{genplot}{
logical. Generate summary plots in R? (TRUE [default] or FALSE)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
one-layer or stacked-layer raster with assignment map
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[isOrigin]{rescale}}
}
\examples{
data(naMap)
data(d2h_world)
d <- subOrigData(taxon = c("Homo sapiens"), mask = naMap)[,2:4]
r <- rescale(orig = d, precip = d2h_world)
id <- "smile"
d2H <- -80
un <- data.frame(id,d2H)
asn <- pdRaster(r$isoscape.rescale,un,mask=naMap)
}
