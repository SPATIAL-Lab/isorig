\name{pdRaster}
\alias{pdRaster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
assign function
}
\description{
Assign unknown tissue isotope to calibrated precipitation isoscape.
}
\usage{
pdRaster(r, unknown, mask = NULL, genplot=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
RasterStack or RasterBrick with two layers. The tissue-calibrated isoscape raster (mean and sd) obtained from rescale function, or user provides
}
  \item{unknown}{
data.frame. First column should be unique sample ID, and second column should be tissue d2H values of the individuals for which the assignments will be made.
}
\item{mask}{
SpatialPolygonsDataFrame. This mask will constrain the assignment area. If this is not provided, a default of mask of rescaled.mean is used.
}
\item{genplot}{
logical. Generate summary plots in R? (TRUE [default] or FALSE)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns an object of class RasterStack or RasterLayer that includes a probability density raster for each individual in "unknown". pdRaster writes individual rasters in GTiff format to disk in the working directory for each individual in "unknown".
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
\code{\link[isOrigin]{calRaster}}
}
\examples{
data(naMap)
data(d2h_world)
d <- subOrigData(taxon = c("Homo sapiens"), mask = naMap)
r <- calRaster(known = d, isoscape = d2h_world, mask = naMap)
id <- "smile"
d2H <- -80
un <- data.frame(id,d2H)
asn <- pdRaster(r$isoscape.rescale, un,mask=naMap)
}
