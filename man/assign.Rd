\name{assign}
\alias{assign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
assign function
}
\description{
Assign unknown tissue isotope to calibrated precipitation isoscape
}
\usage{
assign(rescaled_raster, unknown, mask = NULL, genplot=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rescaled_raster}{
raster. the tissue-clibrated isoscape raster (mean and sd) obtained from rescale function, or user provided
}
  \item{unknown}{
data.frame. First column should be ID, second column should be tissue d2H values of the individuals for which the assignments will be made.
}
\item{mask}{
SpatialPolygonsDataFrame. constrain the investigated area. If this is not provided, default of mask rescaled_raster is used.
}
\item{genplot}{
logical. Generate summary plots in R? (T [default] or F)
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
\code{\link[isorig]{rescale}}
}
\examples{
data(naMap)
data(d2h_world)
d = subOrigData(taxon = c("hair"), mask = naMap)[,2:4]
r = rescale(orig = d, precip = d2h_world)
id = "smile"
d2H = -80
un = data.frame(id,d2H)
asn = assign(r$precip.rescale,un,mask=naMap)
}
