\name{assignment}
\alias{assignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
assignment function
}
\description{
Assign unknown tissue isotope to calibrated precipitation isoscape
}
\usage{
assignment(rescaled_raster, unknown, d2Hcol, IDcol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rescaled_raster}{
the tissue-clibrated isoscape raster (mean and sd) obtained from rescale function
}
  \item{unknown}{
unknown tissue isotope table: Could be a .csv filename (and directory, if applicable). First column should be ID, second column should be tissue d2H values of the individuals for which the assignments will be made. Or be a data.frame table with information mentioned above.
}
\item{mask}{
a SpatialPolygonsDataFrame that is used to constrain the investigated area. If this is not provided, default of whole world is used.
}
\item{genplot}{
Generate summary plots? (T or F)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(naMap)
data(d2h_world)
d = subOrigData(taxon = c("hair"), mask = naMap)
r = rescale(orig = d, precip = d2h_world)
id = "smile"
d2H = -80
un = data.frame(id,d2H)
asn = assignment(r$precip.rescale,un,mask=naMap)
}
