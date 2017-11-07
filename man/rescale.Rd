\name{rescale}
\alias{rescale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rescale
}
\description{
Function to calibrate between precipitation water isotope and tissue isotope
}
\usage{
rescale(orig, precip)
}
\arguments{
  \item{orig}{
known origin.Two options: (1) If user provid known origin data, this should be a filename (with directory, if applicable) of known origin. It should be as a csv file with 3 columns: 1st is longitude, 2nd is latitude, 3rd is tissue isotope value. The first row should be header of "Longitude", "Latitude", "d2H". If you have a shape file, you should select these 3 attributions and save it as a csv file with format mentioned above. (2) known origin data provided by this pacakge. This should be output from function: subOrigData, with format of data.frame
}
  \item{precip}{
 raster of isoscape
}
  \item{mask}{
 a SpatialPolygonsDataFrame that is used to constrain the investigated area. If this is not provided, default of whole world is used.
}
  \item{interpMethod}{
interpMethod = numeric. 1 or 2. This is method for extracting
values from precipitation raster based on know origin position. If 1 values
for the cell a point fall s in are returned. If 2 (default) the returned values are
interpolated from the values of the four nearest raster cells.
}
\item{NA.value}{
 (NA [default] or numeric) what is the value for no data in isoscape raster map, could be NA, -9999, 9999, etc. If you are not sure about it, you can read the raster in R and run \code{\link[raster]{click}}(raster) to click the place without value.
}
\item{ignore.NA}{
 (T or F [default]) If NA value are extracted at location of known origin, do you want to ignore these value and proceed with rescale function?
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
r = rescale(orig = d, precip = d2h_world, mask = naMap)
}
