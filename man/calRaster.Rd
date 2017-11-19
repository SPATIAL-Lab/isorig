\name{calRaster}
\alias{calRaster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to rescale isoscape
}
\description{
This function relies on known-origin tissue data to calibrate between d2H values in precipitation and those in the tissue using simple regression to rescale the precipitation raster into tissue values and generate mean and standard deviation rasters. This will output plots as pdf including regression model and rescaled isoscape mean and standard deviation in output folder under working directory.
}
\usage{
calRaster(known, isoscape, mask = NULL, interpMethod = 2, NA.value = NA, ignore.NA = T)
}
\arguments{
  \item{known}{
data.frame. known-origin data that should contain the following 3 columns: 1st is longitude, 2nd is latitude, 3rd is tissue isotope value. Optional known-origin data can be queried using the \code{\link[isOrigin]{knownOrig}} function.
}
  \item{isoscape}{
 raster. isoscape raster
}
  \item{mask}{
 SpatialPolygonsDataFrame. This mask will constrain the isoscape area. If this is not provided, default of mask of isoscape in this argument is used.
}
  \item{interpMethod}{
numeric. 1 or 2. Designate oen of two methods for extracting
values from isoscape raster based on known-origin position. If 1, values for the cell
in which a point falls are returned. If 2 (default) the returned values are
interpolated (bilinear interpolation) from the values of the four nearest raster cells.
}
\item{NA.value}{
 NA [default] or numeric. Identify the value that was used to represent an absence of data in the isoscape raster map. Commonly used values include NA, -9999, 9999, etc. If you are not sure what value was used in your isoscape, you can load the isoscape raster in R and run \code{\link[raster]{click}}(raster) to click the place without value.
}
\item{ignore.NA}{
 logical. (TRUE or FALSE [default]). If NA values are extracted from isoscape at the location of known-origin, do you want to ignore these values and proceed with calRaster function?
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{isoscape.rescale }{Rescaled isoscape based on known origin data and isosape, two layer raster, one is mean the other is standard deviation which is a univeral value for the map from error of linear regression model.}
\item{lm.data }{known origin data and correlated isoscape value used for linear regression model}
\item{lm.model}{linear regression model}

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
\code{\link[isOrigin]{assign}}
}
\examples{
data("naMap")
data("d2h_world")
# the columns of d should be longitude, latitude, and tissue isotope
d <- subOrigData(taxon = c("hair"), mask = naMap)[,2:4]
r <- calRaster(known = d, isoscape = d2h_world, mask = naMap)
}