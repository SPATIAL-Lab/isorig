\name{d2h_world}
\alias{d2h_world}
\docType{data}
\title{
data: d2h_world
}
\description{
This file is the global growing season d2H precipitation isoscape available on waterisotopes.org and reported in Bowen et al., 2005.
}
\usage{data("d2h_world")}
\format{
  The format is raster.
[class]       RasterLayer;
[dimensions]   419, 1080, 452520  (nrow, ncol, ncell);
[resolution]   0.333333, 0.333333  (x, y);
[extent]       -180, 179.9996, -56, 83.66652  (xmin, xmax, ymin, ymax);
[coord. ref.]  +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0;
[data source]  in memory;
[names]        gs_d2h;
[values]       -169.26, 34.8  (min, max);
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
data(d2h_world)
plot(d2h_world)
}

