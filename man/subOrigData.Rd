\name{subOrigData}
\alias{subOrigData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
subOrigData function
}
\description{
This function will permit subset the known-origin stable isotope dataset included in this package
}
\usage{
subOrigData(taxon = NULL, group = NULL, mask = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxon}{
one string or a string vector. Select species names of known origin data. Details can be seen in \code{\link[isOrigin]{knownOrig}}
}
  \item{group}{
one string or a string vector. Select group names of known origin data. Details can be seen in \code{\link[isOrigin]{knownOrig}}
}
  \item{mask}{
SpatialPolygonsDataFrame. This mask will constrain the investigated area. If this is not provided, default of whole world is used.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Data selected as data frame from known-origin data from this package will be returned. It has three columns: 1. Longitude, 2. Latitude, 3. d2H, which is ready to be fed into calRaster function.
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

\examples{
## WITHOUT mask
# extract data of taxon of Jackdaw, Partridge and Willow Grouse
d1 = subOrigData(taxon = c("Danaus plexippus", "Setophaga ruticilla", "Turdus migratorius"))
summary(d1)

# extract data of group of [1]insect and [2]bird
d2 = subOrigData(group = c("insect","bird"))
summary(d2)

# extract data of group of human which only has hair
d3 = subOrigData(group = "human")
summary(d3)

# extract data of taxon of hair from human
d4 = subOrigData(taxon = "Homo sapiens")
summary(d4)

## Because in knownOrig data base in isOrigin package, human only has hair samples, thus d4 and d3 above are the same
summary(d4-d3)

## WITH mask
data(naMap) # load North America SpatialPolygonsDataFrame data in isOrigin package
d5 = subOrigData(taxon = "Turdus merula", mask = naMap) # No isotope data found in mask you choose!
d6 = subOrigData(taxon = c("Danaus plexippus", "Setophaga ruticilla", "Turdus migratorius"), mask = naMap)
d7 = subOrigData(group = "human", mask = naMap)
}
