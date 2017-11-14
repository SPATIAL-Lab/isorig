\name{subOrigData}
\alias{subOrigData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
subOrigData function
}
\description{
This function is to subset the isotope data from known origin in isorig package
}
\usage{
subOrigData(taxon = NULL, group = NULL, mask = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxon}{
one string or a string vector. Choosing from taxon names of known origin data. Details can be seen in \code{\link[isorig]{knownOrig}}
}
  \item{group}{
one string or a string vector. Choosing from group names of known origin data. Details can be seen in \code{\link[isorig]{knownOrig}}
}
  \item{mask}{
SpatialPolygonsDataFrame. Constrain the investigated area. If this is not provided, default of whole world is used.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
data selected from known origin data from this package
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

\seealso{
\code{\link[isorig]{knownOrig}}
\code{\link[isorig]{origData}}
}

\examples{
## WITHOUT mask
# extract data of taxon of Jackdaw, Partridge and Willow Grouse
d1 = subOrigData(taxon = c("Jackdaw", "Partridge", "Willow Grouse"))
summary(d1)

# extract data of group of [1]insect and [2]bird
d2 = subOrigData(group = c("insect","bird"))
summary(d2)

# extract data of group of human which only has hair
d3 = subOrigData(group = "human")
summary(d3)

# extract data of taxon of hair from human
d4 = subOrigData(taxon = "hair")
summary(d4)

## Because in knownOrig data base in isorig package, human only has hair samples, thus d4 and d3 above are the same
summary(d4-d3)

## WITH mask
data(naMap) # load North America SpatialPolygonsDataFrame data in isorig package
d5 = subOrigData(taxon = "Mallard", mask = naMap) # No isotope data found in mask you choose!
d6 = subOrigData(taxon = c("Monarch","AMRE","AMRO"), mask = naMap)
d7 = subOrigData(group = "human", mask = naMap)
}
