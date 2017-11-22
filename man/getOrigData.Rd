\name{getOrigData}
\alias{getOrigData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getOrigData function
}
\description{
This function is to load isotope data from known origin in isOrigin package and display its summary in console.
}
\usage{
getOrigData()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
data frame. Known-origin data from isOrigin package
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
\code{\link[isOrigin]{knownOrig}}
\code{\link[isOrigin]{subOrigData}}
}
\examples{
d <- getOrigData() # it will show a summary of this data in console and assign data to d
data(wrld_simpl)
plot(wrld_simpl,axes=T)
points(d$Longitude,d$Latitude,col="red",cex=0.5)
}

