\name{getOrigData}
\alias{getOrigData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getOrigData function
}
\description{
This function is to load isotope data from known origin in isorig package and display its summary in console.
}
\usage{
getOrigData()
}
%- maybe also 'usage' for other objects documented here.
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
\code{\link[isorig]{knownOrig}}
\code{\link[isorig]{subOrigData}}
}
\examples{
d = getOrigData() # it will show a summary of this data in console and assign data to d
data(wrld_simpl)
plot(wrld_simpl,axes=T)
points(d$Longitude,d$Latitude,col="red",cex=0.5)
}

