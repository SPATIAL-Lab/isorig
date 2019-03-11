\name{QA_2isotopes}
\alias{QA_2isotopes}
\title{
Test the efficacy of geographic assignment using two certain type of isoscape
}
\description{
What is the power of a certain isoscape used for geographic assignment? Using the known origin data and isoscapes from two kind of isotopes as input to test it. You will get the population accuracy, precision and probability density. Joint probability density of geographic assignment from two isotopes is caculated by the function of \code{\link[isOrigin]{jointP}}.
}
\usage{
QA_2isotopes(isoscape1, isoscape2, known1, known2, valiStation, valiTime, setSeed = T)
}
\arguments{
  \item{isoscape1}{raster. Environmental isoscape. Two layers: the first one is mean and the second one is standard deviation.}
  \item{isoscape2}{raster. Environmental isoscape. Two layers: the first one is mean and the second one is standard deviation.}
  \item{known1}{
data.frame. Known origin tissue isotope values which sould be the same isotope as that in isoscape1. Three columns: 1. longitude, 2. latitude, 3. isotope values.
}
  \item{known2}{
data.frame. Known origin tissue isotope values which sould be the same isotope as that in isoscape2. Three columns: 1. longitude, 2. latitude, 3. isotope values. This should have the same locations as that in known1
}
  \item{valiStation}{
  numeric. How many stations of the known origin with tissue isotope are used for validation.
}
  \item{valiTime}{
numeric. How many times do you want to randomly draw validation stations and run the validation.
}
  \item{setSeed}{
Do you want to set.seed() when you run the randomly draw validation stations? If yes and your input data are the same, the output would be exactely the same.
}
}

\value{
\item{val_stations}{numeric. An X*Y data.frame of validation station IDs for all valiTime. X is the valiTime and Y is the validation station IDs.}
\item{pd_bird_val}{
numeric. An X*Y data.frame containing the posterior probability density for the validation stations. X is the simulation numbers = valiTime and Y is the total number of validation stations (valiStation).
}
  \item{prption_byArea}{
  numeric. An X*Y data.frame shows the population-level accuracy that is measured as the proportion of validation individuals in which the known origin is contained within the top probability ranging from 0.01 to 0.99 with increment of 0.01 (99 probabilities which is Y). X is the valiTime. Higher proportion with lower top probability means higher accuracy.
}
  \item{prption_byProb}{
numeric. An X*Y data.frame shows the population-level accuracy that is measured as the proportion of validation individuals in which the known origin is contained within the top percent of area ranging from 0.01 to 0.99 with increment of 0.01 (99 percentage which is Y). X is the valiTime.
}
  \item{precision}{
list. The length of the list is valiTime which reprsents the population-level precision for each validation. It is assessed by the areal proportion of the total surface area covered by the assignment region for each top percent of probability density (threshold). Thresholds from 0.01 to 0.99 with increment of 0.01 are tested. Each elements of the list is a numerical data.frame with dimention of X*Y. X is the 99 threshold. The length of Y is valiStation. Lower proportion with lower threshold in this assessment suggests higher precision. 
}
  \item{random_prob_density}{
   Random probability density based on the size of the environmental isoscape (i.e. 1 divided by the total number of grid cells of the isoscape). 
}
}

\note{
Please see Ma et al., 2019 for details of these values returned and the methodology. 
}


\references{
Vander Zanden, H.B., Wunder, M.B., Hobson, K.A., Van Wilgenburg, S.L., Wassenaar, L.I., Welker, J.M. and Bowen, G.J., 2014. Contrasting assignment of migratory organisms to geographic origins using long-term versus year-specific precipitation isotope maps. Methods in Ecology and Evolution, 5(9), pp.891-900.
}

\seealso{
\code{\link[isOrigin]{QA}}
\code{\link[isOrigin]{plot.QA}}
}

\examples{
# load data
data(naMap) # North America 
data(d2h_world) # precipitation hydrogen isotope of the world
data(d18o_world) # precipitation oxygen isotope of the world
data(bird_isotope) # oxygen and hydrogen isotopes of known-origin bird

# crop the world hydrogen data to North America
r1 <- crop(d2h_world, naMap)
r2 <- crop(d18o_world, naMap)
plot(r1)
plot(r2)

# convert 2 standard deviation from d2h_world and d18o_world to 1 standard deviation
r1[[2]] <- r1[[2]]/2
r2[[2]] <- r2[[2]]/2

# seperate the hydrogen and oxygen isotope for the known-origin bird
bird_d2h <- bird_isotope[1:20,c("Longitude", "Latitude", "d2H")]
bird_d18o <- bird_isotope[1:20,c("Longitude", "Latitude", "d18O")]

# run quality assessment based hydrogen isotope from precipitation and known-origin bird
d2h_QA <- QA_2isotopes(isoscape1 = r1, isoscape2 = r2, known1 = bird_d2h, kown2 = bird_d18o, valiStation = 2, valiTime = 5, setSeed = T)
                    
# plot the QA result
plot(d2h_QA)
}
