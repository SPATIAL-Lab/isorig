\name{knownOrig}
\alias{knownOrig}
\docType{data}
\title{
Isotope values for tissues of known origin
}
\description{
This dataset consists of hydrogen values for 2178 samples including human hair, insect chitin and birds feathers, with location information.
}
\usage{data("knownOrig")}
\format{
  A data frame with 2151 observations on the following 9 variables.
  \describe{
    \item{\code{ID}}{a numeric vector. 1 to 2151}
    \item{\code{Longitude}}{a numeric vector}
    \item{\code{Latitude}}{a numeric vector}
    \item{\code{d2H}}{a numeric vector}
    \item{\code{Taxon}}{a numeric vector with 80 levels (used for subOrigData function). \code{"Danaus plexippus"            "Setophaga ruticilla"         "Turdus migratorius"
 "Setophaga coronata auduboni" "Poecile atricapillus"        "Thryomanes bewickii"
 "Thryothorus ludovicianus"    "Spizella passerina"          "Geothlypis trichas"
"Setophaga pensylvanica"      "Baeolophus bicolor"          "Vermivora chrysoptera"
"Catharus guttatus"           "Setophaga citrina"           "Geothlypis formosa"
"Geothlypis tolmiei"          "Oreothlypis ruficapilla"     "Cardinalis cardinalis"
"Oreothlypis celata"          "Junco hyemalis oregonus"     "Seiurus aurocapilla"
"Vireo olivaceus"             "Melospiza melodia"           "Catharus ustulatus"
"Catharus fuscescens"         "Vireo griseus"               "Cardellina pusilla"
"Hylocichla mustelina"        "Icteria virens"              "Setophaga petechia"
"Melozone aberti"             "Vermivora cyanoptera"        "Passer domesticus"
"Aimophila ruficeps"          "Poecile carolinensis"        "Troglodytes aedon"
"Dumetella carolinensis"      "Mniotilta varia"             "Lanius ludovicianus"
"Anthus spragueii"            "Euphagus carolinus"          "Empidonax minimus"
"Aythya affinis"              "Oreothlypis peregrina"       "Cyanistes caeruleus"
"Phasianus colchicus"         "Lagopus lagopus"             "Tetrao tetrix"
"Dryocopus maritus"           "Serin serin"                 "Vanellus vanellus"
"Corvus corone"               "Turdus merula"               "Corvus corone "
"Corvus monedula"             "Columba palumbus "           "Turtle Dove"
"Tetrastes bonasia"           "Perdix perdix"               "Anas platyrhyncos"
"Branta canadensis"           "Columba livia"               "Numenius arguata"
"Turdus pilaris"              "Turdus iliacus"              "Turdus philomelos"
"Fringilla coelebs"           "Buteo lagopus"               "Accipiter striatus"
"Falco sparverius"            "Accipiter gentillis"         "Accipiter cooperii"
"Buteo jamaicensis"           "Buteo platypterus"           "Buteo swainsoni"
"Circus cyaneus"              "Falco columbarius"           "Falco mexicanus"
"Falco perigrinus"            "Homo sapiens"   }
}
    \item{\code{Group}}{a factor with 3 levels (used for subOrigData function)
    \code{
"insect"
"bird"
"human"}}
    \item{\code{Source_quality}}{a factor with 2 levels \code{"known source" "presumed known source"}}
    \item{\code{Age_code}}{a factor with 6 levels \code{
 ""         "juvenile" "adult "   "adult"    "nestling" "chick".   }
  " " means unknown
 }
    \item{\code{Reference}}{a factor with 9 levels \code{ "Hobson et al. 1999"        "Hobson et al. 2012"        "Hobson and Wassenaar 1997"
 "Clark et al. 2006"         "Hobson et al., 2004"       "Lott and Smith, 2006"
"Thompson et al., 2010"     "Bowen et al., 2009"        "Ehleringer et al., 2008"  } }

  }
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
\seealso{
\code{\link[isOrigin]{origData}}
\code{\link[isOrigin]{subOrigData}}
}
\examples{
library(maptools)
data("knownOrig")
data("wrld_simpl")
plot("wrld_simpl")
points(knownOrig$Longitude,knownOrig$Latitude,col="red",cex=0.5)
}
\keyword{datasets}
