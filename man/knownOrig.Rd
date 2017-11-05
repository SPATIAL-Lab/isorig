\name{knownOrig}
\alias{knownOrig}
\docType{data}
\title{
Isotope value for tissue of known origin
}
\description{
This data consisits of hydrogen isotope for 2178 samples including human hair, insect and different species of birds, with location information.
}
\usage{data("knownOrig")}
\format{
  A data frame with 2178 observations on the following 10 variables.
  \describe{
    \item{\code{ID}}{a numeric vector}
    \item{\code{Longitude}}{a numeric vector}
    \item{\code{Latitude}}{a numeric vector}
    \item{\code{d2H}}{a numeric vector}
    \item{\code{Taxon}}{a numeric vector with 87 levels ([number] is used for getOrigData function) \code{[1]Monarch
[2]AMRE
[3]AMRO
[4]AUWA
[5]BCCH
[6]BEWR
[7]CARW
[8]CHSP
[9]COYE
[10]CSWA
[11]ETTI
[12]GWWA
[13]HETH
[14]HOWA
[15]KEWA
[16]MGWA
[17]NAWA
[18]NOCA
[19]OCWA
[20]ORJU
[21]OVEN
[22]REVI
[23]SOSP
[24]SWTH
[25]VEER
[26]WEVI
[27]WIWA
[28]WOTH
[29]YBCH
[30]YWAR
[31]ABTO
[32]BWWA
[33]HOSP
[34]RCSP
[35]CACH
[36]HOWR
[37]GRCA
[38]BAWW
[39]LOSH
[40]SPPI
[41]RUBL
[42].
[43]LEFL
[44]LESC
[45]TEWA
[46]Blue Tit
[47]House Sparrow
[48]Pheasant
[49]Red Grouse
[50]Black Grouse
[51]Black Woodpecker
[52]Warbler sp.
[53]Tit sp.
[54]European Serin
[55]Lapwing
[56]Carrion Crow
[57]Eurasian Blackbird
[58]Thrush sp.
[59]Carrion Crow
[60]Jackdaw
[61]Common Wood Pigeon
[62]Turtle Dove
[63]Hazel Grouse
[64]Partridge
[65]Mallard
[66]Canada Goose
[67]Rock Pigeon
[68]Eurasian Curlew
[69]Fieldfare
[70]Redwing
[71]Willow Grouse
[72]Grey Partridge
[73]Song Thrush
[74]Chaffinch
[75]RL
[76]SS
[77]AK
[78]NG
[79]CH
[80]RT
[81]BW
[82]SW
[83]NH
[84]ML
[85]PR
[86]PG
[87]hair} }
    \item{\code{Group}}{a factor with 3 levels ([number] is used for getOrigData function)
    \code{
[1]insect
[2]bird
[3]human}}
    \item{\code{Source_quality}}{a factor with levels \code{[1]known source [2]presumed known source}}
    \item{\code{Age_code}}{a factor with levels \code{
[1]NA
[2]AHY
[3]SY
[4]ASY
[5]HY
[6]unknown
[7]juvenile
[8]adult
[9]adult
[10]pull at nest
[11]flying juvenile
[12]chick
[13]fledging}}
    \item{\code{Reference}}{a factor with levels \code{[1]Hobson et al., 1999
[2]UCLA
[3]MB GWWA
[4]Samuel Hache and Erin Bayne
[5]Amy Chabot LOSH
[6]Katherine Brewster
[7]Rusty Blackbird
[8]Hobson and Wassenaar 1997
[9]Clark et al. 2006
[10]Hobson et al., 2004
[11]Lott and Smith, 2006
[12]Thompson et al., 2010
[13]Bowen et al., 2009
[14]Ehleringer et al., 2008} }
    \item{\code{Note}}{a factor with levels \code{[1]This is treated as Redwing in Hobson et al., 2004} \code{ [2]This is treated as wood pigeon in Hobson et al., 2004}}
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
\code{\link[isorig]{origData}}
\code{\link[isorig]{subOrigData}}
}
\examples{
data(knownOrig)
data(wrld_simpl)
plot(wrld_simpl)
points(knownOrig$Longitude,knownOrig$Latitude,col="red",cex=0.5)
}
\keyword{datasets}
