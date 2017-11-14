##################################################################
## Total Uncertainty estimated from RESCALING FUNCTION
## This uncertainty stem from two sources:
## 1. precipitation isoscape
## 2. rescale function itself
##
## This function select one resample of bootstraping known origin data.
##################################################################


#+++++++++++++++++++++++++#
#+++++Version updates+++++#
# March.9.2017: Sample from known origin data changed to use "sample" function
# with "replacement = TRUE" based on all individual sample data from all sites.
# And added lab uncertainties to tissue.iso
#
# May.3.2017: Boostrapping resample obs.txt in IsoMAP which is used to calculate
# the predictions. These predictions are used to calcualte mean and SD of
# precipitation isomap.
#
# Aug.10.2017-V7: krig bootstraped prediciton is done on IsoMAP. Download it
# from IsoMAP and read it through function of read.isomap
#
# Sep.19.2017-v9: no iteration for rescale
#
# Oct.4.2017-v10: Add different options of known origin input
#
# 10/19/2017-v11: modified according to new functions of subOrigData
#
# 11/06/2017: add feature of checking NA values extracted from raster in lm
#
#+++++++++++++++++++++++++#


rescale <- function(orig, precip, mask = NULL, interpMethod = 2, NA.value = NA, ignore.NA = T) {
  # load libraries
  if (require("raster")) {
    print("raster is loaded correctly")
  } else {
    print("trying to install raster")
    install.packages("raster")
    if (require("raster")) {
      print("raster installed and loaded")
    } else {
      stop("could not install raster")
    }
  }
  if (require("rgdal")) {
    print("rgdal is loaded correctly")
  } else {
    print("trying to install rgdal")
    install.packages("rgdal")
    if (require("rgdal")) {
      print("rgdal installed and loaded")
    } else {
      stop("could not install rgdal")
    }
  }
  if (require("sp")) {
    print("sp is loaded correctly")
  } else {
    print("trying to install sp")
    install.packages("sp")
    if (require("sp")) {
      print("sp installed and loaded")
    } else {
      stop("could not install sp")
    }
  }
  if (require("varhandle")) {
    print("varhandle is loaded correctly")
  } else {
    print("trying to install varhandle")
    install.packages("varhandle")
    if (require("varhandle")) {
      print("varhandle installed and loaded")
    } else {
      stop("could not install varhandle")
    }
  }
  if (require("ggplot2")) {
    print("ggplot2 is loaded correctly")
  } else {
    print("trying to install ggplot2")
    install.packages("ggplot2")
    if (require("ggplot2")) {
      print("ggplot2 installed and loaded")
    } else {
      stop("could not install ggplot2")
    }
  }

  if (class(orig) != "data.frame") {
    stop("orig should be a data.frame, see help page of rescale function")
  }

  # apply mask if provided
  if (!is.null(mask)) {
    if (class(mask) == "SpatialPolygonsDataFrame") {

      s <- SpatialPointsDataFrame(coords = orig[,1:2],
                                                  data = orig, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      o <- over(s, mask)

      orig <- orig[!is.na(o), ]
    } else {
      stop("mask should be a SpatialPolygonsDataFrame")
    }
  }

  nSample <- nrow(orig)
  # create two blank vector for storing the isotope of tissue and precip
  tissue.iso <- vector("numeric", length = nSample)
  precip.iso <- vector("numeric", length = nSample)
  # create a blank vector for storing the location where precip does  not have value
  null.iso <- NULL

  if (class(precip) == "RasterLayer") {
    prediction <- precip
  } else {
    stop("precip should be a RasterLayer")
  }
  ncells <- ncell(prediction)

  # assgin tissue and precipitation isotopic values to the positions of
  # know origin sites
  temp <- apply(orig,2,class)
  if (any(temp!="numeric")){
    stop("orig data should be all numeric. Please check that except header, orig table should only contains numbers.")
  }

  tissue.iso <- as.numeric(orig[, 3])
  if (interpMethod == 1) {
    precip.iso <- raster::extract(prediction, orig[, 1:2], method = "simple")
  } else if (interpMethod == 2) {
    precip.iso <- raster::extract(prediction, orig[, 1:2], method = "bilinear")
  } else {
    stop("interpMethod should be either 1 or 2")
  }

  # if there is no precip value at known origin location
  if (!is.na(NA.value)){
    values(precip)[values(precip)==NA.value] <- NA
  }
  if (any(is.na(precip.iso))){
    na <- which(is.na(precip.iso))
    cat("\n\n----------------------------------------------------------------\n")
    cat("Warning: NO data are found at following locations:\n")
    print(orig[na,1:2])
    if (!ignore.NA){
      stop ("Delete these data in known origin data or use a different isoscape that has values at these locations")
    }
    precip.iso <- precip.iso[-na]
    tissue.iso <- tissue.iso[-na]
  }


  # linear regression between known origin and precipitation data
  lmResult <- lm(tissue.iso ~ precip.iso)
  cat("\n\n---------------------------------------------------------------------------------\n")
  cat("rescale function uses linear regression model, the summary of this model is:\n")
  cat("---------------------------------------------------------------------------------\n")
  print(summary(lmResult))
  intercept <- as.numeric(coef(lmResult)[1])
  slope <- as.numeric(coef(lmResult)[2])
  precip.rescale <- precip * slope + intercept
  # plot linear regression result
  x = precip.iso
  y = tissue.iso
  xy = data.frame(x,y)
  equation <- function(x) {
    lm_coef <- list(a = round(coef(x)[1], digits = 2),
                    b = round(coef(x)[2], digits = 2),
                    r2 = round(summary(x)$r.squared, digits = 2));
    lm_eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,lm_coef)
    as.character(as.expression(lm_eq));
  }

  fit <- lm(y~x, data = xy)
  xmin = max(x)-0.4*(diff(range(x)))
  xmax = max(x)
  ymin = min(y)-0.25*diff(range(y))
  ymax = min(y)-0.1*diff(range(y))

  p11 <- ggplot(xy, aes(x=precip.iso, y=tissue.iso)) + geom_point(shape=1) + geom_smooth(method=lm, se=T) +
    ggtitle("Rescale regression model")+theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(name = "Precipitation isotope") +
    scale_y_continuous(name = "Tissue isotope") +
    annotate("rect", xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill="white", colour="red") +
    annotate("text", x = (xmax+xmin)/2, y = (ymax+ymin)/2,label = equation(fit), parse = TRUE)
  print(p11)

  # calculate mean and sd of rescaled raster, sd of rescaled raster is
  # the error of estimation from the linear regression above.
  sd <- precip
  values(sd) <- summary(lmResult)$sigma  # Same sd is assigned to each cell.
  values(sd)[is.na(values(precip))] <- NA # replace all no data with NA
  precip.rescale <- stack(precip.rescale, sd)
  names(precip.rescale) <- c("mean", "sd")
  if (!is.null(mask)) {
  precip.rescale <- crop(precip.rescale,mask)
  }
  plot(precip.rescale, main = c("rescale mean", "rescale sd"))
  # return stack raster of mean and sd of resale map
  dir.create("output")
  pdf("./output/rescale.result.pdf",width = 8,height = 4)
  print(p11)
  plot(precip.rescale, main = c("rescale mean", "rescale sd"))
  dev.off()

  names(xy) = c("precip.iso","tissue.iso")
  result = list ("precip.rescale" = precip.rescale, "lm.data" = xy, "lm.model" = lmResult)
  return(result)
}


