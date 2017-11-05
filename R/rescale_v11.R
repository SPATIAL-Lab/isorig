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
# Features need to be added:
# 1. generic scale function
# 2. if simulate isoscapes. Then using a loop. if user provide isoscapes then only run resacle and assginment once.
# 3. For assignment:
#    a. Regions with p = 95% CL of origin, where p can be changed.
#    b. user can input a vector of isotope value of tissue then assign.
#    c. multi isotopic assignment.
#+++++++++++++++++++++++++#


# This function conducts simulated regressions between tissue.iso and
# precip.iso. tissue.iso is resampled from known origin data + lab measuring
# uncertainty; precip.iso is from random precip map based on mean and SD. The
# output contains rescaled precipitaiton model with mean and SD. Function
# includes:
#
# 1. orig: known origin. Two senarios: (1)If user provided, this should be a filename (with directory, if applicable) of known origin.
#   It should be as a csv file with 3 columns: 1st is longitude, 2nd is latitude, 3rd is tissue isotope value.
#   The first row should be header of "Longitude", "Latitude", "d2H".
#   If you have a shape file, you should select these 3 attributions and save it as a csv file with format mentioned above.
#   (2) known origin data provided by this pacakge. This should be output from function: subOrigData, with format of data.frame
# 2. precip: raster of isoscape
# 3. mask  : a SpatialPolygonsDataFrame that is used to constrain the investigated area. If this is not provided, default of whole world is used.
# 4. interpMethod = numeric. 1 or 2. This is method for extracting
# values from precipitation raster based on know origin position. If 1 values
# for the cell a point falls in are returned. If 2 (default) the returned values are
# bilinear interpolated from the values of the four nearest raster cells.


rescale <- function(orig, precip, mask = NULL, interpMethod = 2) {
  # load libraries
  if (require("raster")) {
    print("raster is loaded correctly")
  } else {
    print("trying to install raster")
    install.packages("raster")
    if (require(lme4)) {
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
    if (require(lme4)) {
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
    if (require(lme4)) {
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
    if (require(lme4)) {
      print("varhandle installed and loaded")
    } else {
      stop("could not install varhandle")
    }
  }

  # d2h residuals from lab test. Data is given by Mike Wunder Analytical
  # calibration residuals. They are described in Wunder and Norris
  # (2008). They were collected from 150 different runs, each of which
  # had between 2-5 replicates of two different standards.
  data("d2h.resids")

  if (class(orig) == "character") {
    if (substr(orig, nchar(orig) - 3, nchar(orig)) == ".csv") {
      dat <- read.table(orig, nrows = 1, colClasses = "character",
                        sep = ",")
      # titles = is.na(suppressWarnings(as.numeric(dat[1, 1])))
      if (!is.na(suppressWarnings(as.numeric(dat[1, 1])))) {
        cat("\\n * No column titles/headers detected in known origin csv file\\n")
        calibration <- read.table(orig, header = F, sep = ",")
      }
      if (is.na(suppressWarnings(as.numeric(dat[1, 1])))) {
        cat("\\n * Column titles/headers detected in known origin csv file\\n")
        calibration <- read.table(orig, header = T, sep = ",")
      }
    }
  }


  if (class(orig) == "data.frame") {
    cat("known origin data is provided by isorig package \n")
    calibration <- orig[, 2:4]
  }

  # apply mask if provided
  if (!is.null(mask)) {
    if (class(mask) == "SpatialPolygonsDataFrame") {

      s <- SpatialPointsDataFrame(coords = calibration[,1:2],
                                                  data = calibration, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      o <- over(s, mask)

      calibration <- calibration[!is.na(o), ]
    } else {
      stop("mask should be a SpatialPolygonsDataFrame")
    }
  }

  nSample <- nrow(calibration)
  # create two blank vector for storing the isotope of tissue and precip
  tissue.iso <- vector("numeric", length = nSample)
  precip.iso <- vector("numeric", length = nSample)

  # create a raster that has similar length with raster obtained from
  # prediction_i.txt NOTE that the raster obtained from prediction_i.txt
  # using following method has less cells than raster from IsoMAP...

  if (class(precip) == "RasterLayer") {
    prediction <- precip
  } else {
    stop("precip should be a RasterLayer")
  }
  ncells <- ncell(prediction)

  # random pick lab uncertainties
  iso.resids.nSample <- sample(d2h.resids, nSample)

  # assgin tissue and precipitation isotopic values to the positions of
  # know origin sites
  tissue.iso <- as.numeric(calibration[, 3]) + iso.resids.nSample
  if (interpMethod == 1) {
    precip.iso <- raster::extract(prediction, calibration[, 1:2], method = "simple")
  } else if (interpMethod == 2) {
    precip.iso <- raster::extract(prediction, calibration[, 1:2], method = "bilinear")
  } else {
    stop("interpMethod should be either 1 or 2")
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
  values(sd) <- summary(lmResult)$sigma  #Same sd is assigned to each cell.
  values(sd)[is.na(values(precip))] <- NA
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


