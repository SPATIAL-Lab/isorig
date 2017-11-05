##################################################################
## Total Uncertainty estimated from RESCALING FUNCTION
## This uncertainty stem from two sources:
## 1. precipitation isoscape
## 2. rescale function itself
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
# 1. orig = 1: input known origin as a csv file. 3 columns: 1st is longitude, 2nd is latitude, 3rd is tissue isotope value.
# The first row should be header of "longitude", "latitude", "value"
# If you have a shape file, you should select these 3 attributions and save it as a csv file with format mentioned above.
#    orig = 2: input known origin as a raster file. This raster file should have one layer of known origin isotope value.
#    orig = 3: known origin data is from this pacakge. If this is chosen, species and region should be provided.
# 2. file = , should be provided if orgin = 1 or 2, this is the filename (with directory, if applicable) of known origin.
# 3. species = "bird": bird
#    species = "bat": bat
#    species = "insect": insect
# 4. region = "World": whole world
#    region = "NA": North American
#    region = "EU": Europe
#    region = "Asia": Asia
# 5. precip = filename (with directory, if applicable) of prediction files from IsoMAP (NOTE: feature need to be added that users want to use their own IsoMAP)
# This must be a raster file or txt file of prediction from IsoMAP.
# 6. interpMethod = numeric. 1 or 2. This is method for extracting
# values from precipitation raster based on know origin position. If 1 values
# for the cell a point fall s in are returned. If 2 the returned values are
# interpolated from the values of the four nearest raster cells.
#
# 3.species and 4.region can always be polished as database in package increases. If in the future the database is getting too large to store in package, it can be stored online database


rescale <- function(orig, file, species, region, precip, interpMethod = 1)
{
  #load libraries
  if(require("raster")){
    print("raster is loaded correctly")
  } else {
    print("trying to install raster")
    install.packages("raster")
    if(require(lme4)){
      print("raster installed and loaded")
    } else {
      stop("could not install raster")
    }
  }

  if(require("rgdal")){
    print("rgdal is loaded correctly")
  } else {
    print("trying to install rgdal")
    install.packages("rgdal")
    if(require(lme4)){
      print("rgdal installed and loaded")
    } else {
      stop("could not install rgdal")
    }
  }
  if(require("sp")){
      print("sp is loaded correctly")
    } else {
      print("trying to install sp")
      install.packages("sp")
      if(require(lme4)){
        print("sp installed and loaded")
      } else {
        stop("could not install sp")
      }
    }
  if(require("varhandle")){
    print("varhandle is loaded correctly")
  } else {
    print("trying to install varhandle")
    install.packages("varhandle")
    if(require(lme4)){
      print("varhandle installed and loaded")
    } else {
      stop("could not install varhandle")
    }
  }

  #library(raster)
  #library(rgdal)
  #library(sp)
  #library(varhandle)

  # d2h residuals from lab test. Data is given by Mike Wunder
  # Analytical calibration residuals. They are described in Wunder and Norris
  # (2008). They were collected from 150 different runs, each of which had
  # between 2-5 replicates of two different standards.
  data("d2h.resids")

  # error messages
  if (!(orig %in% c(1,2,3)))
    stop("orig in function input must be one number of 1, 2, 3. Detail of this parameter can be seen in help page of this function")
  if (orig %in% c(1,2) && (substr(file, nchar(file)-3, nchar(file))!= ".csv") )
    stop("file name (with directory if applicable and end with .csv) must be provided for the function input of 'file'. Detail of this parameter can be seen in help page of this function")
  if (orig == 3 && !(species %in% c("bird","bat","insect") && region %in% c("World","NA","EU","Asia")))
    stop("species and region must be specified when orig is chosen as 3. Detail of this parameter can be seen in help page of this function")

  if (orig == 1)
  {
    dat <- read.table(file, nrows = 1, colClasses = "character",
                      sep = ",")
    #titles = is.na(suppressWarnings(as.numeric(dat[1, 1])))
    if (!is.na(suppressWarnings(as.numeric(dat[1, 1])))) {
      cat("\\n * No column titles/headers detected in known origin csv file\\n")
      calibration <- read.table(filen, header = F, sep = ",")
    }
    if (is.na(suppressWarnings(as.numeric(dat[1, 1])))) {
      cat("\\n * Column titles/headers detected in known origin csv file\\n")
      calibration <- read.table(filen, header = T, sep = ",")
    }
  }
  if (orig == 2)
  {
    calibration = rasterToPoints(raster(file))
  }
  if (orig == 3)
  {
    if (species == "insect" && region == "NA")
    {
      data("insect_na")
      calibration = insect_na
      print("Known origin data of d2h of insect in North American is selected")
    }
    if (species == "bird" && region == "NA")
    {
      data("bird_na")
      calibration = bird_na
      print("Known origin data of d2h of insect in North American is selected")
    }
  }


  nSample <- nrow(calibration)
  # create two blank vector for storing the isotope of tissue and precip
  tissue.iso <- vector("numeric", length = nSample)
  precip.iso <- vector("numeric", length = nSample)

  # create a raster that has similar length with raster obtained from
  # prediction_i.txt NOTE that the raster obtained from prediction_i.txt
  # using following method has less cells that raster from IsoMAP...
  if (substr(precip, nchar(precip)-3, nchar(precip)) == ".txt")
  {
    pts <- read.table(precip, header = T)
    ptskrig <- cbind(pts[, 3:4], pts[, length(pts[1, ]) - 1])
    prediction <- rasterFromXYZ(ptskrig)
  }else
  {
    prediction = raster(precip)
    ncells <- ncell(prediction)
  }

  #resamaple (bootstrap with raplacement) calibration (known origin) data
  t <- t(calibration)
  t <- as.data.frame(t)
  t1 <- sample(t, replace = TRUE)
  t <- t(t1)
  calibration_resamp <- cbind(as.numeric(t[,1]),as.numeric(t[,2]),as.numeric(t[,3]))

  # random pick lab uncertainties
  iso.resids.nSample <- sample(d2h.resids, nSample)

  # assgin tissue and precipitation isotopic values to the positions of know origin sites
  tissue.iso <- as.numeric(calibration_resamp[, 3]) + iso.resids.nSample
  if (interpMethod == 1)
  {
    precip.iso <- raster::extract(prediction, calibration_resamp[,
                                                                 1:2], method = "simple")
  } else if (interpMethod == 2)
  {
    precip.iso <- raster::extract(prediction, calibration_resamp[,
                                                                 1:2], method = "bilinear")
  }else
     {
       stop("interpMethod should be either 1 or 2")
     }

  # linear regression between known origin and precipitation data
  lmResult <- lm(tissue.iso ~ precip.iso)
  intercept <- as.numeric(coef(lmResult)[1])
  slope <- as.numeric(coef(lmResult)[2])
  precip.rescale <- pts[, length(pts[1, ]) - 1] * slope + intercept
  precip.rescale.XYZ <- cbind(pts[, 3:4], precip.rescale)
  precip.rescale <- rasterFromXYZ(cbind(pts[, 3:4], precip.rescale))

  # calculate mean and sd of rescaled raster, sd of rescaled raster is
  # the error of estimation from the linear regression above.
  precip.rescale.XYZ[,3]=summary(lmResult)$sigma #Same sd is assigned to each cell.
  precip.rescale.sd <- rasterFromXYZ(precip.rescale.XYZ)
  precip.rescale <- stack(precip.rescale, precip.rescale.sd)
  names(precip.rescale) <- c("mean", "sd")
  plot(precip.rescale)

  # return stack raster of mean and sd of resale map
  return(precip.rescale)
}


