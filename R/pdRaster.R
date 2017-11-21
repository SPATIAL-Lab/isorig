###### pdRaster FUNCTION ######
####################################################################################
## pdRaster.R
## 09/21/2017, SLC, Chao Ma
## AIMS: 1. assign origin of animal based on tissue isotope and rescaled isoscape
##
## 11/16/2017: change name from assign to pdRaster
####################################################################################

# This function uses the code from Vander Zanden et al., 2014,
# which uses the likelihood term (Equation S2 in
# Appendix 1 of Vander Zanden et al., 2014) to determine the
# probability that an individual sample was from a particular
# geographic location and writes an tif file to a output directory
#
# Args:
# 1. rescaled_raster = the tissue-specific d2H raster (mean and sd) created in the rescale
# function.
# IsoMAP output.  This is the precip component of the variance term.
# 2. unknown = this a csv filename (and directory, if applicable). First column should be ID,
# second column should be tissue d2H values of the individuals for which the assignments will be made
# 3. mask: a SpatialPolygonsDataFrame that is used to constrain the investigated area. If this is not provided, default of whole world is used.
#
# Return: all the assignment in one raster file

pdRaster <- function(rescaled_raster, unknown, mask = NULL, genplot=T) {

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

    if (!is.null(mask)) {
    if (class(mask) == "SpatialPolygonsDataFrame") {
      rescaled_raster <- crop(rescaled_raster, mask)
    } else {
      stop("mask should be a SpatialPolygonsDataFrame")
    }
  }
  errorV <- getValues(rescaled_raster$sd)
  if (class(unknown) != "data.frame") {
    stop("unknown should be a data.frame, see help page of pdRaster function")
  }
  if (class(unknown) == "data.frame"){
    data <- unknown
  }

  n <- length(data[, 2])
  dir.create("output")
  dir.create("output/pdRaster_Gtif")

  meanV <- getValues(rescaled_raster$mean)
  result <- NULL
  for (i in 1:n) {
    indv.data <- data[i, ]
    indv.id <- indv.data[1, 1]
    assign <- 1/sqrt((2 * pi * errorV^2)) * exp(-1 * (indv.data[1, 2] -
                                                       meanV)^2/(2 * errorV^2))
    assign_norm <- assign/sum(assign[!is.na(assign)])  #normalize so all pixels sum to 1
    assign_norm <- setValues(rescaled_raster$mean,assign_norm)
    if (genplot ==T){
      par(mfrow=c(1,1))
      plot(assign_norm)
    }
    if (i == 1){
      result <- assign_norm
    } else {
      result <- stack(result, assign_norm)
    }
    filename <- paste("output/pdRaster_Gtif/", indv.id, ".like", ".tif", sep = "")
    writeRaster(assign_norm, filename = filename, format = "GTiff",
                overwrite = TRUE)
  }


  if (n > 5){
    pdf("./output/output_pdRaster.pdf", width = 20, height = 40)
    par(mfrow = c(ceiling(length(data.class()[, 1])/5), 5))
  } else {
    pdf("./output/output_pdRaster.pdf", width = 20, height = 20)
  }

  for (i in 1:n) {
    a <- raster(paste("output/pdRaster_Gtif/", data[i, 1], ".like.tif", sep = ""))
    plot(a)
    text(data[i,1], data[i,2], "+", cex = 2)

    if (i == 1){
      result <- a
    } else {
      result <- stack(result, a)
    }
  }
  dev.off()


  return(result)
}
