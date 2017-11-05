###### assignment FUNCTION ######
####################################################################################
## assignment.R
## 09/21/2017, SLC, Chao Ma
## AIMS: 1. assign origin of animal based on tissue isotope and rescaled isoscape
##
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

assignment <- function(rescaled_raster, unknown, mask = NULL, genplot=T) {


  if (!is.null(mask)) {
    if (class(mask) == "SpatialPolygonsDataFrame") {
      rescaled_raster <- crop(rescaled_raster, mask)
    } else {
      stop("mask should be a SpatialPolygonsDataFrame")
    }
  }
  error <- rescaled_raster$sd
  if (class(unknown) == "character") {
    if (substr(unknown, nchar(unknown) - 3, nchar(unknown)) != ".csv") {
      stop("Please use .csv file for unknown. Details refer to help page of assignment function")
    }
    data <- read.table(unknown, sep = ",", header = T)
  }
  if (class(unknown) == "data.frame"){
    data <- unknown
  }

  n <- length(data[, 2])
  dir.create("output")
  dir.create("output/assignment_Gtif")

  for (i in 1:n) {
    indv.data <- data[i, ]
    indv.id <- indv.data[1, 1]
    assign <- 1/sqrt((2 * pi * error^2)) * exp(-1 * (indv.data[1, 2] -
                                                       rescaled_raster$mean)^2/(2 * error^2))
    assign_norm <- assign/cellStats(assign, "sum")  #normalize so all pixels sum to 1



    filename <- paste("output/assignment_Gtif/", indv.id, ".like", ".tif", sep = "")
    writeRaster(assign_norm, filename = filename, format = "GTiff",
                overwrite = TRUE)
  }

  result <- NULL

  for (i in 1:n) {
    a <- raster(paste("output/assignment_Gtif/", data[i, 1], ".like.tif", sep = ""))
    if (genplot ==T){
      par(mfrow=c(1,1))
      plot(a)
    }

    if (i == 1){
      result <- a
    } else {
      result <- stack(result, a)
    }
  }


  if (n > 5){
    pdf("./output/output_assignment.pdf", width = 20, height = 40)
    par(mfrow = c(ceiling(length(data.class()[, 1])/5), 5))
  } else {
    pdf("./output/output_assignment.pdf", width = 20, height = 20)
  }

  for (i in 1:n) {
    a <- raster(paste("output/assignment_Gtif/", data[i, 1], ".like.tif", sep = ""))
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
