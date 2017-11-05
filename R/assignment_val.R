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
# second column should be latitude, third column should be longitude and furth column should be
# tissue d2H values of the individuals for which the assignments will be made
# 3. d2Hcol = column number in the assign_table with the d2H tissue values
# 4. IDcol = column number with individual identifiers
#
# Return: all the assignment in one raster file

assignment_val <- function(rescaled_raster, unknown, d2Hcol, IDcol) {
  error <- rescaled_raster$sd
  if (class(unknown) == "character") {
    if (substr(unknown, nchar(unknown) - 3, nchar(unknown)) != ".csv")
      stop("Please use .csv file for unknown. Details refer to help page of this function")
    data <- read.table(assign_table, sep = ",", header = T)
  }
  if (class(unknown) == "data.frame")
    data <- unknown

  n <- length(data[, d2Hcol])
  dir.create("output")

  for (i in 1:n) {
    indv.data <- data[i, ]
    indv.id <- indv.data[1, IDcol]
    assign <- 1/sqrt((2 * pi * error^2)) * exp(-1 * (indv.data[1, d2Hcol] -
                                                       rescaled_raster$mean)^2/(2 * error^2))
    assign_norm <- assign/cellStats(assign, "sum")  #normalize so all pixels sum to 1
    filename <- paste("output/", indv.id, ".like", ".tif", sep = "")
    writeRaster(assign_norm, filename = filename, format = "GTiff",
                overwrite = TRUE)
  }

  result <- NULL
  dir.create("output_pdf")
  pdf("./output_pdf/output.pdf", width = 20, height = 40)
  # compare prediction and known origin
  par(mfrow = c(ceiling(length(unknown[, 1])/5), 5))
  for (i in 1:length(unknown[, 1])) {
    a <- raster(paste("output/", unknown[i, IDcol], ".like.tif", sep = ""))
    plot(a)
    text(unknown$Longitude[i], unknown$Latitude[i], "+", cex = 2)

    if (i == 1)
      result <- a else result <- stack(result, a)

  }
  dev.off()

  return(result)
}
