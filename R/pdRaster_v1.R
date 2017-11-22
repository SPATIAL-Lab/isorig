pdRaster <- function(rescaled.mean, rescaled.sd, unknown, mask = NULL, genplot=T) {

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
      rescaled.mean <- crop(rescaled.mean, mask)
      rescaled.sd <- crop(rescaled.sd, mask)
    } else {
      stop("mask should be a SpatialPolygonsDataFrame")
    }
  }
  if (class(unknown) != "data.frame") {
    stop("unknown should be a data.frame, see help page of pdRaster function")
  }
  if (class(unknown) == "data.frame"){
    data <- unknown
  }
  n <- length(data[, 2])
  dir.create("output")
  dir.create("output/pdRaster_Gtif")
  errorV <- getValues(rescaled.sd)
  meanV <- getValues(rescaled.mean)
  result <- NULL
  for (i in 1:n) {
    indv.data <- data[i, ]
    indv.id <- indv.data[1, 1]
    assign <- dnorm(indv.data[i, 2], mean = meanV, sd = errorV)
    assign.norm <- assign/sum(assign[!is.na(assign)])
    assign.norm <- setValues(rescaled.mean,assign.norm)
    if (genplot ==T){
      par(mfrow=c(1,1))
      plot(assign.norm)
    }
    if (i == 1){
      result <- assign.norm
    } else {
      result <- stack(result, assign.norm)
    }
    filename <- paste("output/pdRaster_Gtif/", indv.id, ".like", ".tif", sep = "")
    writeRaster(assign.norm, filename = filename, format = "GTiff",
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
