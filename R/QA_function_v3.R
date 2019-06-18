QA <- function(isoscape, known, valiStation, valiTime, setSeed = T){

  #check that isoscape is valid and has defined CRS
  if (class(isoscape) == "RasterStack" | class(isoscape) == "RasterBrick") {
    if (is.na(proj4string(isoscape))) {
      stop("isoscape must have valid coordinate reference system")
    }
  } else {
    stop("isoscape should be a RasterStack or RastrBrick")
  }

  #check that known is valid and has defined, correct CRS
  if (class(known) != "SpatialPointsDataFrame") {
    stop("known should be a SpatialPointsDataFrame, see help page of calRaster function")
  } else if (is.na(proj4string(known))) {
    stop("known must have valid coordinate reference system")
  } else if(proj4string(known) != proj4string(isoscape)){
    stop("known must have same coordinate reference system as isoscape")
  } else if(ncol(known@data) != 1){
    stop("known must include a 1-column data frame containing only the isotope values")
  }

  if(!valiStation<nrow(known)){
    stop("valiStation must be smaller than the number of known-origin stations in known")
  }
  if(setSeed == T){
    set.seed(100)
  }
  rowLength <- nrow(known)
  val_stations <- sort(sample(1:rowLength,valiStation,replace = F))
  for (i in 1:99){
    val_stations <- rbind(val_stations,sort(sample(1:rowLength,valiStation,replace = F)))
  }

  stationNum4model <- rowLength - valiStation
  prption_byProb <- matrix(0, valiTime, 99) # accuracy by checking top percentage by cumulative prob.
  prption_byArea <- matrix(0, valiTime, 99) # accuracy by checking top percentage by area
  pd_bird_val <- matrix(0, valiTime, valiStation) # pd value for each validation location
  precision <- list() # precision
  known <- as.data.frame(known)
  for (i in 1:valiTime){
    bird_val <- known[val_stations[i,],]
    bird_model <- known[-val_stations[i,],]
    bird_model_iso <- SpatialPointsDataFrame(cbind(bird_model[,1],bird_model[,2]), as.data.frame(bird_model[,3]))
    bird_val_iso <- SpatialPointsDataFrame(cbind(bird_val[,1],bird_val[,2]), as.data.frame(bird_val[,3]))
    crs(bird_model_iso) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    crs(bird_val_iso) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    rescale <- isOrigin::calRaster(bird_model_iso, isoscape, sdMethod = 1, genplot = F, savePDF = F)
    pd <- isOrigin::pdRaster(rescale$isoscape.rescale, data.frame(rownames(bird_val),bird_val[,3]), genplot = F, saveFile = F)

    # pd value for each validation location
    for(m in 1:nlayers(pd)){
      pd_bird_val[i, m] <- extract(pd[[m]], data.frame(bird_val[,1],bird_val[,2])[m,])
    }

    xx <- seq(0.01, 0.99, 0.01) ## 0.01 to 0.99

    # total area
    Tarea <- length(na.omit(pd[[1]][]))

    # accuracy and precision by checking top percentage by cumulative prob.
    precision[[i]] <- matrix(0, 99, valiStation) # precision
    for(j in xx){
      qtl <- isOrigin::qtlRaster(pd, threshold = j, pdf = F, thresholdType = 1,genplot = F)
      prption_byProb[i, j*100] <- 0
      for(k in 1:nlayers(qtl)){
        prption_byProb[i, j*100] <- prption_byProb[i, j*100] +
          raster::extract(qtl[[k]], data.frame(bird_val[,1],bird_val[,2])[k,])
        precision[[i]][j*100, k] <- sum(na.omit(qtl[[k]][]))/Tarea # precision
      }
    }

    # result is not right when j = 0.07 in above loop with no reason, so j = 0.07 is rerun
    j = 0.07
    qtl <- isOrigin::qtlRaster(pd, threshold = j, pdf = F, thresholdType = 1,genplot = F)
    prption_byProb[i, j*100] <- 0
    for(k in 1:nlayers(qtl)){
      prption_byProb[i, j*100] <- prption_byProb[i, j*100] +
        raster::extract(qtl[[k]], data.frame(bird_val[,1],bird_val[,2])[k,])
      precision[[i]][j*100, k] <- sum(na.omit(qtl[[k]][]))/Tarea # precision
    }
    for(k in 1:nlayers(qtl)){
      precision[[i]][j*100, k] <- sum(na.omit(qtl[[k]][]))/Tarea # precision
    }

    # accuracy by checking top percentage by cumulative area
    for(n in xx){
      qtl <- isOrigin::qtlRaster(pd, threshold = n, pdf = F, thresholdType = 2,genplot = F)
      prption_byArea[i, n*100] <- 0
      for(k in 1:nlayers(qtl)){
        prption_byArea[i, n*100] <- prption_byArea[i, n*100] +
          raster::extract(qtl[[k]], data.frame(bird_val[,1],bird_val[,2])[k,])
      }
    }
    # result is not right when n = 0.07 in above loop with no reason, so n = 0.07 is rerun
    n = 0.07
    qtl <- isOrigin::qtlRaster(pd, threshold = n, pdf = F, thresholdType = 2,genplot = F)
    prption_byArea[i, n*100] <- 0
    for(k in 1:nlayers(qtl)){
      prption_byArea[i, n*100] <- prption_byArea[i, n*100] +
        raster::extract(qtl[[k]], data.frame(bird_val[,1],bird_val[,2])[k,])
    }
  }

  random_prob_density=1/length(na.omit(getValues(isoscape[[1]])))

  result <- list(val_stations, pd_bird_val, prption_byArea, prption_byProb, precision, random_prob_density)
  names(result) <- c("val_stations", "pd_bird_val", "prption_byArea", "prption_byProb", "precision", "random_prob_density")
  class(result) <- "QA"
  return(result)
}
