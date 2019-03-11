QA_2isotopes <- function(isoscape1, isoscape2, known1, known2, valiStation, valiTime, setSeed = T){
  if (class(isoscape1) == "RasterStack" | class(isoscape1) == "RasterBrick") {
    if (is.na(proj4string(isoscape1))){
      stop("isoscape1 must have coord. ref.")
    } else {
      isoscape1 <- projectRaster(isoscape1, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    }
  } else {
    stop("isoscape1 should be a RasterStack or RastrBrick")
  }
  if (class(isoscape2) == "RasterStack" | class(isoscape2) == "RasterBrick") {
    if (is.na(proj4string(isoscape2))){
      stop("isoscape must have coord. ref.")
    } else {
      isoscape2 <- projectRaster(isoscape2, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    }
  } else {
    stop("isoscape should be a RasterStack or RastrBrick")
  }
  if(!valiStation<nrow(known1)){
    stop("valiStation must be smaller than the number of known-origin stations in known1 or known2")
  }
  
  isoscape2 <- setValues(isoscape1, extract(isoscape2, coordinates(isoscape1)))
  if(setSeed == T){
    set.seed(100)
  }
  rowLength <- nrow(known1)
  val_numbers <- sort(sample(1:rowLength,valiStation,replace = F))
  for (i in 1:99){
    val_numbers <- rbind(val_numbers,sort(sample(1:rowLength,valiStation,replace = F)))
  }
  
  stationNum4model <- rowLength - valiStation
  
  prption_byProb <- matrix(0, valiTime, 99) # accuracy by checking top percentage by cumulative prob. 
  prption_byArea <- matrix(0, valiTime, 99) # accuracy by checking top percentage by area 
  pd_bird_val <- matrix(0, valiTime, valiStation) # pd value for each validation location
  precision <- list() # precision
  for (i in 1:valiTime){
    bird_val1 <- known1[val_numbers[i,],]
    bird_model1 <- known1[-val_numbers[i,],]
    bird_model_iso1 <- SpatialPointsDataFrame(cbind(bird_model1[,1],bird_model1[,2]), 
                                              as.data.frame(bird_model1[,3]))
    bird_val_iso1 <- SpatialPointsDataFrame(cbind(bird_val1[,1],bird_val1[,2]), 
                                            as.data.frame(bird_val1[,3]))
    crs(bird_model_iso1) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    crs(bird_val_iso1) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    rescale1 <- isOrigin::calRaster(bird_model_iso1, isoscape1, 
                                    sdMethod = 1, genplot = F, savePDF = F)
    pd1 <- isOrigin::pdRaster(rescale1$isoscape.rescale, 
                             data.frame(rownames(bird_val1),bird_val1[,3]), genplot = F, 
                             saveFile = F)
    
    bird_val2 <- known2[val_numbers[i,],]
    bird_model2 <- known2[-val_numbers[i,],]
    bird_model_iso2 <- SpatialPointsDataFrame(cbind(bird_model2[,1],bird_model2[,2]), 
                                              as.data.frame(bird_model2[,3]))
    bird_val_iso2 <- SpatialPointsDataFrame(cbind(bird_val2[,1],bird_val2[,2]), 
                                            as.data.frame(bird_val2[,3]))
    crs(bird_model_iso2) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    crs(bird_val_iso2) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    rescale2 <- isOrigin::calRaster(bird_model_iso2, isoscape2, 
                                    sdMethod = 1, genplot = F, savePDF = F)
    pd2 <- isOrigin::pdRaster(rescale2$isoscape.rescale, 
                              data.frame(rownames(bird_val2),bird_val2[,3]), genplot = F, 
                              saveFile = F)
    pd <- pd1
    for(j in 1:nlayers(pd)){
      pd[[j]] <- jointP(stack(pd1[[j]], pd2[[j]]))
    }
    
    # pd value for each validation location
    for(m in 1:nlayers(pd)){
      pd_bird_val[i, m] <- extract(pd[[m]], data.frame(bird_val1[,1],bird_val1[,2])[m,])
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
  
  result <- list(val_numbers, pd_bird_val, prption_byArea, prption_byProb, precision, random_prob_density)
  names(result) <- c("val_numbers", "pd_bird_val", "prption_byArea", "prption_byProb", "precision", "random_prob_density")
  class(result) <- "QA"
  return(result)
}