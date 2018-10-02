qtlRaster <- function(pdR, threshold, thresholdType){
  if(class(pdR) != "RasterLayer" & class(pdR) != "RasterStack" & class(pdR) != "RasterBrick"){
    stop("input probability density map (pdR) should be one of the following class: RasterLayer, RasterStack or RasterBrick")
  }
  if(class(threshold) != "numeric"){
    stop("threshold must be one number between 0 and 1 ")
  }
  if(length(threshold) != 1){
    stop("threshold must be one number between 0 and 1 ")
  }
  if(threshold < 0 | threshold > 1){
    stop("threshold must be one number between 0 and 1")
  }
  if(length(thresholdType) != 1){
    stop("thresholdType must be 1 or 2. See help page for further information")
  }
  if(thresholdType != 1 & thresholdType != 2){
    stop("thresholdType must be 1 or 2. See help page for further information")
  }
  result <- pdR
  if(thresholdType == 1){
    for(i in 1:nlayers(pdR)){
      pdR.values <- na.omit(getValues(pdR[[i]]))
      pdR.values <- sort(pdR.values)
      k <- length(pdR.values)
      left <- 1
      right <-  k
      while((right-left) > 2){
        start <- round((left+right)/2)
        sum <- sum(pdR.values[start:k])
        if(sum > threshold){
          left <- start
        }
        if(sum < threshold){
          right <- start
        }
      }
      result[[i]] <- pdR[[i]] > pdR.values[start]
    }
  }
  if(thresholdType == 2){
    for(i in 1:nlayers(pdR)){
      pdR.values <- na.omit(getValues(pdR[[i]]))
      k <- length(pdR.values)
      cut <- sort(pdR.values)[round((1-threshold)*k)]
      result[[i]] <- pdR[[i]]>cut
    }
  }
  names(result) <- names(pdR)
  print(plot(result))
  return(result)
}
