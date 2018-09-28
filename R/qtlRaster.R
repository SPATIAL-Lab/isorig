qtlRaster <- function(pdR, threshold, thresholdType){
  if(class(pdR) != "RasterLayer" & class(pdR) != "RasterStack"){
    stop("input probability density map (pdR) should be either RasterLayer or RasterStack")
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
      for(j in 9999:1){
        a <- cellStats(overlay(pdR[[i]]>quantile(pdR, j/10000)[i], pdR[[i]], fun=function(x,y){return(x*y)}),sum)- threshold
        if(a>0){break}
        result[[i]] <- pdR[[i]]>quantile(pdR, j/10000)[i]
      }
    }
  }
  if(thresholdType == 2){
    totCell <- cellStats(pdR[[1]]>0, sum)
    for(i in 1:nlayers(pdR)){
      for(j in 1:1000){
        if(abs(cellStats(pdR[[i]]>quantile(pdR, j/1000)[i], sum) - totCell*threshold) < 10){
          break
        }
      result[[i]] <- pdR[[i]]>quantile(pdR, j/1000)[i]
      }
    }
  }
  print(plot(result))
  return(result)
}
