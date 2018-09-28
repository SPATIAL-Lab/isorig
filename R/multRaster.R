multRaster <- function(pdR){
  if(class(pdR) != "RasterStack"){
    stop("input probability density map (pdR) should be a RasterStack")
  }
  n <- nlayers(pdR)
  result <- pdR[[1]]*pdR[[2]]
  if(n > 2){
    for(i in 3:n){
      result <- result*pdR[[i]]
    }
  }
  result <- result/cellStats(result,sum)
  print(plot(result))
  return(result)
}
