oddsRatio <- function(pdR, inputP){
  if(class(pdR) != "RasterLayer" & class(pdR) != "RasterStack"){
    stop("input probability density map (pdR) should be either RasterLayer or RasterStack")
  }
  if(class(inputP) != "data.frame" & class(inputP) != "SpatialPolygons"){
    stop("inputP should be either data.frame or SpatialPolygons")
  }

  if(class(inputP) == "data.frame"){
    n <- nrow(inputP)
    extrVals <- extract(pdR, inputP)
    if(n == 1){
      result1 <- (extrVals/(1-extrVals))/(maxValue(pdR)/(1-maxValue(pdR)))
    }
    else if(n == 2){
      if(class(pdR) == "RasterStack"){
        result1 <- (extrVals[1,]/(1-extrVals[1,]))/(extrVals[2,]/(1-extrVals[2,]))
      } else {
        result1 <- (extrVals[1]/(1-extrVals[1]))/(extrVals[2]/(1-extrVals[2]))
      }
    }
    else{
      stop("input points (inputP) should be one or two points with longitude and latitude")
    }
    result2 <- data.frame(ratioToMax = extrVals/maxValue(pdR), ratioToMin = extrVals/minValue(pdR))
    result <- list(oddsRatio = result1, ratioToMaxMin = result2)
  }

  if(class(inputP) == "SpatialPolygons"){
    if(dimensions(inputP) != 2){
      stop("If inputP are SpatialPolygons, there should be SpatialPolygons with two polygons")
    }
    extrVals <- extract(pdR, inputP)
    extrVals.p1 <- colSums(extrVals[[1]])
    extrVals.p2 <- colSums(extrVals[[2]])
    result1 <- (extrVals.p1/(1-extrVals.p1))/(extrVals.p2/(1-extrVals.p2))
    result2 <- ncell(crop(pdR, inputP[1]))/ncell(crop(pdR, inputP[2]))
    result <- list(oddsRatio = result1, polygonCellRatio = result2)
  }

  return(result)
}
