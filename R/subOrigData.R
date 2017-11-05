###### subOrigData FUNCTION ######
####################################################################################
## origData.R
## 11/01/2017, SLC, Chao Ma
## AIMS:
## query data from known origin tissue d2H from isorig pacakge, and apply mask if necessary.
####################################################################################
#
# Args:
#   taxon: one string or string vector, choosing from taxon in known origin data Details of these names can be seen in knownOrig
#   group: one string or a string vector, choosing from group names of known origin data. Details can be seen in knownOrig
#   mask : a SpatialPolygonsDataFrame that is used to constrain the investigated area. If this is not provided, default of whole world is used.
#
# Returns: known origin data from isorig package

subOrigData <- function(taxon = NULL, group = NULL, mask = NULL) {
  # load data knownOrig
  data("knownOrig")
  result <- NULL

  # define taxon
  if (!is.null(taxon)){
    if (all(taxon %in% unique(knownOrig$Taxon))) {
      result <- knownOrig[knownOrig$Taxon %in% taxon, ]
    } else {
      stop("taxon should be string or string vector given from Taxon column in knownOrig. Please see knownOrig help page!")
    }
  }

  # define group
  if (!is.null(group)){
    if (all(group %in% unique(knownOrig$Group))) {
      result <- knownOrig[knownOrig$Group %in% group, ]
    } else {
      stop("group should be string or string vector given from Group column in knownOrig. Please see knownOrig help page!")
    }
  }


  if (!is.null(taxon) && !is.null(group))
    stop("Please either choose taxon or group")

  if (!is.null(mask)) {
    if (class(mask) == "SpatialPolygonsDataFrame") {

      s <- SpatialPointsDataFrame(coords = cbind(result$Longitude,
                                                 result$Latitude), data = result, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      o <- over(s, mask)

      overlap <- result[!is.na(o), ]
    } else {
      stop("mask should be a SpatialPolygonsDataFrame")
    }

    if (length(overlap[, 1]) != 0) {
      s1 <- SpatialPointsDataFrame(coords = cbind(overlap$Longitude,
                                                  overlap$Latitude), data = overlap, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      plot(mask, axes = T)
      plot(s1, add = T, col = "red")
    } else {
      cat("No isotope data found in mask you choose!\n")
    }
    result <- overlap
  } else {
    require(maptools)
    data(wrld_simpl)
    plot(wrld_simpl, axes = T)
    points(result$Longitude, result$Latitude, col = "red", cex = 0.5)
  }
  cat(length(result[,1]),"data points are found\n")
  return(result)
}
