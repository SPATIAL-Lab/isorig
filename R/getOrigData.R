###### origData FUNCTION ######
####################################################################################
## origData.R
## 11/01/2017, SLC, Chao Ma
## AIMS:
## get data from known origin tissue d2H from isOrigin pacakge.
## and print all catagories and unique value for each catagory
####################################################################################
#
# Args: no args
#
# Returns: known origin data from isOrigin package

getOrigData <- function() {
  data(knownOrig)  # load knownOrig data from isOrigin package
  print("Known origin data has following catagories:")
  colnames(knownOrig)

  # summary(knownOrig)
  for (i in 1:length(knownOrig[1, ])) {
    if (class(knownOrig[, i]) == "numeric" || class(knownOrig[, i]) ==
        "integer") {
      cat(paste("Category:", names(knownOrig)[i]), "\n")
      cat("max , min\n")
      cat(max(knownOrig[, i]), ",", min(knownOrig[, i]), "\n")
      cat("\n")
    }
    if (class(knownOrig[, i]) == "factor") {
      cat(paste("Category:", names(knownOrig)[i]), "\n")
      print(as.character(unique(knownOrig[, i])))
      cat("\n")
    }
  }
  invisible(knownOrig)
}

