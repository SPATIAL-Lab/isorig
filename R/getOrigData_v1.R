getOrigData <- function() {

  data(knownOrig)
  print("Known origin data has following catagories:")
  colnames(knownOrig)
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

