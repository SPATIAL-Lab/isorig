calRaster <- function(known, isoscape, mask = NULL, interpMethod = 2, NA.value = NA, ignore.NA = T) {

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
  if (require("rgdal")) {
    print("rgdal is loaded correctly")
  } else {
    print("trying to install rgdal")
    install.packages("rgdal")
    if (require("rgdal")) {
      print("rgdal installed and loaded")
    } else {
      stop("could not install rgdal")
    }
  }
  if (require("sp")) {
    print("sp is loaded correctly")
  } else {
    print("trying to install sp")
    install.packages("sp")
    if (require("sp")) {
      print("sp installed and loaded")
    } else {
      stop("could not install sp")
    }
  }
  if (require("varhandle")) {
    print("varhandle is loaded correctly")
  } else {
    print("trying to install varhandle")
    install.packages("varhandle")
    if (require("varhandle")) {
      print("varhandle installed and loaded")
    } else {
      stop("could not install varhandle")
    }
  }
  if (require("ggplot2")) {
    print("ggplot2 is loaded correctly")
  } else {
    print("trying to install ggplot2")
    install.packages("ggplot2")
    if (require("ggplot2")) {
      print("ggplot2 installed and loaded")
    } else {
      stop("could not install ggplot2")
    }
  }

  if (class(known) != "data.frame") {
    stop("known should be a data.frame, see help page of rescale function")
  }

  if (!is.null(mask)) {
    if (class(mask) == "SpatialPolygonsDataFrame") {

      s <- SpatialPointsDataFrame(coords = known[,1:2],
                                                  data = known, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      o <- over(s, mask)

      known <- known[!is.na(o), ]
    } else {
      stop("mask should be a SpatialPolygonsDataFrame")
    }
  }

  nSample <- nrow(known)
  tissue.iso <- vector("numeric", length = nSample)
  isoscape.iso <- vector("numeric", length = nSample)
  null.iso <- NULL

  if (class(isoscape) == "RasterLayer") {
    prediction <- isoscape
  } else {
    stop("isoscape should be a RasterLayer")
  }
  ncells <- ncell(prediction)
  temp <- apply(known,2,class)
  if (any(temp!="numeric")){
    stop("known data should be all numeric. Please check that except header, known table should only contains numbers.")
  }

  tissue.iso <- as.numeric(known[, 3])
  if (interpMethod == 1) {
    isoscape.iso <- raster::extract(prediction, known[, 1:2], method = "simple")
  } else if (interpMethod == 2) {
    isoscape.iso <- raster::extract(prediction, known[, 1:2], method = "bilinear")
  } else {
    stop("interpMethod should be either 1 or 2")
  }
  if (!is.na(NA.value)){
    values(isoscape)[values(isoscape)==NA.value] <- NA
  }
  if (any(is.na(isoscape.iso))){
    na <- which(is.na(isoscape.iso))
    cat("\n\n----------------------------------------------------------------\n")
    cat("Warning: NO data are found at following locations:\n")
    print(known[na,1:2])
    if (!ignore.NA){
      stop ("Delete these data in known origin data or use a different isoscape that has values at these locations")
    }
    isoscape.iso <- isoscape.iso[-na]
    tissue.iso <- tissue.iso[-na]
  }

  lmResult <- lm(tissue.iso ~ isoscape.iso)
  cat("\n\n---------------------------------------------------------------------------------\n")
  cat("rescale function uses linear regression model, the summary of this model is:\n")
  cat("---------------------------------------------------------------------------------\n")
  print(summary(lmResult))
  intercept <- as.numeric(coef(lmResult)[1])
  slope <- as.numeric(coef(lmResult)[2])
  temp <- getValues(isoscape)
  temp1 <- temp * slope + intercept
  isoscape.rescale <- setValues(isoscape,temp1)
  x = isoscape.iso
  y = tissue.iso
  xy = data.frame(x,y)
  equation <- function(x) {
    lm_coef <- list(a = round(coef(x)[1], digits = 2),
                    b = round(coef(x)[2], digits = 2),
                    r2 = round(summary(x)$r.squared, digits = 2));
    lm_eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,lm_coef)
    as.character(as.expression(lm_eq));
  }

  fit <- lm(y~x, data = xy)
  xmin = max(x)-0.4*(diff(range(x)))
  xmax = max(x)
  ymin = min(y)-0.25*diff(range(y))
  ymax = min(y)-0.1*diff(range(y))

  p11 <- ggplot(xy, aes(x=isoscape.iso, y=tissue.iso)) + geom_point(shape=1) + geom_smooth(method=lm, se=T) +
    ggtitle("Rescale regression model")+theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(name = "Precipitation isotope") +
    scale_y_continuous(name = "Tissue isotope") +
    annotate("rect", xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill="white", colour="red") +
    annotate("text", x = (xmax+xmin)/2, y = (ymax+ymin)/2,label = equation(fit), parse = TRUE)
  print(p11)

  sd <- isoscape
  values(sd) <- summary(lmResult)$sigma
  values(sd)[is.na(values(isoscape))] <- NA
  isoscape.rescale <- stack(isoscape.rescale, sd)
  names(isoscape.rescale) <- c("mean", "sd")
  if (!is.null(mask)) {
  isoscape.rescale <- crop(isoscape.rescale,mask)
  }
  print(spplot(isoscape.rescale$mean, scales = list(draw = TRUE), main = "rescale mean"))
  print(spplot(isoscape.rescale$sd,scales = list(draw = TRUE), main="rescale sd"))

  dir.create("output")
  pdf("./output/rescale.result.pdf",width = 8,height = 4)
  print(p11)
  print(spplot(isoscape.rescale$mean, scales = list(draw = TRUE), main = "rescale mean"))
  print(spplot(isoscape.rescale$sd,scales = list(draw = TRUE), main="rescale sd"))
  dev.off()

  names(xy) = c("isoscape.iso","tissue.iso")
  result = list ("isoscape.rescale" = isoscape.rescale, "lm.data" = xy, "lm.model" = lmResult)
  return(result)
}


