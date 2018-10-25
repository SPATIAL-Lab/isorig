calRaster <- function(known, isoscape, mask = NULL, interpMethod = 2, NA.value = NA, ignore.NA = T) {

  if (class(known) != "SpatialPointsDataFrame") {
    stop("known should be a SpatialPointsDataFrame, see help page of rescale function")
  } else {
    s <- known
  }
  if(is.na(proj4string(known))){
    stop("known must have coord. ref. of +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  }

  if (!is.null(mask)) {
    if (class(mask) == "SpatialPolygonsDataFrame") {
      if (is.na(proj4string(mask))){
        stop("mask must have coord. ref.")
      } else {
        mask <- spTransform(mask, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      }

      #s <- SpatialPointsDataFrame(coords = known[,1:2],
      #                                            data = known, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      o <- over(s, mask)
      known <- as.data.frame(known)
      known <- known[!is.na(o), ]
    } else {
      stop("mask should be a SpatialPolygonsDataFrame")
    }
  } else {
    known <- as.data.frame(known)
  }

  nSample <- nrow(known)
  tissue.iso <- vector("numeric", length = nSample)
  isoscape.iso <- vector("numeric", length = nSample)
  null.iso <- NULL

  if (class(isoscape) == "RasterLayer") {
    if (is.na(proj4string(isoscape))){
      stop("isoscape must have coord. ref.")
    } else {
      isoscape <- projectRaster(isoscape, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    }
    prediction <- isoscape
  } else {
    stop("isoscape should be a RasterLayer")
  }
  ncells <- ncell(prediction)
  temp <- apply(known,2,class)
  if (any(temp!="numeric")){
    stop("known data should be all numeric. Please check that except header, known table should only contains numbers.")
  }

  tissue.iso <- as.numeric(known[, 1])
  if (interpMethod == 1) {
    isoscape.iso <- raster::extract(prediction, known[, 2:3], method = "simple")
  } else if (interpMethod == 2) {
    isoscape.iso <- raster::extract(prediction, known[, 2:3], method = "bilinear")
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
    lm_coef <- list(a = as.numeric(round(coef(x)[1], digits = 2)),
                    b = as.numeric(round(coef(x)[2], digits = 2)),
                    r2 = round(summary(x)$r.squared, digits = 2));
    lm_eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,lm_coef)
    as.character(as.expression(lm_eq))
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

  crs(isoscape.rescale) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  names(xy) = c("isoscape.iso","tissue.iso")
  result = list ("isoscape.rescale" = isoscape.rescale, "lm.data" = xy, "lm.model" = lmResult)
  return(result)
}


