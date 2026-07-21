library(terra)
library(sf)

# Reference raster
refras <- rast("Z:/Late_blight/SOD_OR/sod_PoPSpred_resampledtomatchSDMWRM26710.tif")

# Observed
obspt <- vect("Z:/Late_blight/helpful_shapes_rasters/blightptsmovedforSDM.gpkg")

obspt$id <- paste0(obspt$x,"_",obspt$y)
obspt <- obspt[!duplicated(obspt$id),]
obspt <- obspt[obspt$Year %in% 2008:2021,]
obspt <- project(obspt, y =crs(refras))
obr1 <- rasterize(obspt, refras)
refras[refras > 0] = 0
valras <- obr1
valras[is.na(valras)] = 0
valras <- crop(valras, refras, mask = T)

nofil <- c(796,716,637,557,478,398,318,239,159)
filtred <- c(381, 343, 305, 267, 229, 190, 152, 114, 76)

nfdfl <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/SOD", 
                                      pattern = "sodnof", full.names = T))

ftdfl <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/SOD", 
                                      pattern = "sodftr", full.names = T))

# create positive raster
newvalraf <- function(x){
  obspt <- read.csv(x)
  obspt <- vect(obspt, geom=c("x","y"), crs = "epsg:4326")
  obspt <- project(obspt, y =crs(refras))
  obr1 <- rasterize(obspt, refras)
  refras[refras > 0] = 0
  valras <- obr1
  valras[is.na(valras)] = 0
  valras <- crop(valras, refras, mask = T)
  return(valras)
}

# create psa raster
newngraf <- function(x){
  obspt <- read.csv(x)
  obspt <- vect(obspt, geom=c("x","y"), crs = "epsg:4326")
  obspt <- project(obspt, y =crs(refras))
  obr1 <- rasterize(obspt, refras)
  refras[refras > 0] = 0
  valras <- obr1
  valras[is.na(valras)] = 0
  valras <- crop(valras, refras, mask = T)
  valras[valras > 0] <- 10
  return(valras)
}

# Compare observation and prediction and calculate metrics
metric_calcf <- function(obsdf, ngrdf, pred1, x) {
  print(pred1)
  pred1 = rast(pred1)
  #obr1 = newvalraf(obsdf)
  # Raster of observation/truth
  truras = newvalraf(obsdf)
  # Raster of pseudoabsences
  psaras = newngraf(ngrdf)
  # Combine obs/pseudo
  obs1 <- truras + psaras
  obs1[obs1 < 1] <- NA
  obs1[obs1 == 10] <- 0
  obr1 = obs1
  print(freq(obr1))
  # Check for no values
  if (any(values(obr1) > 0) == FALSE) 
    print("Observations don't exist!")
  
  if (any(values(pred1) > 0) == FALSE)
    print("Predictions don't exist!")
  
  # True positives
  trpo1 <- sum(values(obr1) > 0 & values(pred1) > 0, na.rm = T)
  
  # False positives
  fapo1 <- sum(values(obr1) == 0 & values(pred1) >0, na.rm = T)
  
  # True negatives
  trne1 <- sum(values(obr1) == 0 & values(pred1) == 0, na.rm = T)
  
  # False negatives
  fane1 <- sum(values(obr1) > 0 & values(pred1) == 0, na.rm = T)
  
  #  Calculations for TSS
  sensit <- trpo1 / (trpo1 + fane1)
  
  specif<- trne1 / (trne1 + fapo1)
  
  mtss <- (sensit + specif ) - 1
  
  preci1 <- trpo1 / (trpo1 + fapo1)
  
  outdf <- data.frame(tru_pos = trpo1, fal_pos = fapo1, tru_neg = trne1, 
                      fal_neg = fane1, sensitivity = sensit, specificity = specif,
                      tru_skil = mtss, precision = preci1)
  #return(outdf)
  write.csv(outdf, paste0("Z:/Late_blight/SDM/SOD/nofilter_metrs_stable",yrlis[x],".csv"), row.names = F)
  #write.csv(outdf, paste0("Z:/Late_blight/SDM/SOD/filtertd_metrics_",yrlis[x],".csv"), row.names = F)
}

predlist <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/SOD", full.names = T, 
                                         pattern = "relcls_unf_"))
predlist <- predlist[10:19]

predlist <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/SOD", full.names = T, 
                                         pattern = "relcls_ftd_"))

predlist <- predlist[10:19]

yrlis <- seq(10, 100, by = 10)

#lapply(1:10, \(x) metric_calcf(nfdfl[x], predlist[x], x))
#lapply(1:10, \(x) metric_calcf(ftdfl[x], predlist[x], x))
# Fix the obs and psa to 100%
lapply(1:10, \(x) metric_calcf(nfdfl[10], psalis[10], predlist[x], x))

# No filter
metlis1 <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/SOD/",
                                        pattern="nofilter.*metrs.*stabl", full.names = T))

metlis2 <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/SOD/",
                                        pattern = "lohict_unf", full.names = T))


# Filter then reduce
metlis1 <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/SOD/", 
                                        pattern = "filtert.*metr", full.names = T))

metlis2 <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/SOD/",
                                        pattern = "recldf_ftd", full.names = T))

