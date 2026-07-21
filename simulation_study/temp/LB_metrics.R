library(terra)
library(sf)

# Reference raster
refras <- rast("Z:/Late_blight/helpful_shapes_rasters/econus_mask_3857.tif")

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

# List of positives from "no filter"
nfdfl <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/lateblight", 
                                      pattern = "lbnof", full.names = T))

# List of positives from "filter then reduce"
ftdfl <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/lateblight", 
                                      pattern = "lbftr", full.names = T))

# List of pseudoabsences
psalis <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/lateblight/",
                                       pattern = "psa", full.names = T))

# create pos raster
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
  # Switch between these for either no filter or filter then reduce. 
  write.csv(outdf, paste0("Z:/Late_blight/SDM/lateblight/nofilter_mets_stable_",yrlis[x],".csv"), row.names = F)
  #write.csv(outdf, paste0("Z:/Late_blight/SDM/lateblight/filtertd_metrics_stable_",yrlis[x],".csv"), row.names = F)
}

# No filter predictions
predlist <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/lateblight", full.names = T, 
                                         pattern = "relcls_unf_"))

# Filter then reduce predictions
predlist <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/lateblight", full.names = T, 
                                         pattern = "relcls_ftd_"))

# Years
yrlis <- seq(20, 100, by = 10)

#lapply(1:9, \(x) metric_calcf(nfdfl[x], predlist[x], x))
#lapply(1:9, \(x) metric_calcf(ftdfl[x], predlist[x], x))
#lapply(1:9, \(x) metric_calcf(valras, predlist[x], x))
#lapply(1:9, \(x) metric_calcf(valras, predlist[x], x))
# Fix to 100% -- [9]
lapply(1:9, \(x) metric_calcf(nfdfl[9], psalis[9], predlist[x], x))


# List of stats by %, No filter 
metlis1 <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/lateblight/",
                                        pattern="nofilter.*mets.*stab", full.names = T))
# Make dataframe from list of stats
metdf2 <- do.call(rbind,lapply(rev(metlis1), read.csv))


# Filter then reduce
metlis1 <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/lateblight/", 
                                        pattern = "filtert.*metr.*stab", full.names = T))

metlis2 <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/lateblight/",
                                        pattern = "recldf_ftd", full.names = T))

