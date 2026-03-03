library(terra)
library(stringr)

cor1 <- vect("Z:/Late_blight/helpful_shapes_rasters/lb_centralregion.geojson")
templ1k <- rast("Z:/Late_blight/hosts/no_na_no_mask_revsd/redchost_no_na_rv_2020.tif")
cor2 <- project(cor1, y = crs(templ1k))
cor2 <- crop(templ1k, cor2, mask = T)

l1 <- list.files("Z:/Late_blight/Weather_data/", pattern = "core", full.names = T)
l1 <- l1[str_detect(l1, "b.tif")]
l2 <- gsub("b.tif",".tif", l1)

l1 <- list.files("Z:/Data/Original/Daymet/VaporPressure", full.names = T,
                 pattern = "2021|2022|2023")
l2 <- paste0("Z:/Late_blight/Weather_data/core_vp_",2021:2023,".tif")

pcwrf <- function(x){
  r1 = rast(l1[x])
  print(x)
  print(crs(r1))
  nm1 = l2[x]
  r2 = project(r1, cor2)
  print(crs(r2))
  r3 = crop(r2, cor2, mask =T)
  r4 = resample(r3, cor2, "near")
  r4 = round(r4, 2)
  writeRaster(r4, nm1)
}

for (i in 1:length(l1)){
  pcwrf(i)
}

for (i in 1:3){
  r1 = rast(l1[i])
  r2 = 24 - r1
  writeRaster(r2, l2[i])
}

dl1 <- list.files('Z:/Late_blight/Weather_data/', full.names = T, pattern = "dayl[0-9]")
dl2 <- gsub("yl","ylhr",dl1)

pcrrwf <- function(x){
  r1 = rast(dl1[x])
  print(x)
  print(crs(r1))
  nm1 = dl2[x]
  r2 = project(r1, cor2)
  print(crs(r2))
  r3 = crop(r2, cor2, mask =T)
  r4 = resample(r3, cor2, "near")
  r4 = round(r4, 2)
  # sunrise
  r5 = 12 - 0.5*r4
  nm2 = gsub("daylhr","srs",nm1)
  # sunset
  r6 = 12 + 0.5*r4
  nm3 = gsub("daylhr","sst",nm1)
  writeRaster(r4, nm1, overwrite=T)
  writeRaster(r5, nm2)
  writeRaster(r6, nm3)
}

for (i in 1:3){
  pcrrwf(i)
}
