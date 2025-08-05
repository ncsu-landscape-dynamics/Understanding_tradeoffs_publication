library(terra)

cdll1 <- list.files("Z:/Data/Raster/USA/landcover/", recursive=T, full.names=T)

cdll1 <- cdll1[grepl("tif$|img$", cdll1)]
cdll1 <- cdll1[14:16]
cdlr1 <- rast(cdll1)

# Core region
cor1 <- vect("Z:/Late_blight/helpful_shapes_rasters/lb_centralregion.geojson")

templ1k <- rast("Z:/Late_blight/hosts/no_na_no_mask_revsd/redchost_no_na_rv_2020.tif")
cor2 <- project(cor1, y = crs(templ1k))
cor2 <- crop(templ1k, y = cor2, mask = T)

cor3 <- project(cor1, y = crs(cdlr1))


# Change the filter to remove big crops
cdlcrrcmt3 <- matrix(c(0,0,NA,1,24,0,25,39,100,40,40,0,41,61,100,62,204,0,205,
                       227,100,228,228,0,229,250,100,251,
                       254,0,255,255,100), ncol = 3, byrow = T)


# Change the filter to remove big crops
cdlcrrcmt4 <- matrix(c(0,0,NA,1,40,0,41,61,100,62,204,0,205,
                       227,100,228,228,0,229,250,100,251,
                       254,0,255,255,100), ncol = 3, byrow = T)

yseq <- 2021:2023

# This is more succinct
newrv <- function(x){
  r1 = rast(cdll1[x])
  print(x)
  # crop to area
  cr1 = crop(r1, cor3, mask = T)
  # Mask has to be done separately
  #  ma1 = mask(cr1, us3)
  print("crop")
  cl1 = classify(cr1, cdlcrrcmt4, right=NA, overwrite=T)
  cl1[cl1 == 100] = 0.1
  ag1 = terra::aggregate(cl1, 33, fun="sum")
  ag2 = clamp(ag1, lower=0, upper=100)
  pr1 = project(ag2, crs(cor2), overwrite=T)
  re1 = resample(pr1, cor2, method="near", overwrite=T)
  print("resample")
  # Hosts
  re1[is.na(re1)] = 0
  re1 = round(re1, 0)
  # Total pops
  tp1 = re1
  tp1[tp1 > 0] = 100
  writeRaster(re1, paste0("Z:/Late_blight/hosts/no_na_no_mask_revsd/host_no_na_rv_B_",yseq[x],".tif"))
  writeRaster(tp1, paste0("Z:/Late_blight/all_populations/no_na_no_mask_revsd/tot_pop_no_na_rv_B_",yseq[x],".tif"))
}

for(i in 1:length(cdll1)){
  newrv(i)
}

# This function is relevant
crreag <- function(x){
  r1 = rast(cdll1[x])
  print("start crop")
  print(Sys.time())
  # Do mask separately
  cr1 = crop(x=r1, y=us3)
  ma1 = mask(cr1, y=us3)
  print("start classify")
  print(Sys.time())
  rec1 = classify(ma1, cdlcrrcmt3, right=NA)
  # Continue for total populations
  print("start aggregate")
  print(Sys.time())
  agg1 = aggregate(rec1, 33, fun="max")
  
  print("start project")
  print(Sys.time())
  pr1 = project(agg1, crs(r0))
  print("start resample")
  print(Sys.time())
  res1 = resample(pr1, r0, method="near")
  res1[res1 > 0] = 100
  #mk1 = mask(res1, hz1, maskvalue=0, updatevalue=0)
  #mk1[is.na(mk1)] = 0
  res1[is.na(res1)] = 0
  writeRaster(res1, paste0("Z:/Late_blight/all_populations/tot_pop_",
                           yseq[x],"_nona_nomask.tif"), overwrite=T)
  #return(res1)
  ## For hosts
  host1 = rec1
  host1[host1 == 100] = 0.1
  print("start host aggregate")
  print(Sys.time())
  agg2 = aggregate(host1, 33, fun="sum")
  agg2 = clamp(agg2, lower=0, upper=100)
  print("start host project")
  print(Sys.time())
  pr2 = project(agg2, crs(r0))
  print("start host resample")
  print(Sys.time())
  res2 = resample(pr2, r0, method="near")
  #mk2 = mask(res2, hz1, maskvalue=0, updatevalue=0)
  #mk2[is.na(mk2)] = 0
  res2[is.na(res2)] = 0
  writeRaster(res2, paste0("Z:/Late_blight/hosts/host_",
                           yseq[x],"_nona_nomask.tif"), overwrite=T)
  
}
