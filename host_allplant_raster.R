library(terra)
library(tigris)
library(sf)

cdll1 <- list.files("Z:/Late_blight/temp/cdls", recursive=T, full.names=T)

cdll1 <- cdll1[grepl("tif$|img$", cdll1)]


r0 <- rast("Z:/Late_blight/Weather_data/weather_us_2009.tif")
r0 <- r0[[1]]

r0[!is.na(r0)] = 100
r0[r0 == 100] = 0
rmask <- r0
rmask[is.na(rmask)] <- 1
rmask[rmask != 1] <- NA

crmask <- rast("Z:/Late_blight/helpful_shapes_rasters/cropmaskforCDL_albers.tif")

cdlcrrcmt3 <- matrix(c(0,0,NA,1,6,100,7,9,0,10,14,100,15,20,0,21,39,100,
                       40,40,0,41,61,100,62,62,0,63,72,0,73,73,0,74,204,0,
                       205,227,100,228,228,0,229,250,100,251,254,0,255,255,100), ncol = 3,
                     byrow = T)

# Change the filter to remove big crops
cdlcrrcmt3 <- matrix(c(0,0,NA,1,24,0,25,39,100,40,40,0,41,61,100,62,204,0,205,
                       227,100,228,228,0,229,250,100,251,
                       254,0,255,255,100), ncol = 3, byrow = T)




#us2 <- st_read("data/gis_data/all_states/lower48_bound/state_nrcs_aea.shp")
#us2 <- us2[!us2$STPO %in% c("AZ", "HI", "VI", "MP", "GU", "AS", "PR", 
#                           "IA", "KS","MN","ND","ID","NE","OR","SD",
#                           "AK","LA","MT","CA","CO","NB","UT","WA","AR",
#                           "MO","NM","OK","WY","NV","MS","TX"),]
#us2 <- st_union(us2)
#us2 <- st_transform(us2, crs="epsg:3857")
#us2 <- vect(us2)
#us3 <- project(us2, y=crs(rast(cdll1[2])))
us3 <- vect("Z:/Late_blight/helpful_shapes_rasters/crmask_for_host_processing.gpkg")
us3 <- project(us3, y=crs(rast(cdll1[1])))

yseq <- 2009:2020

# This is more succinct
newrv <- function(x){
  r1 = rast(cdll1[x])
  print(x)
  # crop to area
  cr1 = crop(r1, us3)
  # Mask has to be done separately
  ma1 = mask(cr1, us3)
  print("crop")
  cl1 = classify(ma1, cdlcrrcmt3, right=NA, overwrite=T)
  cl1[cl1 == 100] = 0.1
  ag1 = terra::aggregate(cl1, 33, fun="sum")
  ag2 = clamp(ag1, lower=0, upper=100)
  pr1 = project(ag2, crs(r0), overwrite=T)
  re1 = resample(pr1, r0, method="near", overwrite=T)
  print("resample")
  # Hosts
  re1[is.na(re1)] = 0
  re1 = round(re1, 0)
  # Total pops
  tp1 = re1
  tp1[tp1 > 0] = 100
  writeRaster(re1, paste0("Z:/Late_blight/hosts/no_na_no_mask_revsd/host_no_na_rv_",yseq[x],".tif"))
  writeRaster(tp1, paste0("Z:/Late_blight/all_populations/no_na_no_mask_revsd/tot_pop_no_na_rv_",yseq[x],".tif"))
}

for(i in 2:length(cdll1)){
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
