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
