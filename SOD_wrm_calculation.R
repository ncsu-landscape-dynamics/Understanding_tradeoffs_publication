# Early in the project, we called the WRM "DSV". We decided WRM would be the acronym for the publication.
# So everywhere it says "DSV", that's the WRM.

library(terra)
library(tigris)
library(tidyverse)

# Function to assign DSV score based on precipiation and mean T.
dsv <- function(precip, tmean) {
  prec_val = ifel(precip >= 2, 1, 0)
  
  t_val = ifel(tmean > 14.5 & tmean < 24.5, 4, 
        ifel(tmean > 11.75 & tmean <= 14.5 | tmean >=24.5 & tmean < 25.75, 3,
             ifel(tmean <= 11.75 & tmean > 9 | tmean >= 25.75 & tmean < 27.0, 2,
                 ifel(tmean <= 9 & tmean > 5.75 | tmean <= 28.25 & tmean >= 27 , 1, 0))))
  
  dsv_val = prec_val * t_val
  
  return(dsv_val)
}

# Function for use in for loop
dsvf <- function(i){
  mx = rast(txl[i])
  mn = rast(tnl[i])
  pr = rast(prl[i])
  mx = crop(mx, co1, mask=T)
  mn = crop(mn, co1, mask=T)
  precip = crop(pr, co1, mask=T)
  
  tmean = (mx+mn)/2
  print(tmean)
  
  codsv = dsv(precip, tmean)
  return(codsv)
}

# Geographical focus
co1 <- counties(state = c("OR"))
co1 <- co1[co1$NAME %in% c("Coos","Curry","Douglas","Josephine"),]

# Weather data for DSV
txl <- list.files("Q:/Shared drives/Data/Original/Daymet/tmax", full.names = T)
tnl <- list.files("Q:/Shared drives/Data/Original/Daymet/tmin", full.names = T)
prl <- list.files("Q:/Shared drives/Data/Original/Daymet/precip/", full.names = T)

# The lists above grab all years. Select correct years 
txl <- txl[...:...]
tnl <- tnl[...:...]
prl <- prl[...:...]

# Set the focus to matching projection
r0 <- rast(txl[1])
co1 <- vect(co1)
co1 <- project(co1, y=crs(r0))

tdsv <- list()

for(i in 1:length(txl)){
  tdsv[i] = dsvf(i)
}

# Project the DSVs to EPSG 3857 to match all other files.
dct1 = list()

for(i in 1:8){
  #indq = seq(1, 4380, by=365)
  #temr = r1[[indq[i]:(indq[i]+364)]]
  dct1[i] = project(rast(tdsv[i]), crs(r2))
}

# The resolution was wrong after projecting. 
# Used one of the weather rasters as a template for resampling in order
# to have right resolution.
res1 = list()
for(i in 1:20){
  res1[i] = resample(rast(dct1[i]), r2, method="near")
}

# Get the list of files for revising further.
l1 <- list.files("Q:/My Drive/temp/",pattern="soddsvdaily.*tif$", full.names = T)

# Fix some values that resampled as decimal
rm1 <- matrix(c(0,0.5,0,0.5,1,1,1,1.5,1,1.5,2,2,2,2.5,2,2.5,3,3,3,3.5,3,3.5,4,4),
              ncol=3, byrow = T)

yseq <- 2001:2020

for(i in 1:20){
  temr = rast(res1[i])
  temr2 = classify(temr, rm1, right=T)
  writeRaster(temr2, paste0("D:/data/ncsu/googledrop/newareas/SOD/soddsvdaily_", yseq[i],".tif"), overwrite=T)
}

# Cumulative sum of the daily values in DSV.
cd1 <- list()

l1 <- list.files("D:/data/ncsu/googledrop/newareas/SOD/", pattern="soddsvdaily.*tif$", full.names = T)

for(i in 1:20){
  cd1[i] = app(rast(l1[i]), cumsum)
}

# Cumulative counts of DSV/18
cd2 <- list()

for(i in 1:12){
  cd2[i] = app(rast(cd1[i]), function(x) floor(x / 18))
}

yseq <- 2009:2020

l1 <- list.files("Q:/My Drive/temp/", pattern = "cumu", full.names = T)

for(i in 1:20){
  l1 <- list.files(paste0("Q:/My Drive/popblightdata/newareas/DSV/",yseq[i],"/daily/"), full.names = T)
  l2 <- l1[c(1, 112, 223, 300, 311, 322, 333, 344, 355, 2, 13, 24, 35, 46, 57, 68, 79, 90, 101, 113, 124, 135, 146, 157, 168, 179, 190, 201, 212, 224, 235, 246, 257, 268, 279, 290, 297:299, 301:310, 312:321, 323:332, 334:343, 345:354, 356:365, 3:12, 14:23, 25:34, 36:45, 47:56, 58:67, 69:78, 80:89, 91:100, 102:111, 114:123, 125:134, 136:145, 147:156, 158:167, 169:178, 180:189, 191:200, 202:211, 213:222, 225:234, 236:245, 247:256, 258:267, 269:278, 280:289, 291:296)]
  cd1 <- app(rast(l1), cumsum)
  cd2 <- app(cd1, function(x) floor(x/18))
  for(j in 1:365){
    jpeg(paste0("Q:/My Drive/popblightdata/newareas/DSV/",yseq[i],"/cumulthreshold/cuthresh",yseq[i],"_",j,".jpg"))
    plot(cd2[[j]])
    dev.off()
  }
}
