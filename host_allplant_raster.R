library(terra)
library(tidyverse)
library(tigris)
library(lubridate)
library(sf)

l1 <- list.files("Q:/My Drive/popblightdata/newareas/ctyhostcbs/combo/2010/", pattern = "tif$",
                 full.names = T)
r0 <- rast("Q:/My Drive/popblightdata/newareas/infections/2009/usa_2009_infections.tif")
r0 <- r0[[1]]
values(r0 > 0) <- 0

for(n in 2011:2020){
ln1 <- list.files(paste0("Q:/My Drive/popblightdata/newareas/ctyhostcbs/combo/",n,"/")
                  , pattern="tif$", full.names=T)
tem1 <- mosaic(rast(l1[1]), rast(l1[2]), fun="min")
for(i in 3:length(l1)){
  temr1 = rast(l1[i])
  tem1 = mosaic(tem1, temr1, fun = "min")
}

tem2 <- extend(tem1, r0)
txy1 <- xyFromCell(tem2, which(tem2[]==100))
r1 <- r0
r1[cellFromXY(r1, txy1)] <- 100
r2 = c(r1,r1)
for(i in 3:52){
  r2 = c(r2, r1)
}
writeRaster(r2, paste0("Q:/My Drive/temp/host_",n,".tif"))
}


l2 <- list.files("Q:/My Drive/relatedtoCDL/cdl1000/", pattern="tif$", full.names = T)
l2 <- l2[c(1:13, 15:27)]
# Set the <1 to NA

cdlcrrcmt3 <- matrix(c(0,0,NA,1,7,100,8,10,NA,11,15,100,16,21,NA,22,40,100,
                       41,41,NA,42,43,100, 44,44,100,45,54,100,55,55,100,56,62,100,
                       63,63,NA,64,73,100,74,74,NA,75,78,100,79,81,NA,82,84,0,
                       85,87,NA,88,89,0,93,93,0,94,111,NA,112,113,0,114,121,NA,
                       122,125,0,126,131,NA,132,132,0,131,141,NA,142,144,100,
                       145,152,NA,153,153,100,154,176,NA,177,177,100,178,190,NA,
                       191,191,100,192,195,NA,196,196,100,197,204,NA,205,
                       228,100,229,229,NA,230,251,100,252,255,NA,256,256,100), ncol = 3,
                     byrow = T)

cdlrcpt <- matrix(c(0,42,0,43,43,100,44,53,0,54,54,100,55,256,0), ncol = 3,
                  byrow = T)

indseq <- seq(3,20, by=2)

for(i in 1:length(indseq)){
  sr1 <- rast(l2[1])[[indseq[i]]]
  sr1[sr1 < 1] <- NA
  sr2 <- rast(l2[2])[[indseq[i]]]
  sr2[sr2 < 1] <- NA
  
  tem3 <- mosaic(sr1, sr2)
  
  for(j in 3:length(l2)){
    temr2 = rast(l2[j])[[indseq[i]]]
    temr2[temr2 < 1] <- NA
    tem3 = mosaic(tem3, temr2)
  }
  
  col1 <- c(tem3, tem3)
  
  for(k in 3:52){
    col1 = c(col1, tem3)
  }
  
  yseq <- 2011:2020
  
  pott <- classify(col1, cdlrcpt, right=NA)
  writeRaster(pott, paste0("Q:/My Drive/temp/potto_",yseq[i],".tif"), overwrite=T)
  
  allp <- classify(col1, cdlcrrcmt3, right=NA)
  writeRaster(pott, paste0("Q:/My Drive/temp/allplant_",yseq[i],".tif"), overwrite=T)
}

