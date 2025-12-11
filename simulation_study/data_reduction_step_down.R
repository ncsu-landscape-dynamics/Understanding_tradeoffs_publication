library(terra)
library(tidyverse)

# Read in raster
# 2021
#r1 <- rast("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs21/locs10x2/upwdrev/rasts/pops_mean_Year_2021.tif")
#r1b <- rast("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs21/locs10x2/upwdrev/rasts/mean2022redczero.tif")

# 2022
#r1 <- rast("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs22/locs10x2/upwdrev/rasts/pops_mean_Year_2022.tif")
#r1b <- rast("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs22/locs10x2/upwdrev/rasts/mean2022redczero.tif")

# 2023
#r1 <- rast("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs23new/rasts/pops_mean_Year_2023.tif")
#r1b <- rast("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs22/locs10x2/upwdrev/rasts/mean2022redczero.tif")


# Make reduction repeatable
set.seed(42)

# Ensure whole numbers
r1 <- round(r1)

# Which cells to consider
cwv1 <- which(values(r1)>=1)

# Function
setcll0 <- function(x, inst) {
  r_tem = r1
  # Pick X% of cells
  smlset = sample(cwv1, round(length(cwv1)*x))
  print(paste0("x: ", x))
  print(length(smlset))
  # Set selected cells to same value as background
  r_tem[smlset] = 0
  
  pct = 100 - (x*100)
  #print(paste0("pct: ", pct))
  #return(r_tem)
  writeRaster(r_tem, paste0("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs23new/rasts/r",
    pct,"/mean_loc10_",pct, "-", inst, "round.tif"), overwrite=T)
  #writeRaster(r_tem, paste0("data/temp/new21/r",pct,"/mean_loc10_",
  #    pct, "-", inst, ".tif"), overwrite=T)
}

pseq <- seq(0.1, 0.8, by = 0.1)

lapply(pseq, \(x) lapply(1:10, \(y) setcll0(x, y)))

# Change the number of zeroes to match number of infected locations
# Get a list of the reduced rasters
l2 <- list.files("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs23new/rasts/",
     full.names = T, pattern = "round", recursive = T)

host1 <- rast("Z:/Late_blight/Manuscript_1_Data/simulation/LB/inputs/host/host_no_na_rv_2023.tif")

setzf <- function(x) {
  r1 = rast(x)
  nm1 = gsub("\\.tif", "", str_sub(x, 82)) # 96)) for '21 and '22
  pct1 = str_sub(nm1, 1,2)
  print(paste0("name ",nm1))
  print(paste0("pct ", pct1))
  # Get counts of locations
  td1 = as.data.frame(freq(r1))
  td1 = td1[-1,]
  print(td1)
  n_locs = sum(td1$count)
  print(n_locs)
  # New raster for assigning NA and reassigning 0
  r2 = round(r1)
  r2 = crop(r2, host1, mask = T)  
  nonzer = which(r2[] > 0)
  zerloc = which(r2[] == 0)
  print(length(zerloc))
  r2[r2 == 0] = NA
  newzero = sample(zerloc, n_locs)
  print(length(newzero))
  
  r2[newzero] = 0
  #print(paste0("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs23new/rasts/r",
  #             pct1,"/mean10loc_",nm1,"_redc.tif"))
  writeRaster(r2, paste0("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs23new/rasts/r",
                         pct1,"/mean10loc_",nm1,"_redc.tif"), overwrite =T)
}

lapply(l2, \(x) setzf(x))

# Create a dataframe or file of the counts for each reduced iteration.
l3 <- gtools::mixedsort(list.files("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs23new/rasts/",
                 full.names = T, pattern = "round_redc", recursive = T))
                              

fret1 <- function(x) {
  rtm1 = rast(x)
  df1 = as.data.frame(freq(rtm1))
  df1$tot = df1$value*df1$count
  ttot = colSums(df1)[4]
  return(ttot)
}

raslis <- lapply(l3, fret1)
df1 <- data.frame(pct = rep(seq(10,90,by=10), each = 10), index=rep(seq(10,1,by=-1), 9), n_infx = unlist(raslis), n_loc = unlist(lapply(l3, \(x) length(which(values(rast(x) > 0))))))

write.csv(df1, "Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs23new/rasts/infxlocsmeanrndrasts2023.csv", row.names=F)



