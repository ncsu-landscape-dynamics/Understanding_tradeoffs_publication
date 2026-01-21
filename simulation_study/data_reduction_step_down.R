library(terra)
library(tidyverse)

# Set path as needed
path21 <- "Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs21/locs10x2/upwdrev/rasts/"
path22 <- "Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs22/locs10x2/upwdrev/rasts/"
path23 <- "Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs23new/rasts/"

source("code/blight_related/reduce_data_function.R")

# Read in raster
# Raster to reduce
red_raster <- rast("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs22/locs10x2/upwdrev/rasts/pops_mean_Year_2022.tif")
# Raster for original locations
#orig_raster <- rast("Z:/Late_blight/Manuscript_1_Data/simulation/LB/inputs/infection/infection_2023_add10.tif")
orig_raster <- rast("Z:/Late_blight/Manuscript_1_Data/simulation/LB/inputs/infection/infection_2022_rev_10locs.tif")

# Make reduction repeatable
set.seed(42)

# Ensure whole numbers
red_raster <- round(red_raster)

# Original infection cells
cells2avoid <- which(values(orig_raster)>=1)

# Cells in raster that will be reduced
endingcells <- which(values(red_raster)>=1)

# Create set of cells that have difference of cells of raster to reduce
# and the set of original infections.
cells2samplefrom <- endingcells[!endingcells %in% cells2avoid]

# Optional
#outputpath <- "Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs21/locs10x2/upwdrev/rasts/r"

# Set up the percentages to remove
pseq <- seq(0.1, 0.9, by = 0.1)

# Use pseq as amount to reduce with 1 - 10 as index for iteration with cells2samplefrom
lapply(pseq, \(x) lapply(1:10, \(y) location_reduction(cells2samplefrom, x, y)))

##
##
## The rest is probably only necessary for late blight

# Get a list of the reduced rasters
l2 <- list.files("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs23new/rasts/",
     full.names = T, pattern = "round", recursive = T)

# Reduce N of zeros.
setzf <- function(x) {
  r1 = rast(x)
  nm1 = gsub("\\.tif", "", str_sub(x, 82)) # 96)) for '21 and '22
  pct1 = str_sub(nm1, 1,2)
  print(paste0("name ",nm1))
  print(paste0("pct ", pct1))

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

write.csv(df1, ..., row.names=F)



