library(terra)
library(tidyverse)

# Read in rast
r1 <- rast("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs21/locs10x2/upwdrev/rasts/pops_median_Year_2021.tif")

# Frequency
freq(r1)

# Make reduction repeatable
set.seed(42)

# Which cells to consider
cwv1 <- which(values(r1)>=1)

# Function
setcll0 <- function(x) {
  r_tem = r1
  # Pick X% of cells
  smlset = sample(cwv1, floor(length(cwv1)*x))
  # Set selected cells to same value as background
  r_tem[smlset] = 0
  pct = 100 - (x*100)
  writeRaster(r_tem, paste0("Z:/Late_blight/Manuscript_1_Data/simulation/LB/outputs21/locs10x2/upwdrev/rasts/med_loc10_",
              pct, ".tif"))
}

pseq <- seq(0.5, 0.7, by = 0.1)

raslis <- lapply(pseq, setcll0)

