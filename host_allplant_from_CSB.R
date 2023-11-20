library(tigris)
library(sf)
library(tidyverse)
library(terra)

# Read in the CSB
csb0916 <- read_sf(<path to crop sequence boundary database>)

# Get the study area geography needed for late blight
st1 <- states()
stl1 <- c("CT", "DE", "FL", "GA", "IL", "IN", "KY", "MA", "MD", "ME", "MI",
                  "NC", "NH", "NJ", "NY", "OH", "PA", "RI", "SC",
                  "TN", "VA", "VT", "WI", "WV", "AL")
st1 <- st1[st1$STATEFP %in% stl1,]

# Subset to the host and study area
csbhost09 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R09 %in% c(43,54),c()]
csbhost10 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R10 %in% c(43,54),c()]
csbhost11 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R11 %in% c(43,54),c()]
csbhost12 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R12 %in% c(43,54),c()]
csbhost13 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R13 %in% c(43,54),c()]
csbhost14 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R14 %in% c(43,54),c()]
csbhost15 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R15 %in% c(43,54),c()]
csbhost16 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R16 %in% c(43,54),c()]
csbhost17 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R17 %in% c(43,54),c()]
csbhost18 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R18 %in% c(43,54),c()]
csbhost19 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R19 %in% c(43,54),c()]
csbhost20 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R20 %in% c(43,54),c()]

# Subset to the non-host plants and study area
csballp09 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R09 %in% c(1:6,
             10:14,21:39,41:42,44,53,55:61,63,141:143,176,190,195,204:227,229:250,255),]
csballp10 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R10 %in% c(1:6,
             10:14,21:39,41:42,44,53,55:61,63,141:143,176,190,195,204:227,229:250,255),]
csballp11 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R11 %in% c(1:6,
             10:14,21:39,41:42,44,53,55:61,63,141:143,176,190,195,204:227,229:250,255),]
csballp12 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R12 %in% c(1:6,
             10:14,21:39,41:42,44,53,55:61,63,141:143,176,190,195,204:227,229:250,255),]
csballp13 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R13 %in% c(1:6,
             10:14,21:39,41:42,44,53,55:61,63,141:143,176,190,195,204:227,229:250,255),]
csballp14 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R14 %in% c(1:6,
             10:14,21:39,41:42,44,53,55:61,63,141:143,176,190,195,204:227,229:250,255),]
csballp15 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R15 %in% c(1:6,
             10:14,21:39,41:42,44,53,55:61,63,141:143,176,190,195,204:227,229:250,255),]
csballp16 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R16 %in% c(1:6,
             10:14,21:39,41:42,44,53,55:61,63,141:143,176,190,195,204:227,229:250,255),]
csballp17 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R17 %in% c(1:6,
             10:14,21:39,41:42,44,53,55:61,63,141:143,176,190,195,204:227,229:250,255),]
csballp18 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R18 %in% c(1:6,
             10:14,21:39,41:42,44,53,55:61,63,141:143,176,190,195,204:227,229:250,255),]
csballp19 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R19 %in% c(1:6,
             10:14,21:39,41:42,44,53,55:61,63,141:143,176,190,195,204:227,229:250,255),]
csballp20 <- csb0916[csb0916$STATEFP %in% st1$STATEFP & csb0916$R20 %in% c(1:6,
             10:14,21:39,41:42,44,53,55:61,63,141:143,176,190,195,204:227,229:250,255),]

# Raster template
r1 <- rast("meantemp_2009.tif")

# Function to convert the shapefiles to rasters
# X is the shapefile, y is the column for the relevant year's crop
host2ras <- function(x, y){
  h2v = vect(x)
  # Change the pixels from crop code to 100 or 0. 
  h2v[h2v > 0] = 100
  h2v = project(h2v, y=crs(r1))
  h2r = rasterize(h2v, r1, field=x[,])
  writeRaster(hr2, <path to save tif>)
}

l1 <- ls(pattern = "host")

# Loop through the hosts in environment
for(i in 1:length(l1)){
  host2ras(l1[i])
}

l2 <- ls(pattern = "allp")

# Loop through the hosts in environment
for(i in 1:length(l2)){
  host2ras(l2[i])
}


