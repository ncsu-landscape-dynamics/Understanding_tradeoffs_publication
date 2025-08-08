library(terra)
library(tidyverse)

# Coordinates from updated plant aid 
usdf <- read.csv("Z:/Late_blight/ms3/usabl2025.csv")
usdf$observation_date <- as_date(usdf$observation_date)
usdf$year <- year(usdf$observation_date)
usdf <- usdf %>% filter(year %in% 2021:2023)
usve <- vect(usdf, geom = c("longitude","latitude"), crs = "epsg:4326", keepgeom =T)

# New area for the simulation, the "core"
lb_core <- vect("Z:/Late_blight/helpful_shapes_rasters/lb_centralregion.geojson")
lb_core <- project(lb_core, y = crs(usve))

lb_locs_sim <- crop(usve, lb_core)

# Get a raster similar to established core
temperature_coefficient_file <- "Z:/Late_blight/Weather_data/lb_2021_coef_avg.tif"
templ1 <- rast(temperature_coefficient_file)
templ2 <- templ1[[1]]
templ2[!is.na(templ2)] <- 0
lb_locs_sim <- project(lb_locs_sim, y = crs(templ2))
lb_locs_sim$longitude <- crds(lb_locs_sim)[,1]
lb_locs_sim$latitude <- crds(lb_locs_sim)[,2]

# Year 1
lb21 <- lb_locs_sim[lb_locs_sim$year == 2021,]
mx1 <- matrix(c(lb21$longitude, lb21$latitude), ncol = 2, byrow = F)

# Year 2
lb22 <- lb_locs_sim[lb_locs_sim$year == 2022,]
mx2 <- matrix(c(lb22$longitude, lb22$latitude), ncol = 2, byrow = F)

# Year 3
lb23 <- lb_locs_sim[lb_locs_sim$year == 2023,]
mx3 <- matrix(c(lb23$longitude, lb23$latitude), ncol = 2, byrow = F)

# Add points to template
loc21 <- templ2
loc22 <- templ2
loc23 <- templ2
loc21[cellFromXY(loc21, mx1)] <- 1
loc22[cellFromXY(loc22, mx2)] <- 1
loc23[cellFromXY(loc23, mx3)] <- 1

writeRaster(loc21, "Z:/Late_blight/simulation/LB/infection_2021.tif")
writeRaster(loc22, "Z:/Late_blight/simulation/LB/infection_2022.tif")
writeRaster(loc23, "Z:/Late_blight/simulation/LB/infection_2023.tif")




