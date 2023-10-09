library(terra)
library(doParallel)

dayl <- rast("~/Google Drive/Shared drives/John_Polo_Dissertation/MS_1/data/lateblight/rasterinputs/dayl09.tif")
rast_size <- nlyr(dayl)
rm(dayl)

temp_matrix <-
  matrix(c(-100, 3, 0.0,
           3, 4.33, 0.1,
           4.33, 5.66, 0.2,
           5.66, 7.00, 0.3,
           7.00, 8.33, 0.4,
           8.33, 9.66, 0.5,
           9.66, 11.00, 0.6,
           11.00, 12.33, 0.7,
           12.33, 13.66, 0.8,
           13.66, 15.00, 0.9,
           15.0, 18, 1.0,
           18, 18.77, 0.9,
           18.77, 19.55, 0.8,
           19.55, 20.33, 0.7,
           20.33, 21.11, 0.6,
           21.11, 21.88, 0.5,
           21.88, 22.66, 0.4,
           22.66, 23.44, 0.3,
           23.44, 24.22, 0.2,
           24.22, 25, 0.1,
           25, 100, 0.0),
         ncol = 3, byrow = TRUE)

rh_threshold <- 90
rh_matrix <- matrix(c(0, rh_threshold, 0, rh_threshold, Inf, 1), nrow = 2, ncol = 3, byrow = TRUE)

outpath <- "~/Google Drive/Shared drives/John_Polo_Dissertation/MS_1/data/lateblight/rasterinputs/outputs/"
year <- 2009

cl <- makeCluster(7) # on Darwin
registerDoParallel(cl)
# Not sure if I'm supposed to use the .combine = c... tried a smaller example without and it created list of rasters.
temdel <- foreach (i = 1:7, .combine = c, .packages = c("terra")) %dopar% {
  dayl <- rast("~/Google Drive/Shared drives/John_Polo_Dissertation/MS_1/data/lateblight/rasterinputs/dayl09_hour.tif")[[i]]
  tmin <- rast("~/Google Drive/Shared drives/John_Polo_Dissertation/MS_1/data/lateblight/rasterinputs/tmin09.tif")[[i]]
  tmax <- rast("~/Google Drive/Shared drives/John_Polo_Dissertation/MS_1/data/lateblight/rasterinputs/tmax09.tif")[[i]]
  vp <- rast("~/Google Drive/Shared drives/John_Polo_Dissertation/MS_1/data/lateblight/rasterinputs/vp09.tif")[[i]]
  tmaxb <- rast("~/Google Drive/Shared drives/John_Polo_Dissertation/MS_1/data/lateblight/rasterinputs/tmaxb09.tif")[[i]]
  tmina <- rast("~/Google Drive/Shared drives/John_Polo_Dissertation/MS_1/data/lateblight/rasterinputs/tmina09.tif")[[i]]


  zero_rast <- dayl
  values(zero_rast) <- 0

  #create time_bin raster for use in the ifel statement
  time_bin1 <- zero_rast
  values(time_bin1) <- 13.5

  # computer nightlenght, sunrise and sunset
  nightl <- 24 - dayl
  sunrise <- 12 - 0.5 * dayl
  sunset <- 12 + 0.5 * dayl

  number_hours_above_rh_threshold <-
    wc_and_dsv(dayl, tmin, tmax, vp, tmaxb, tmina, time_bin1, nightl, sunrise,
                sunset, zero_rast, rh_matrix, temp_matrix, outpath, year, day = i)
  list(i, number_hours_above_rh_threshold)
}

stopCluster(cl)
