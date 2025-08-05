# This is the script that calculates the weather coefficent and WRM for late blight.
# It relies on two separately defined functions. It is parallelized for speed, but can
# still take a day or two to finish, depending on extent and resolution. 

library(terra)
library(doParallel)
source("late_blight_WRM_calc.R")
source("late_blight_temp_hour.R")

rast_size <- nlyr(dayl)

temp_matrix <-
  matrix(c(-100, 3, 0.0,
           3, 4.33, 0.0,
           4.33, 5.66, 0.0,
           5.66, 7.00, 0.0,
           7.00, 8.33, 0.0,
           8.33, 9.66, 0.0,
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

# can change RH threshold depending on disease
rh_threshold <- 90
rh_matrix <- matrix(c(0, rh_threshold, 0, rh_threshold, Inf, 1), 
                    nrow = 2, ncol = 3, byrow = TRUE)

## filepaths and files
inpath <- "..."
outpath <- "..."

tmin_file <- paste0(inpath, "tmin2021.tif")
tmax_file <- paste0(inpath, "tmax2021.tif")

dayl_file    <- paste0(inpath, "daylhr21.tif")
vp_file      <- paste0(inpath, "core_vp_2021.tif")
tmaxb_file   <- paste0(inpath, "weather_core_mxb21.tif")
tmaxb_file   <- paste0(inpath, "weather_core_mxb21.tif")
tmina_file   <- paste0(inpath, "weather_core_mna21.tif")
nightl_file  <- paste0(inpath, "nlhr21.tif")
sunrise_file <- paste0(inpath, "srs21.tif")
sunset_file  <- paste0(inpath, "sst21.tif")
tmean_file   <- paste0(inpath, "weather_avgt_core_2021.tif")
year         <- 2021

filelist <- list.files(outpath, pattern = "wc_sd_")

if (length(filelist) > 0) {
  # creat days from missed days
  splitfront <- strsplit(filelist, paste0("wc_sd_", year, "_"))
  back <- matrix(unlist(splitfront), ncol = 2, byrow = TRUE)[,2]
  mid <- as.matrix(unlist(strsplit(back, ".tif")))
  mid <- as.numeric(mid)
  dayslist <- c(1:365)
  
  "%notin%" <- Negate("%in%")
  
  days <- dayslist[dayslist %notin% mid]
} else {
  days <- c(365)
}

cl <- makeCluster(12)
registerDoParallel(cl)
start_time <- Sys.time()

temdel <- foreach (i = days, .combine = c, .packages = c("terra")) %dopar% {
  dayl <- rast(dayl_file)[[i]]
  tmin <- rast(tmin_file)[[i]]
  tmax <- rast(tmax_file)[[i]]
  vp <- rast(vp_file)[[i]]
  tmaxb <- rast(tmaxb_file)[[i]]
  tmina <- rast(tmina_file)[[i]]
  nightl <- rast(nightl_file)[[i]]
  sunrise <- rast(sunrise_file)[[i]]
  sunset <- rast(sunset_file)[[i]]
  day <- i

  rast_stack <- c(dayl, tmin, tmax, vp, tmaxb, tmina, nightl, sunrise, sunset)
  names(rast_stack) <- c("dayl", "tmin", "tmax", "vp", "tmaxb", "tmina", "nightl", "sunrise", "sunset")
  hourly_temp <- app(rast_stack, function(i, ff) ff(i), ff=tempI)
  hourly_svp <- 611 * 10**(7.5 * hourly_temp / (237.7 + hourly_temp))
  hourly_rh <- vp * 100 / hourly_svp
  hourly_rh_binary <- terra::classify(hourly_rh, rh_matrix)
  hourly_temp_rcl <- terra::classify(hourly_temp, temp_matrix)
  rh_hours <- app(hourly_rh_binary, sum)
  hourly_wc <- hourly_rh_binary*hourly_temp_rcl
  wc <- app(hourly_wc, mean)
  wc_sd <- app(hourly_wc, sd)
  
  tmean <- rast(tmean_file)[[i]]
  calc_rast <- c(tmean, rh_hours)
  names(calc_rast) <- c("tmean", "rh_hours")

  dsv_score <- app(calc_rast, function(i, ff) ff(i), ff=dsv)
  writeRaster(wc, paste0(outpath, "wc_", year, "_", day, ".tif"))
  writeRaster(wc_sd, paste0(outpath, "wc_sd_", year, "_", day, ".tif"))
  writeRaster(dsv_score, paste0(outpath, "dsv_score_", year, "_", day, ".tif"))
  writeRaster(calc_rast, paste0(outpath, "calc_rast_", year,"_",day,".tif"))
  return(1)
}

end_time <- Sys.time()
time_taken <- end_time - start_time
stopCluster(cl)
