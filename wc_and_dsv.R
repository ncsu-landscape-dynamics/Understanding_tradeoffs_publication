
wc_and_dsv <-
  function(dayl, tmin, tmax, vp, tmaxb, tmina, time_bin1, nightl, sunrise, sunset, zero_rast, rh_matrix, temp_matrix, outpath, year, day) {
  TC <- 4.0
  P <- 1.5
  hours <- seq(1:24)

  hourly_temp <- rast(zero_rast, nlyr = 24)
  for (hour in hours) {
    hourly_temp[[hour]] <-
      terra::ifel(hour < sunrise,
                  # hours between 0 and sunrise calc for TSunset is (tmin + (tmaxb - tmin) * sin(pi*(dayl / (dayl + 2 * P)))) before sunrise
                  ((tmin - (tmin + (tmaxb - tmin) * sin(pi*(dayl / (dayl + 2 * P)))) * exp( -nightl / TC) + ((tmin +(tmaxb - tmin) * sin(pi * (dayl / (dayl+ 2*P)))) - tmin) *
                      exp( -(hour + 24 - sunset) / TC)) / (1.0 - exp( -nightl / TC))),
                  terra::ifel(hour < time_bin1,
                              # hour between sunrise and time_bin1 which is 13.5
                              (tmin + (tmax- tmin) * sin( pi * (hour - sunrise) / (dayl + 2 * P))),
                              terra::ifel(hour < sunset,
                                          # hour between time_bin1 which is 13.5 and sunset
                                          (tmina + (tmax - tmina) * sin( pi * (hour - sunrise) / (dayl + 2 * P))),
                                          #hours between sunset and 24 TSunset after sunset is (tmina + (tmax - tmina) * sin(pi * (dayl / (dayl + 2 * P))))
                                          ((tmina - (tmina + (tmax - tmina) * sin(pi * (dayl / (dayl + 2 * P)))) * exp(-nightl / TC) +
                                              ((tmin + (tmax - tmin) * sin(pi * (dayl / (dayl + 2 * P))))
                                               - tmin) * exp(-(hour - sunset) / TC)) / (1 - exp(-nightl / TC))))))
  }

  hourly_svp <- 611 * 10**(7.5 * hourly_temp / (237.7 + hourly_temp))
  hourly_rh <- vp * 100 / hourly_svp
  hourly_rh_binary <- terra::classify(hourly_rh, rh_matrix)
  hourly_temp_rcl <- terra::classify(hourly_temp, temp_matrix)
  number_hours_above_rh_threshold <- app(hourly_rh_binary, sum)
  hourly_wc <- hourly_rh_binary*hourly_temp_rcl
  wc <- app(hourly_wc, mean)
  wc_sd <- app(hourly_wc, sd)

  tmean <- (tmax + tmin) / 2
  dsv_score <- dsv(tmean, number_hours_above_rh_threshold)

  writeRaster(wc, paste0(outpath, "wc_", year, "_", day, ".tif"))
  writeRaster(wc_sd, paste0(outpath, "wc_sd_", year, "_", day, ".tif"))
  writeRaster(dsv_score, paste0(outpath, "dsv_score_", year, "_", day, ".tif"))
  # writeRaster(hourly_temp, paste0(outpath, "hourly_temp_", year, "_", day, ".tif"))
  # writeRaster(hourly_rh, paste0(outpath, "hourly_rh_", year, "_", day, ".tif"))
  success <- TRUE
  return(success)
}
