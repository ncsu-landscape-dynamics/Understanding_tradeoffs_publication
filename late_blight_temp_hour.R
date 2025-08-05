# This is a supporting function for the weather coefficient processing.
# It uses all the weather and time related inputs with a couple of constants
# to calculate hourly temperature values. The time related inputs are based on
# latitude. For more information see Goudriaan and van Laar 1994.

tempI <- function(rast_stack) {
  dayl <- rast_stack["dayl"]
  tmin <- rast_stack["tmin"]
  tmax <- rast_stack["tmax"]
  vp <- rast_stack["vp"]
  tmaxb <- rast_stack["tmaxb"]
  tmina <- rast_stack["tmina"]
  nightl <- rast_stack["nightl"]
  sunrise <- rast_stack["sunrise"]
  sunset <- rast_stack["sunset"]
  # hour <- rast_stack["hour"]
  hours <- seq(1:24)
  out <- rep(0, 24)
  
  TC <- 4.0
  P <- 1.5
  time_bin1 <- 13.5
  for (hour in hours) {
    if (is.na(sunrise) || is.na(sunset)) {
      out[hour] <- NA
    } else if (hour <= sunrise) {
      # hours between 0 and sunrise calc for TSunset is (tmin + (tmaxb - tmin) * sin(pi*(dayl / (dayl + 2 * P)))) before sunrise
      out[hour] <- ((tmin - (tmin + (tmaxb - tmin) * sin(pi*(dayl / (dayl + 2 * P)))) * exp( -nightl / TC) + ((tmin +(tmaxb - tmin) * sin(pi * (dayl / (dayl + 2*P)))) - tmin) *
                 exp( -(hour + 24 - sunset) / TC)) / (1.0 - exp( -nightl / TC)))
    } else if (hour <= time_bin1) {
      # hour between sunrise and time_bin1 which is 13.5
      out[hour] <- (tmin + (tmax - tmin) * sin( pi * (hour - sunrise) / (dayl + 2 * P)))
    } else if (hour <= sunset) {
      # hour between time_bin1 which is 13.5 and sunset
      out[hour] <- (tmina + (tmax - tmina) * sin( pi * (hour - sunrise) / (dayl + 2 * P)))
    } else {
      # hours between sunset and 24 TSunset after sunset is (tmina + (tmax - tmina) * sin(pi * (dayl / (dayl + 2 * P))))
      out[hour] <- ((tmina - (tmina + (tmax - tmina) * sin(pi * (dayl / (dayl + 2 * P)))) * exp(-nightl / TC) +
                 ((tmin + (tmax - tmin) * sin(pi * (dayl / (dayl + 2 * P))))
                  - tmin) * exp(-(hour - sunset) / TC)) / (1 - exp(-nightl / TC)))
    }
  }

  return(out)
}
