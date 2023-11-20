dsv <- function(calc_rast) {
  tmean <- calc_rast["tmean"]
  rh_hours <- calc_rast["rh_hours"]
  
  # return 0 if arguments are NA
  if (is.na(tmean) || is.na(rh_hours)) {
    return(0)
  } else if (tmean > 13 && tmean < 18) {
    if (rh_hours <= 6) return(0)
    if (rh_hours <= 15) return(1)
    if (rh_hours <= 20) return(2)
    return(3)
  } else if (tmean >= 18 && tmean < 21) {
    if (rh_hours <= 3) return(0)
    if (rh_hours <= 8) return(1)
    if (rh_hours <= 15) return(2)
    if (rh_hours <= 22) return(3)
    return(4)
  } else if (tmean >= 21 && tmean < 26) {
    if (rh_hours <= 2) return(0)
    if (rh_hours <= 5) return(1)
    if (rh_hours <= 12) return(2)
    if (rh_hours <= 20) return(3)
    return(4)
  } else if (tmean >= 26) {
    if (rh_hours <= 3) return(0)
    if (rh_hours <= 8) return(1)
    if (rh_hours <= 15) return(2)
    if (rh_hours <= 22) return(3)
    return(4)
  }
  
  return(0)
}
