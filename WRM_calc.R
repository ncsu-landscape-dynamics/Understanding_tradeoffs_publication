# This script is called for support in the WC_WRM_cal.R script.

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

# This is the version that Chris wrote. I think it was used at some point, but the most recent coefficient and WRM calculations used
# the version above
#dsv <- function(tmean, rh_hours) {
#  # nested ifel statement
#  # first tmean
#  dsv_val <-
#    ifel(tmean <= 13, 0,
#       # second tmean
#       ifel(tmean > 13 & tmean < 18,
#            ifel(rh_hours <= 6, 0,
#                 ifel(rh_hours > 6 & rh_hours <= 15, 1,
#                      ifel(rh_hours > 15 & rh_hours <= 20, 2, 3))),
#            # third tmean
#            ifel(tmean >= 18 & tmean < 21,
#                 ifel(rh_hours <= 3, 0,
#                      ifel(rh_hours > 3 & rh_hours <= 8, 1,
#                           ifel(rh_hours > 8 & rh_hours <= 15, 2,
#                                ifel(rh_hours > 15 & rh_hours <= 22, 3, 4)))),
#                 # fourth tmean
#                 ifel(tmean >= 21 & tmean < 26,
#                      ifel(rh_hours <= 2, 0,
#                           ifel(rh_hours > 2 & rh_hours <= 5, 1,
#                                ifel(rh_hours > 5 & rh_hours <= 12, 2,
#                                     ifel(rh_hours > 12 & rh_hours <= 20, 3, 4)))),
#                      # fifth tmean
#                      ifel(tmean >= 26,
#                           ifel(rh_hours <= 3, 0,
#                                ifel(rh_hours > 3 & rh_hours <= 8, 1,
#                                     ifel(rh_hours > 5 & rh_hours <= 12, 2,
#                                          ifel(rh_hours > 15 & rh_hours <= 22, 3, 4)))), NA
#                           )))))
#  return(dsv_val)
#}
