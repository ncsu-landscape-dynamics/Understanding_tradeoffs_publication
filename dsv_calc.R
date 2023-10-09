dsv <- function(tmean, rh_hours) {
  # nested ifel statement
  # first tmean
  dsv_val <-
    ifel(tmean <= 13, 0,
       # second tmean
       ifel(tmean > 13 & tmean < 18,
            ifel(rh_hours <= 6, 0,
                 ifel(rh_hours > 6 & rh_hours <= 15, 1,
                      ifel(rh_hours > 15 & rh_hours <= 20, 2, 3))),
            # third tmean
            ifel(tmean >= 18 & tmean < 21,
                 ifel(rh_hours <= 3, 0,
                      ifel(rh_hours > 3 & rh_hours <= 8, 1,
                           ifel(rh_hours > 8 & rh_hours <= 15, 2,
                                ifel(rh_hours > 15 & rh_hours <= 22, 3, 4)))),
                 # fourth tmean
                 ifel(tmean >= 21 & tmean < 26,
                      ifel(rh_hours <= 2, 0,
                           ifel(rh_hours > 2 & rh_hours <= 5, 1,
                                ifel(rh_hours > 5 & rh_hours <= 12, 2,
                                     ifel(rh_hours > 12 & rh_hours <= 20, 3, 4)))),
                      # fifth tmean
                      ifel(tmean >= 26,
                           ifel(rh_hours <= 3, 0,
                                ifel(rh_hours > 3 & rh_hours <= 8, 1,
                                     ifel(rh_hours > 5 & rh_hours <= 12, 2,
                                          ifel(rh_hours > 15 & rh_hours <= 22, 3, 4)))), NA
                           )))))
  return(dsv_val)
}