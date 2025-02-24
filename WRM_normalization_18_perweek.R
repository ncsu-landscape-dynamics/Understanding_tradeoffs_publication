# Create a new mosaic of counts of accumulations of 18. 
# Divide by week to see how many go over 1/week in a season&zone.
# These are rasters with cumulative daily sums, in order to go period by period.
cl1 <- list.files("Z:/Late_blight/lbdsveconu/subs", pattern = "cs", full.names = T)
cl1 <- gtools::mixedsort(cl1)


# Consumes too much memory
cnt18f <- function(x) {
  rt1 = rast(x)
  rn0 = rt1
  rn0[rn0 > 0] <- 0
  nl1 = nlyr(rt1)
  for (n in 1:nl1) {
    rn0[[n]] = floor(rt1[[n]] / 18)
  }
  return(rn0)
}

# USE this
cnt18f <- function(x) {
  rt1 = rast(x)
  nl1 = nlyr(rt1)
  # Get the zone from the file name
  yn1 = gsub("\\.", "", str_sub(x, 44, 45))
  weekseqmax = floor(nl1/7)
  # Sequence of starting days for each week for given no. of layers
  stseq = seq(1, length.out = weekseqmax, by = 7)
  # Sequence to conclude each
  # Change 7 to 14 for values for every two weeks.
  enseq = seq(7, length.out = weekseqmax, by = 7)
  rt2 = rt1[[1:weekseqmax]]
  rt2[rt2 > 0] = 0
  for (i in 1:weekseqmax) {
    print(stseq[i])
    print(enseq[i])
    rs1 = rt1[[stseq[i]:enseq[i]]]
    rt2[[i]] = floor(rt1[[enseq[i]]] / 18) - floor(rt1[[stseq[i]]] / 18)
    rt2[rt2 <= 1] = 0
  }
  writeRaster(rt2, paste0("Z:/Late_blight/lbdsveconu/subs/wrm_count18week_2021_hz",
                          yn1,".tif"))
}
