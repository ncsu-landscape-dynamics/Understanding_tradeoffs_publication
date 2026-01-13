library(terra)

 
location_reduction <- function(cells, reduce_percent, inst) {
  r_tem = red_raster
  # Pick reduce %age of cells
  smlset = sample(cells, round(length(cells)*reduce_percent))
  
  print(paste0("reduce percent: ", reduce_percent))
  print(length(smlset))
  # Set selected cells to same value as background
  r_tem[smlset] = 0
  
  pct = 100 - (reduce_percent*100)
  print(paste0(outputpath, pct, "/TEST_", pct, "-", inst, "round.tif"))

  writeRaster(r_tem, paste0(outputpath, pct, "/TEST_", pct, "-", inst, 
                            "round.tif"), overwrite=T)

}


