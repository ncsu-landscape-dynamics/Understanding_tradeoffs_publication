library(terra)

 
location_reduction <- function(cells, x, inst) {
  r_tem = red_raster
  # Pick X% of cells
  smlset = sample(cells, round(length(cells)*x))
  
  print(paste0("x: ", x))
  print(length(smlset))
  # Set selected cells to same value as background
  r_tem[smlset] = 0
  
  pct = 100 - (x*100)
  print(paste0(outputpath, pct, "/TEST_", pct, "-", inst, "round.tif"))

  writeRaster(r_tem, paste0(outputpath, pct, "/TEST_", pct, "-", inst, 
                            "round.tif"), overwrite=T)

}


