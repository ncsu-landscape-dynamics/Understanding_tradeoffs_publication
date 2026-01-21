library(terra)

 
location_reduction <- function(cells, reduce_percent, inst) {
  r_tem = red_raster
  # Set original infections to background
  r_tem[cells2avoid] = NA
  # Pick reduce_percent% of cells
  smlset = sample(cells, round(length(cells)*reduce_percent))
  print(smlset)
  print(paste0("reduce_percent: ", reduce_percent))
  print(length(smlset))
  # Set selected cells to same value as background
  r_tem[smlset] = NA
  
  pct = 100 - (reduce_percent*100)
  #print(paste0(outputpath, pct, "/mean22_", pct, "-", inst, "round.tif"))
  print(paste0(outputpath, pct, "/TEST2_", pct, "-", inst, ".tif"))
  #writeRaster(r_tem, paste0(outputpath, pct, "/mean22_", pct, "-", inst, 
  #                          "round.tif"), overwrite=T)
  writeRaster(r_tem, paste0(outputpath, "/TEST2_", pct, "-", inst, #add "pct, 
                            ".tif"), overwrite=T)
                            
}

