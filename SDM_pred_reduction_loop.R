library(flexsdm)
library(terra)
library(data.table)
library(parallel)
library(spatialEco)

set.seed(42)

sdmpath <- </path/to/sdm/predictors>

sdm_files <- list.files(sdmpath, full.names=T)
sdm_files <- sdm_files[-1] #remove a folder from the list of files

env_layer <- rast(sdm_files)

# Create a data frame that has the layer name, the group, and column for inserting the TSS from first, solo model
env_group_layer <- data.frame(layer = names(env_layer)[c(1,12:19,2:11,20,37,
                                                         21:22,24,23,36,25:35)])
env_group_layer$group <- c(rep("Temp", 11), rep("Prec", 8), rep("GDD",2), 
                           rep("VI", 3), rep("Human",2), rep("Soil",7), "RH", 
                           rep("Topographic",3))

# Match order of rows in data frame to the order of layers in raster
nm1 <- names(env_layer)
env_group_layer <- env_group_layer[match(nm1, env_group_layer$layer),]
rownames(env_group_layer) <- 1:nrow(env_group_layer)

env_group_layer$TSS <- NA
env_group_layer$indx <- 1:nrow(env_group_layer)

pred_selc <- list()

# lb_all is created in the runSDM_<model> script.
# It terra::extracts all the predictors to using the presence and absence points created in the script.
lb_all <- read.csv("Z:/Late_blight/SDM/lateblight/presenceabsenceallpredictors.csv")
# Sometimes the NLCD layer causes issues. I just drop it when it does.
lb_all <- lb_all[, names(lb_all) != "NLCD.Land.Cover.Class"]

# Dropped several lines. Reusing a data frame that may or may not be appropriate for this package's process.
sdmf <- function(ee) {
  Ind = ee
  
  cols2use <- c(1:5, 5+ee)
  
  lb_temp <- lb_all[,c(cols2use)]
  
  glm_run <- fit_glm(data = lb_temp, response = "pr_ab", predictors = names(env_layer[[ee]]),
                     partition = ".part", thr = c("max_sens_spec"),
                     select_pred = F, poly =0, inter_order = 0)
  
  
  return(list(paste0("Ind: ", Ind, "; Layer(s): ", names(env_layer[[ee]]), 
                     "; TSS mean: ", glm_run$performance$TSS_mean)))
}

for (i in 1:nlyr(env_layer)) {
  pred_selc[i] = sdmf(i)
}  
