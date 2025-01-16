library(flexsdm)
library(terra)
library(data.table)
library(parallel)
library(spatialEco)

set.seed(42)

path <- "Z:/pops_pesthostuse/pops.sdm/Data/"
setwd("Z:/Late_blight/SDM/lateblight/")
domain <- "USA"
extent <- "Z:/Late_blight/helpful_shapes_rasters/lb_extent.gpkg"

res <- 1000

lb2 <- vect("Z:/Late_blight/helpful_shapes_rasters/blightptsmovedforSDM.gpkg")
lbdf <- as.data.frame(lb2[,c(1,16:17)])
lbdf$id <- paste0(lbdf$x,"_",lbdf$y)
lbdf <- lbdf[!duplicated(lbdf$id),]
lb2sel <- lb2[lb2$Rec_ID %in% lbdf$Rec_ID,]

# This splits the USA blight points up for testing and validating. Gonna skip
# for now.
lbtn <- sample(lb2sel, nrow(lb2sel)*.8)
lbvl <- lb2sel[!lb2sel$Rec_ID %in% lbtn$Rec_ID,]

lb3 <- as.data.frame(lbtn[,c(1,16,17)])

#' 2. Set up SDM directories
flexsdm::sdm_directory(
  algorithm = TRUE
)

subset_files <- list.files("Z:/Late_blight/SDM/lateblight/flexsdm_results/1_Inputs/2_Predictors/1_Current/cropped/transformed/", full.names=T)
subset_files <- subset_files[-c(1)]

cropped_predictors <- lapply(subset_files, rast)

env_layer <- rast(cropped_predictors)
env_layer <- env_layer[[-25]]

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

lb_all <- read.csv("Z:/Late_blight/SDM/lateblight/presenceabsenceallpredictors.csv")

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

for (i in 1:2) {
  pred_selc[i] = sdmf(i)
}  
