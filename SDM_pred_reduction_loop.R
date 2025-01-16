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

# My stuff!!
# Create a directory for predictors
pred_copies_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/2_Predictors/1_Current/copies/") # nolint
# Crop later
pred_cropped_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/2_Predictors/1_Current/cropped/") # nolint
dir.create(pred_copies_path, showWarnings = FALSE, recursive = TRUE)
dir.create(pred_cropped_path, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(pred_cropped_path, "/transformed"), showWarnings = FALSE,
           recursive = TRUE)

cropped_predictors <- lapply(cropped_predictors, rast)

env_layer <- rast(cropped_predictors)
env_layer <- env_layer[[-25]]

env_group_layer <- data.frame(layer = names(env_layer)[c(1,12:19,2:11,20,37,
                                                         21:22,24,23,36,25:35)])
env_group_layer$group <- c(rep("Temp", 11), rep("Prec", 8), rep("GDD",2), 
                           rep("VI", 3), rep("Human",2), rep("Soil",7), "RH", 
                           rep("Topographic",3))

env_group_layer$TSS

pred_selc <- list()

sdmf <- function(ee) {
  Ind = ee

  # This step creates the partitions for validation. It takes >30 minutes
  prt_lb <- part_sblock(env_layer = env_layer[[ee]], data = lb4, x = "x", 
        y = "y", pr_ab = "pr_ab", n_part = 6, min_res_mult = 4, max_res_mult = 25)
  
  bll2 <- get_block(env_layer, prt_lb$grid)
  
  lb4 <- prt_lb$part
  
  # Background
  bgpt1 <- lapply(1:6, \(x) {
    sample_background(
      data = lb4, x= "x", y ="y", n = (sum(lb4$.part == x)*2) , method = c("thickening",
     width = 6000), rlayer = bll2
    )
  }) %>% bind_rows()
  
  bgpt1$.part <- as.vector(unlist(lapply(1:6, \(x) rep(x, 
       times = (nrow(lb4[lb4$.part == x,]) *2)))))
  
  bgpt1$sp <- "Phytophthora_infestans"
  
  regions1 <- env_layer[[1]]
  regions1[!is.na(regions1)] <- 1
  
  lb5pa <- lapply(1:6, \(x) {
    sample_pseudoabs(
      data = lb4, x= "x", y ="y", n = sum(lb4$.part == x)*2, method = c("geo_const",
                width = "3000"), rlayer = regions1, sp_name = "Phytophthora_infestans"#, 
    )
  }) %>% bind_rows()
  
  lb5pa$.part <- unlist(lapply(1:6, function(x) rep(x, sum(lb4$.part == x)*2)))
  lb5pa <- lb5pa[,c(2:5,1)]
  
  lb4sm$sp <- "Phytophthora_infestans"
  # Combine the presence and pseudo-absence
  lb_allsm <- rbind(lb4sm, lb5pas)
  
  lb_allsm2 <-lb_allsm %>% sdm_extract(data = ., x = "x", y = "y", env_layer = env_cutsm, 
                                       filter_na = T)
  #write.csv(lb_all, "Z:/Late_blight/SDM/lateblight/presenceabsenceallpredictors.csv", row.names = F)
  
  bgpt1s <- bgpt1s[,1:4]
  bgpt1s <- bgpt1s %>% sdm_extract(data = ., x = "x", y = "y", env_layer = env_cutsm,
                                   filter_na = T)
  
  lbnoa <- complete.cases(lb_all)
  bgnoa <- complete.cases(bgpt1)
  lb_all <- lb_all[lbnoa,]
  bgpt1 <- bgpt1[bgnoa,]
  
  glm_run <- fit_glm(data = lb_all, response = "pr_ab", predictors = names(env_layer),
                    partition = ".part", thr = c("max_sens_spec", "equal_sens_spec"),
                    select_pred = F, poly =0, inter_order = 0)
  
  
  return(list(paste0("Ind: ", Ind, "; Layer(s): ", names(env_layer[[ee]]), 
    "; TSS mean: ", glm_run$performance$TSS_mean)))
}

for (i in 1:2) {
  pred_selc[i] = sdmf(i)
}  
