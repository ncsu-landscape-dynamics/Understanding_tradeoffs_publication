# Biolcim variables
#BIO1 = Annual Mean Temperature
#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO3 = Isothermality (BIO2/BIO7) (×100)
#BIO4 = Temperature Seasonality (standard deviation ×100)
#BIO5 = Max Temperature of Warmest Month
#BIO6 = Min Temperature of Coldest Month
#BIO7 = Temperature Annual Range (BIO5-BIO6)
#BIO8 = Mean Temperature of Wettest Quarter
#BIO9 = Mean Temperature of Driest Quarter
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter
#BIO12 = Annual Precipitation
#BIO13 = Precipitation of Wettest Month
#BIO14 = Precipitation of Driest Month
#BIO15 = Precipitation Seasonality (Coefficient of Variation)
#BIO16 = Precipitation of Wettest Quarter
#BIO17 = Precipitation of Driest Quarter
#BIO18 = Precipitation of Warmest Quarter
#BIO19 = Precipitation of Coldest Quarter

theta_s 	saturated soil water content 	m3/m3
theta_r 	residual soil water content 	m3/m3

#' 1. Load the required packages & set path
library(flexsdm)
library(terra)
library(data.table)
library(parallel)
library(spatialEco)

set.seed(42)

path <- "//pops.sdm/Data/"
domain <- "USA"
# My late blight extent
extent <- "Z:/Late_blight/helpful_shapes_rasters/lb_extent.gpkg"
# Another extent with which to crop
crex1 <- vect(ext(-82,-65.9042, 27, 47.5), crs="epsg:4326")
#crex1 <- project(crex1, y = "epsg:3857")
extent2 <- terra::intersect(vect(extent), crex1)

# The relative humidity data for late blight. It's the mean of the USDA growing
# zones during the tomato growing seasons from multiple years. Add this to 
# env_layers
rhmean <- rast("Z:/Late_blight/SDM/lateblight/rhmeanall.tif")

res <- 1000

lb2 <- vect("Z:/Late_blight/helpful_shapes_rasters/blightptforSDM.gpkg")
lbdf <- as.data.frame(lb2[,c(1,16:17)])
lbdf$id <- paste0(lbdf$x,"_",lbdf$y)
lbdf <- lbdf[!duplicated(lbdf$id),]
lb2sel <- lb2[lb2$Rec_ID %in% lbdf$Rec_ID,]

# Break up the lb2sel into two regions for calibration and validation.
crex1 <- project(crex1, y = crs(lb2sel))
lbcr1 <- terra::intersect(lb2sel, crex1)
lbv1 <- lb2sel[!lb2sel$Rec_ID %in% lbcr1$Rec_ID,]

# This splits the USA blight points up for testing and validating. Gonna skip
# for now.
lbtn <- sample(lb2sel, nrow(lb2sel)*.8)
lbvl <- lb2sel[!lb2sel$Rec_ID %in% lbtn$Rec_ID,]

lb3 <- as.data.frame(lbtn[,c(1,16,17)])
# Buffers for calibration area
lbbuf5k <- buffer(lb2, 5000) 
lbbuf5k <- calib_area(data = lbtn, x = "x", y = "y", method = c("buffer", width =5000), crs = "epsg:4326")

#' 2. Set up SDM directories
flexsdm::sdm_directory(
  algorithm = TRUE
)

#' 3. Create base raster for the study extent
source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/sdm_helpers.R")
source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/get_envi_chunked.R") # get_topo_global # nolint
source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/raster_base.R") # base_raster at correct resolution and extent # nolint
source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/part_sblock.R") # nolint
source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/sample_background.R") # nolint

res <- fix_resolution(res, domain)
base <- crop_base_raster(domain, res, path, extent2)
domain <- "USA"

# Copy predictors, crop to the study extent, and transform to z-scores
predictors <- subset_rasters(path, domain = domain, all = FALSE, gdd = TRUE, 
  roads = TRUE, rails = FALSE, pop = TRUE, downscaled = FALSE)

files <- as.vector(unlist(lapply(predictors, function(x) {
  get_filename(x, "USA", path)
})))


# Retain only the rasters at the correct resolution
subset_files <- files[grep(paste0(res, "m"), files)]
subset_files <- list.files("Z:/Late_blight/SDM/lateblight/flexsdm_results/1_Inputs/2_Predictors/1_Current/cropped/transformed/", 
  full.names=T)
subset_files <- subset_files[-c(1,26,34)]

# Create a directory for predictors
pred_copies_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/2_Predictors/1_Current/copies/") # nolint
# Crop later
pred_cropped_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/2_Predictors/1_Current/cropped/") # nolint
dir.create(pred_copies_path, showWarnings = FALSE, recursive = TRUE)
dir.create(pred_cropped_path, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(pred_cropped_path, "/transformed"), showWarnings = FALSE,
           recursive = TRUE)

# Copy and crop the predictors (SDM helper)
cropped_predictors <- list()

for (i in seq_along(subset_files)) {
  cropped_predictors[[i]] <- crop_predictor(
    file = subset_files[i],
    pred_copies_path = pred_copies_path,
    pred_cropped_path = pred_cropped_path,
    extent = extent
  )
}

# Combine predictors and RH
cropped_predictors <- c(cropped_predictors, rhm2)

crp_prd <- rast(cropped_predictors)

# Load the species occurrence data
#species <- "Phytophthora infestans"
#species <- format_species_name(species)


# Geographic filtering of occurrence data to address spatial autocorrelation - READ flxdm
cellsizes <- get_geo_cellsizes(domain = "USA", res = 1000) #<- Brit's helper
filtgeo_files <- paste0(getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/filtered/",
                        species, "_occ_pa_geofilter_", cellsizes, "m2.csv")

# Isolate numerical predictors for spatial block partitioning
env_layer <- list()

for (i in seq_along(cropped_predictors)) {
  env_layer[[i]] <- get_numerical_rasters(cropped_predictors[[i]])
}

env_layer <- rast(env_layer)
names(env_layer)[32] <- "RH"

regions1 <- env_layer[[1]]
regions1[!is.na(regions1)] <- 1

lb3b <- lb3[,c(1,2,4)]

psu1 <- sample_pseudoabs(data = lb4, x = "x", y = "y", n = nrow(lb4), 
                         method = c("geo_const", width = "5000"),#, env = crp_prd), 
                         rlayer = regions1, sp_name = "Phytophthora_infestans")

set.seed(42)
lb5pa <- lapply(1:6, \(x) {
  sample_pseudoabs(
    data = lb4, x= "x", y ="y", n = sum(lb4$.part == x)*2, method = c("geo_const",
                                                                      width = "3000"), rlayer = regions1, sp_name = "Phytophthora_infestans"#, 
    #calibarea = 
  )
}) %>% bind_rows()

lb5pa$.part <- unlist(lapply(1:6, function(x) rep(x, sum(lb4$.part == x)*2)))
lb5pa <- lb5pa[,c(2:5,1)]

lb4$sp <- "Phytophthora_infestans"
# Combine the presence and pseudo-absence
lb_all <- rbind(lb4, lb5pa)

# Extract environmental data
lb_all <-lb_all %>% sdm_extract(data = ., x = "x", y = "y", env_layer = env_layer, 
                                filter_na = T)

bgpt1 <- bgpt1 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = env_layer,
                               filter_na = T)

lbnoa <- complete.cases(lb_all)
bgnoa <- complete.cases(bgpt1)
lb_all <- lb_all[lbnoa,]
bgpt1 <- bgpt1[bgnoa,]

tmax1 <- tune_max(data = lb_all, response = "pr_ab", predictors = names(env_layer),
  background = bgpt1, partition = ".part", grid = expand.grid(regmult = seq(0.1, 
  3, 0.5), classes = c("l", "lq", "lqhpt"), thr = c("max_sens_spec", 
  "equal_sens_spec", "max_sorensen"), metric = "TSS", clamp = TRUE, pred_type = "cloglog"))

glm1 <- fit_glm(data = lb_all, response = "pr_ab", predictors = names(env_layer),
  partition = ".part", thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
  select_pred = F, poly =0, inter_order = 0)

lb_ens <- fit_ensemble(models = list(tmax1, glm1),
                       ens_method = "meanw",
                       thr = c("equal_sens_spec", "max_sens_spec"),
                       thr_model = "max_sens_spec",
                       metric = "TSS")

# Another take on the MaxEnt
em1 <- esm_max(data = lb_all, response = "pr_ab", predictors = names(env_layer), 
               partition = ".part", thr = c("equal_sens_spec", "max_sens_spec"), background = bgpt1)

val_area <- rast(lapply(subset_files, rast))
prd2 <- sdm_predict(models = lb_ens, pred = env_layer, thr = "max_sens_spec",
                    predict_area = val_area)

