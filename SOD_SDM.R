# 1. Load the required packages & set path
library(flexsdm)
library(terra)
library(data.table)
library(parallel)
library(spatialEco)

path <- "/path/to/sdm/predictor/data/"
domain <- "USA"
# Extent
or1 <- vect("Z:/Late_blight/helpful_shapes_rasters/usa2/USA_adm2.shp")
or1 <- or1[or1$NAME_1 == "Oregon"& or1$NAME_2 %in% c("Coos","Curry","Douglas","Josephine"),]
or1w <- project(or1, y = "epsg:4326")
sodpp <- rast("Z:/Late_blight/SOD_OR/mean5thyear.tif")
sodpg <- project(sodpp, y = "epsg:4326", method = "near")
extent <- as.polygons(crop(sodpg, or1w, mask =T))
#writeVector(extent, "Z:/Late_blight/SDM/SOD/sodextent.gpkg")
#extent <- "Z:/Late_blight/helpful_shapes_rasters/sodextent.gpkg"

res <- 100

# Infections
sod1 <- vect("Z:/Late_blight/helpful_shapes_rasters/SOD_all.gpkg")

# Set up SDM directories
flexsdm::sdm_directory(
  algorithm = TRUE
)

# Create base raster for the study extent
source("sdm_helpers.R")
source("get_envi_chunked.R")
source("raster_base.R")
source("part_sblock.R")
source("sample_background.R")

res <- fix_resolution(res, domain)
base <- crop_base_raster(domain, res, path, extent)

predictors <- subset_rasters(path =path, domain = domain, all = FALSE, gdd = TRUE, 
  roads = TRUE, rails = FALSE, pop = TRUE, downscaled = FALSE)

files <- as.vector(unlist(lapply(predictors, function(x) {
  get_filename(x, "USA", path)
})))

# Retain only the rasters at the correct resolution
subset_files <- files[grep(paste0(res, "m"), files)]
subset_files <- list.files("Z:/Late_blight/SDM/SOD/flexsdm_results/1_Inputs/2_Predictors/1_Current/cropped/transformed/", 
  full.names = T)

# Create a directory for predictors
pred_copies_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/2_Predictors/1_Current/copies/") # nolint
# Crop later
pred_cropped_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/2_Predictors/1_Current/cropped/") # nolint
dir.create(pred_copies_path, showWarnings = FALSE, recursive = TRUE)
dir.create(pred_cropped_path, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(pred_cropped_path, "/transformed"), showWarnings = FALSE,
           recursive = TRUE)

cropped_predictors <- lapply(subset_files, rast)
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

crp_prd <- rast(cropped_predictors)

# Load the species occurrence data
species <- "Phytophthora ramorum"
species <- format_species_name(species)

# Geographic filtering of occurrence data to address spatial autocorrelation - READ flxdm
cellsizes <- get_geo_cellsizes(domain = "USA", res = 1000) #<- Brit's helper
filtgeo_files <- paste0(getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/filtered/",
                        species, "_occ_pa_geofilter_", cellsizes, "m2.csv")

# Create a spatial-block partition for the data based on continuous predictors READ flxdm
# Isolate numerical predictors for spatial block partitioning
env_layer <- list()

for (i in seq_along(cropped_predictors)) {
  env_layer[[i]] <- get_numerical_rasters(cropped_predictors[[i]])
}

env_layer <- rast(env_layer)
env_layer <- rast(cropped_predictors)

# Best predictors based on iterative GLM tests.
env_rd <- env_layer[[names(env_layer) %in% c("bio1", "bio15","par", "gdd",
   "dist2roads","ndvi","theta_r_0_5","clay_0_5", "ph_0_5", "elevation")]]

# Get the infections to right projection
sod1 <- project(sod1, y = crs(rast(cropped_predictors))) # The y is from later..
sod1$x <- crds(sod1)[,1]
sod1$y <- crds(sod1)[,2]
sod2 <- as.data.frame(sod1[,c(1,6,7)])
sod3 <- sod2[,c(2:3)]
sod3$pr_ab <- 1
sod3$id <- 1:nrow(sod3)

# Partition the filtered data into spatial blocks & write out
part_dt_files <- paste0(getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/partitioned/",
                            species, "_part_data_filtgeo", cellsizes, "m2.csv")
part_raster_files <- paste0(getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/partitioned/",
                            species, "_part_grid_filtgeo", cellsizes, "m2.tif")

partitioned_dt <- list()
partitioned_raster <- list()
partitioned_data <- list()

sod3b <- sod3[,1:3]

prt_sod <- part_random(data = sod3b, pr_ab = "pr_ab", method = c(method = "kfold",
  folds = 4))

sod4 <- prt_sod

rly1 <- env_layer[[1]]
rly1[!is.na(rly1)] <- 1

bgpt1 <- sample_background(
    data = sod3b, x= "x", y ="y", n = nrow(sod3b), method = "random", 
    rlayer = rly1#, calibarea = lbbufmcp#, 
    #calibarea = 
  ) %>% bind_rows()

prt_sodb <- part_random(data = bgpt1, pr_ab = "pr_ab", method = c(method="kfold", 
  folds = 4))

sod4$species <- "Phytophthora_ramorum"

psu1 <- sample_pseudoabs(data = sod4, x = "x", y = "y", n = nrow(sod3), 
                         method = c("geo_const", width = "5000"),#, env = crp_prd), 
                         rlayer = rly1, sp_name = "Phytophthora_ramorum")

psu2 <- part_random(data = psu1[,2:4], pr_ab = "pr_ab", method = c(method = "kfold",
  folds =4))

psu2$species <- "Phytophthora_ramorum"
# Combine the presence and pseudo-absence
sod_all <- rbind(sod4, psu2)

sod_all <- sod_all[,1:4]
# Extract environmental data
sod_all <- sod_all %>% sdm_extract(data = ., x = "x", y = "y", env_layer = env_rd, 
  filter_na = T)

bgpt1$species <- "Phytophthora_ramorum"
bgpt1 <- bgpt1[,1:4]
bgpt1 <- bgpt1 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = env_rd,
  filter_na = T)

sodnoa <- complete.cases(sod_all)
bgnoa <- complete.cases(bgpt1)
sod_all <- sod_all[sodnoa,]
bgpt1 <- bgpt1[bgnoa,]

###

set.seed(42)
tmax2 <- tune_max(data = sod_all, response = "pr_ab", predictors = names(env_rd),
  #predictors_f = c("NLCD Land Cover Class")                 ,
  background = bgpt1, partition = ".part", grid = expand.grid(regmult = seq(0.1, 
  3, 0.5), classes = c("l", "lq", "lqhpt")), thr = c("max_sens_spec", 
  "equal_sens_spec", "max_sorensen"), metric = "TSS", clamp = TRUE, 
  pred_type = "cloglog")

glm1 <- fit_glm(data = sod_all, response = "pr_ab", predictors = names(env_rd),
  partition = ".part", thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
  select_pred = F, poly =0, inter_order = 0#,predictors_f = c("NLCD Land Cover Class"))
)

sod_ens <- fit_ensemble(models = list(tmax2, glm1), ens_method = "meanw", 
                        thr = c("equal_sens_spec", "max_sens_spec"),
                        thr_model = "max_sens_spec", metric = "TSS")


# Predict
pred_area <- vect("Z:/Late_blight/helpful_shapes_rasters/sodextent.gpkg")
pred_area <- project(pred_area, y = crs(rast(cropped_predictors)))
# Threshold maximizes TSS
prd1 <- sdm_predict(models = sod_ens, pred = env_rd, thr = "max_sens_spec",
                    predict_area = env_rd)
