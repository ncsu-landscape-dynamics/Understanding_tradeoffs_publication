library(flexsdm)
library(terra)
library(data.table)
library(parallel)
library(spatialEco)
library(tidyverse)

path <- "Z:/pops_pesthostuse/pops.sdm/Data/"

setwd("Z:/Late_blight/SDM/SOD/")
domain <- "USA"

extent <- "Z:/Late_blight/helpful_shapes_rasters/sodextent.gpkg"
extent <- vect(extent)
qboun <- vect("Z:/Late_blight/helpful_shapes_rasters/SOD_Quarantine_Boundary.gpkg")
res <- 100

# N of cells not NA in 26710 100 x 100 m
nnacl <- 620530 
# Infections
sod1 <- vect("Z:/Late_blight/helpful_shapes_rasters/SOD_all.gpkg")

# Predictors
cropped_predictorlis <- list.files("flexsdm_results/1_Inputs/2_Predictors/1_Current/cropped/transformed/", full.names = T)
subset_files <- cropped_predictorlis[c(1,7,20,22,24,25,27,30,33,36)]
cropped_predictors <- rast(lapply(subset_files, rast))

sod1 <- project(sod1, y = crs(rast(cropped_predictors))) # The y is from later..
sod1 <- sod1[sod1$YEAR %in% 2008:2021,]
sod1$x <- crds(sod1)[,1]
sod1$y <- crds(sod1)[,2]
sod2 <- as.data.frame(sod1[,c(1,6,7)])

source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/sdm_helpers.R")
source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/get_envi_chunked.R") # get_topo_global # nolint
source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/raster_base.R") # base_raster at correct resolution and extent # nolint
source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/part_sblock.R") # nolint
source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/sample_background.R") # nolint

res <- fix_resolution(res, domain)
base <- crop_base_raster(domain, res, path, extent)
domain <- "USA"

crp_prd <- cropped_predictors

sod3 <- sod2[,c(2:3)]
sod3$pr_ab <- 1
sod3$id <- 1:nrow(sod3)

set.seed(42)


sod3b <- sod3[,1:3]
sod3b$id <- paste0(sod3b$x,"_",sod3b$y)

# 100
# Random partitioning for calib/valid
prt_sod100 <- part_random(data = sod3b, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                 folds = 4))
# Observations/infections
sod5 <- prt_sod100

# Background
bgpt100 <- sample_background(
  data = sod5, x= "x", y ="y", n = nrow(sod5)*2, method = "random", 
  rlayer = crp_prd#, calibarea = lbbufmcp#, 
  #calibarea = 
) %>% bind_rows()

prt_sodb <- part_random(data = bgpt100, pr_ab = "pr_ab", method = c(method="kfold", 
                                                                  folds = 4))

bgpt100 <- prt_sodb

# Pseudo absences, 2 steps
psu1 <- sample_pseudoabs(data = sod5, x = "x", y = "y", n = nrow(sod5), 
                         method = c("geo_const", width = "5000"),#, env = crp_prd), 
                         rlayer = crp_prd, sp_name = "Phytophthora_ramorum")

psu2 <- part_random(data = psu1[,2:4], pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds =4))

# Combine the presence and pseudo-absence
sod_all <- rbind(sod5[,-1], psu2)

sod_all <- sod_all[,1:4]
# Extract environmental data
sod_all <- sod_all %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                   filter_na = T)

bgpt100 <- bgpt100[,1:4]
bgpt100 <- bgpt100 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd,
                               filter_na = T)

sodnoa <- complete.cases(sod_all)
bgnoa <- complete.cases(bgpt100)
sod_all <- sod_all[sodnoa,]
bgpt100 <- bgpt100[bgnoa,]

# Maxent
tmax2 <- tune_max(data = sod_all, response = "pr_ab", predictors = names(crp_prd),
                  #predictors_f = c("NLCD Land Cover Class")                 ,
                  background = bgpt100, partition = ".part", grid = expand.grid(regmult = seq(0.1, 
                                                                                            3, 0.5), classes = c("l", "lq", "lqhpt")), thr = c("max_sens_spec", 
                                                                                                                                             "equal_sens_spec", "max_sorensen"), metric = "TSS", clamp = TRUE, 
                  pred_type = "cloglog")

# GLM
glm1 <- fit_glm(data = sod_all, response = "pr_ab", predictors = names(crp_prd),
                partition = ".part", thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
                select_pred = F, poly =0, inter_order = 0#,predictors_f = c("NLCD Land Cover Class"))
)

maxe_pr <- sdm_predict(tmax2, pred = env_rd, thr = "max_sens_spec", predict_area = env_rd)

glm_pr <- sdm_predict(glm1, pred= env_rd, thr = "max_sens_spec", predict_area = env_rd)

sod_ens <- fit_ensemble(models = list(tmax2, glm1), ens_method = "meanw", 
                        thr = c("equal_sens_spec", "max_sens_spec"),
                        thr_model = "max_sens_spec", metric = "TSS")
saveRDS(sod_ens, "Z:/Late_blight/SDM/SOD/sod100_redo_filterthenreduc.RDS")

prd1 <- sdm_predict(models = sod_ens, pred = crp_prd, thr = "max_sens_spec",
                    predict_area = crp_prd)

# Same ext and res for predictions...
rastemp <- rast("Z:/Late_blight/SOD_OR/sod_PoPSpred_resampledtomatchSDMWRM26710.tif")
sodsdmmm <- project(rast(prd1), y = crs(rastemp))
sodsdmmm <- resample(sodsdmmm$meanw$meanw, rastemp)
sodsdmmm <- crop(sodsdmmm, y = rastemp, mask=T)
# Crop to quarantine area
qboun <- project(qboun, y = crs(rastemp))
sodsdmq <- crop(sodsdmmm, y = qboun, mask=T)
writeRaster(sodsdmmm, "Z:/Late_blight/SDM/SOD/ensrast100_redo_nofilter.tif", overwrite=T)
writeRaster(sodsdmq, "Z:/Late_blight/SDM/SOD/ensrast100_redo_filterthenreduc.tif", overwrite=T)


# 90%
sod90 <- sod5[sample(nrow(sod5), round(nrow(sod5) * 0.9)),]

prt_sod90 <- part_random(data = sod90, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                    folds = 4))

# Background
bgpt90 <- sample_background(
  data = prt_sod90, x= "x", y ="y", n = nrow(prt_sod90)*2, method = "random", 
  rlayer = crp_prd#, calibarea = lbbufmcp#, 
  #calibarea = 
) %>% bind_rows()

bgp90pt <- part_random(data = bgpt90, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                  folds = 4))

bgpt90 <- bgp90pt

# Pseudo absences
psu1 <- sample_pseudoabs(data = prt_sod90, x = "x", y = "y", n = nrow(prt_sod90)*2, 
                         method = c("geo_const", width = "5000"),#, env = crp_prd), 
                         rlayer = crp_prd, sp_name = "Phytophthora_ramorum")

psu2 <- part_random(data = psu1[,2:4], pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds =4))


# Combine the presence and pseudo-absence
sod_all90 <- rbind(prt_sod90[,-1], psu2)

sod_all90 <- sod_all90[,1:4]
# Extract environmental data
sod_all90 <- sod_all90 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                   filter_na = T)

bgpt90 <- bgpt90[,1:4]
bgpt90 <- bgpt90 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd,
                                   filter_na = T)

sodnoa <- complete.cases(sod_all90)
bgnoa <- complete.cases(bgpt90)
sod_all90 <- sod_all90[sodnoa,]
bgpt90 <- bgpt90[bgnoa,]

tmax90 <- tune_max(data = sod_all90, response = "pr_ab", predictors = names(crp_prd),
                  #predictors_f = c("NLCD Land Cover Class")                 ,
                  background = bgpt90, partition = ".part", grid = expand.grid(regmult = seq(0.1, 
                  3, 0.5), classes = c("l", "lq", "lqhpt")), thr = c("max_sens_spec", 
                 "equal_sens_spec", "max_sorensen"), metric = "TSS", clamp = TRUE, 
                  pred_type = "cloglog")

glm90 <- fit_glm(data = sod_all90, response = "pr_ab", predictors = names(crp_prd),
                partition = ".part", thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
                select_pred = F, poly =0, inter_order = 0#,predictors_f = c("NLCD Land Cover Class"))
)

sod_ens90 <- fit_ensemble(models = list(tmax90, glm90), ens_method = "meanw", 
                        thr = c("equal_sens_spec", "max_sens_spec"),
                        thr_model = "max_sens_spec", metric = "TSS")

saveRDS(sod_ens90, "Z:/Late_blight/SDM/SOD/ensmodel90pctdata_filterthenreduc.RDS")

prd90 <- sdm_predict(models = sod_ens90, pred = crp_prd, thr = "max_sens_spec",
                    predict_area = crp_prd)

sodsdmmm <- project(rast(prd90), y = crs(rastemp))
sodsdmmm <- resample(sodsdmmm$meanw$meanw, rastemp)
sodsdmmm <- crop(sodsdmmm, y = rastemp, mask=T)
qboun <- project(qboun, y = crs(rastemp))
sodsdmq <- crop(sodsdmmm, y = qboun, mask=T)
writeRaster(sodsdmmm, "Z:/Late_blight/SDM/SOD/ensrast90pctdata_filterthenreduc.tif", overwrite=T)
writeRaster(sodsdmq, "Z:/Late_blight/SDM/SOD/ensrast90pctqurdata_filterthenreduc.tif", overwrite=T)


# 80%
sod80 <- sod5[sample(nrow(sod5), round(nrow(sod5)*0.8)),]

prt_sod80 <- part_random(data = sod80, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds = 4))

# Background
bgpt80 <- sample_background(
  data = prt_sod80, x= "x", y ="y", n = nrow(prt_sod80)*2, method = "random", 
  rlayer = crp_prd#, calibarea = lbbufmcp#, 
  #calibarea = 
) %>% bind_rows()

bgp80pt <- part_random(data = bgpt80, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                  folds = 4))

bgpt80 <- bgp80pt

# Pseudo absences
psu1 <- sample_pseudoabs(data = prt_sod80, x = "x", y = "y", n = nrow(prt_sod80)*2, 
                         method = c("geo_const", width = "5000"),#, env = crp_prd), 
                         rlayer = crp_prd, sp_name = "Phytophthora_ramorum")

psu2 <- part_random(data = psu1[,2:4], pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds =4))


# Combine the presence and pseudo-absence
sod_all80 <- rbind(prt_sod80[,-1], psu2)

sod_all80 <- sod_all80[,1:4]
# Extract environmental data
sod_all80 <- sod_all80 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                       filter_na = T)

bgpt80 <- bgpt80[,1:4]
bgpt80 <- bgpt80 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd,
                                 filter_na = T)

sodnoa <- complete.cases(sod_all80)
bgnoa <- complete.cases(bgpt80)
sod_all80 <- sod_all80[sodnoa,]
bgpt80 <- bgpt80[bgnoa,]

tmax80 <- tune_max(data = sod_all80, response = "pr_ab", predictors = names(crp_prd),
                   #predictors_f = c("NLCD Land Cover Class")                 ,
                   background = bgpt80, partition = ".part", grid = expand.grid(regmult = seq(0.1, 
                                                                                              3, 0.5), classes = c("l", "lq", "lqhpt")), thr = c("max_sens_spec", 
                                                                                                                                                 "equal_sens_spec", "max_sorensen"), metric = "TSS", clamp = TRUE, 
                   pred_type = "cloglog")

glm80 <- fit_glm(data = sod_all80, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
                 select_pred = F, poly =0, inter_order = 0#,predictors_f = c("NLCD Land Cover Class"))
)

sod_ens80 <- fit_ensemble(models = list(tmax80, glm80), ens_method = "meanw", 
                          thr = c("equal_sens_spec", "max_sens_spec"),
                          thr_model = "max_sens_spec", metric = "TSS")

saveRDS(sod_ens80, "Z:/Late_blight/SDM/SOD/ensmodel80pctdata_filterthenreduc.RDS")

# Threshold maximizes TSS
prd80 <- sdm_predict(models = sod_ens80, pred = crp_prd, thr = "max_sens_spec",
                     predict_area = crp_prd)

sodsdmmm <- project(rast(prd80), y = crs(rastemp))
sodsdmmm <- resample(sodsdmmm$meanw$meanw, rastemp)
sodsdmmm <- crop(sodsdmmm, y = rastemp, mask=T)
qboun <- project(qboun, y = crs(rastemp))
sodsdmq <- crop(sodsdmmm, y = qboun, mask=T)
writeRaster(sodsdmmm, "Z:/Late_blight/SDM/SOD/ensrast80pctdata_filterthenreduc.tif", overwrite=T)
writeRaster(sodsdmq, "Z:/Late_blight/SDM/SOD/ensrast80pctqurdata_filterthenreduc.tif", overwrite=T)


# 70%
sod70 <- sod5[sample(nrow(sod5), round(nrow(sod5)*0.7)),]
prt_sod70 <- part_random(data = sod70, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds = 4))

# Background
bgpt70 <- sample_background(
  data = prt_sod70, x= "x", y ="y", n = nrow(prt_sod70)*2, method = "random", 
  rlayer = crp_prd#, calibarea = lbbufmcp#, 
  #calibarea = 
) %>% bind_rows()

bgp70pt <- part_random(data = bgpt70, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                  folds = 4))

bgpt70 <- bgp70pt

# Pseudo absences
psu1 <- sample_pseudoabs(data = prt_sod70, x = "x", y = "y", n = nrow(prt_sod70)*2, 
                         method = c("geo_const", width = "5000"),#, env = crp_prd), 
                         rlayer = crp_prd, sp_name = "Phytophthora_ramorum")

psu2 <- part_random(data = psu1[,2:4], pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds =4))


# Combine the presence and pseudo-absence
sod_all70 <- rbind(prt_sod70[,-1], psu2)

sod_all70 <- sod_all70[,1:4]
# Extract environmental data
sod_all70 <- sod_all70 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                       filter_na = T)

bgpt70 <- bgpt70[,1:4]
bgpt70 <- bgpt70 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd,
                                 filter_na = T)

sodnoa <- complete.cases(sod_all70)
bgnoa <- complete.cases(bgpt70)
sod_all70 <- sod_all70[sodnoa,]
bgpt70 <- bgpt70[bgnoa,]

tmax70 <- tune_max(data = sod_all70, response = "pr_ab", predictors = names(crp_prd),
                   #predictors_f = c("NLCD Land Cover Class")                 ,
                   background = bgpt70, partition = ".part", grid = expand.grid(regmult = seq(0.1, 
                                                                                              3, 0.5), classes = c("l", "lq", "lqhpt")), thr = c("max_sens_spec", 
                                                                                                                                                 "equal_sens_spec", "max_sorensen"), metric = "TSS", clamp = TRUE, 
                   pred_type = "cloglog")

glm70 <- fit_glm(data = sod_all70, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
                 select_pred = F, poly =0, inter_order = 0#,predictors_f = c("NLCD Land Cover Class"))
)

sod_ens70 <- fit_ensemble(models = list(tmax70, glm70), ens_method = "meanw", 
                          thr = c("equal_sens_spec", "max_sens_spec"),
                          thr_model = "max_sens_spec", metric = "TSS")

saveRDS(sod_ens70, "Z:/Late_blight/SDM/SOD/ensmodel70pctdata_filterthenreduc.RDS")

prd70 <- sdm_predict(sod_ens70, , pred = crp_prd, thr = "max_sens_spec",
                     predict_area = crp_prd)

sodsdmmm <- project(rast(prd70), y = crs(rastemp))
sodsdmmm <- resample(sodsdmmm$meanw$meanw, rastemp)
sodsdmmm <- crop(sodsdmmm, y = rastemp, mask=T)
sodsdmq <- crop(sodsdmmm, y = qboun, mask=T)
writeRaster(sodsdmmm, "Z:/Late_blight/SDM/SOD/ensrast70pctdata_filterthenreduc.tif", overwrite=T)
writeRaster(sodsdmq, "Z:/Late_blight/SDM/SOD/ensrast70pctqurdata_filterthenreduc.tif", overwrite=T)

# 60%
sod60 <- sod5[sample(nrow(sod5), round(nrow(sod5)*0.6)),]
prt_sod60 <- part_random(data = sod60, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds = 4))

# Background
bgpt60 <- sample_background(
  data = prt_sod60, x= "x", y ="y", n = nrow(prt_sod60)*2, method = "random", 
  rlayer = crp_prd#, calibarea = lbbufmcp#, 
  #calibarea = 
) %>% bind_rows()

bgp60pt <- part_random(data = bgpt60, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                  folds = 4))

bgpt60 <- bgp60pt

# Pseudo absences
psu1 <- sample_pseudoabs(data = prt_sod60, x = "x", y = "y", n = nrow(prt_sod60)*2, 
                         method = c("geo_const", width = "5000"),#, env = crp_prd), 
                         rlayer = crp_prd, sp_name = "Phytophthora_ramorum")

psu2 <- part_random(data = psu1[,2:4], pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds =4))

# Combine the presence and pseudo-absence
sod_all60 <- rbind(prt_sod60[,-1], psu2)

sod_all60 <- sod_all60[,1:4]
# Extract environmental data
sod_all60 <- sod_all60 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                       filter_na = T)

bgpt60 <- bgpt60[,1:4]
bgpt60 <- bgpt60 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd,
                                 filter_na = T)

sodnoa <- complete.cases(sod_all60)
bgnoa <- complete.cases(bgpt60)
sod_all60 <- sod_all60[sodnoa,]
bgpt60 <- bgpt60[bgnoa,]

tmax60 <- tune_max(data = sod_all60, response = "pr_ab", predictors = names(crp_prd),
                   #predictors_f = c("NLCD Land Cover Class")                 ,
                   background = bgpt60, partition = ".part", grid = expand.grid(regmult = seq(0.1, 
                                                                                              3, 0.5), classes = c("l", "lq", "lqhpt")), thr = c("max_sens_spec", 
                                                                                                                                                 "equal_sens_spec", "max_sorensen"), metric = "TSS", clamp = TRUE, 
                   pred_type = "cloglog")

glm60 <- fit_glm(data = sod_all60, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
                 select_pred = F, poly =0, inter_order = 0#,predictors_f = c("NLCD Land Cover Class"))
)

sod_ens60 <- fit_ensemble(models = list(tmax60, glm60), ens_method = "meanw", 
                          thr = c("equal_sens_spec", "max_sens_spec"),
                          thr_model = "max_sens_spec", metric = "TSS")

saveRDS(sod_ens60, "Z:/Late_blight/SDM/SOD/ensmodel60pctdata_filterthenreduc.RDS")

prd60 <- sdm_predict(models = sod_ens60, pred = crp_prd, thr = "max_sens_spec",
                     predict_area = crp_prd)

sodsdmmm <- project(rast(prd60), y = crs(rastemp))
sodsdmmm <- resample(sodsdmmm$meanw$meanw, rastemp)
sodsdmmm <- crop(sodsdmmm, y = rastemp, mask=T)
sodsdmq <- crop(sodsdmmm, y = qboun, mask=T)
writeRaster(sodsdmmm, "Z:/Late_blight/SDM/SOD/ensrast60pctdata_filterthenreduc.tif", overwrite=T)
writeRaster(sodsdmq, "Z:/Late_blight/SDM/SOD/ensrast60pctqurdata_filterthenreduc.tif", overwrite=T)


# 50%
sod50 <- sod5[sample(nrow(sod5), round(nrow(sod5)*0.5)),]
prt_sod50 <- part_random(data = sod50, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds = 4))

# Background
bgpt50 <- sample_background(
  data = prt_sod50, x= "x", y ="y", n = nrow(prt_sod50)*2, method = "random", 
  rlayer = crp_prd#, calibarea = lbbufmcp#, 
  #calibarea = 
) %>% bind_rows()

bgp50pt <- part_random(data = bgpt50, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                  folds = 4))

bgpt50 <- bgp50pt

psu1 <- sample_pseudoabs(data = prt_sod50, x = "x", y = "y", n = nrow(prt_sod50)*2, 
                         method = c("geo_const", width = "5000"),#, env = crp_prd), 
                         rlayer = crp_prd, sp_name = "Phytophthora_ramorum")

psu2 <- part_random(data = psu1[,2:4], pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds =4))

# Combine the presence and pseudo-absence
sod_all50 <- rbind(prt_sod50[,-1], psu2)

sod_all50 <- sod_all50[,1:4]
# Extract environmental data
sod_all50 <- sod_all50 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                       filter_na = T)

bgpt50 <- bgpt50[,1:4]
bgpt50 <- bgpt50 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd,
                                 filter_na = T)

sodnoa <- complete.cases(sod_all50)
bgnoa <- complete.cases(bgpt50)
sod_all50 <- sod_all50[sodnoa,]
bgpt50 <- bgpt50[bgnoa,]

tmax50 <- tune_max(data = sod_all50, response = "pr_ab", predictors = names(crp_prd),
                   #predictors_f = c("NLCD Land Cover Class")                 ,
                   background = bgpt50, partition = ".part", grid = expand.grid(regmult = seq(0.1, 
                                                                                              3, 0.5), classes = c("l", "lq", "lqhpt")), thr = c("max_sens_spec", 
                                                                                                                                                 "equal_sens_spec", "max_sorensen"), metric = "TSS", clamp = TRUE, 
                   pred_type = "cloglog")

glm50 <- fit_glm(data = sod_all50, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
                 select_pred = F, poly =0, inter_order = 0#,predictors_f = c("NLCD Land Cover Class"))
)

sod_ens50 <- fit_ensemble(models = list(tmax50, glm50), ens_method = "meanw", 
                          thr = c("equal_sens_spec", "max_sens_spec"),
                          thr_model = "max_sens_spec", metric = "TSS")

saveRDS(sod_ens50, "Z:/Late_blight/SDM/SOD/ensmodel50pctdata_filterthenreduc.RDS")

prd50 <- sdm_predict(models = sod_ens50, pred = crp_prd, thr = "max_sens_spec",
                     predict_area = crp_prd)

sodsdmmm <- project(rast(prd50), y = crs(rastemp))
sodsdmmm <- resample(sodsdmmm$meanw$meanw, rastemp)
sodsdmmm <- crop(sodsdmmm, y = rastemp, mask=T)
sodsdmq <- crop(sodsdmmm, y = qboun, mask=T)
writeRaster(sodsdmmm, "Z:/Late_blight/SDM/SOD/ensrast50pctdata_filterthenreduc.tif", overwrite=T)
writeRaster(sodsdmq, "Z:/Late_blight/SDM/SOD/ensrast50pctqurdata_filterthenreduc.tif", overwrite=T)

# 40%
# For filtering then doing the reduction. The prior line has to ignored/commentedout
sod40 <- sod5[sample(nrow(sod5), round(nrow(sod5)*0.4)),]
prt_sod40 <- part_random(data = sod40, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds = 4))

# Background
bgpt40 <- sample_background(
  data = prt_sod40, x= "x", y ="y", n = nrow(prt_sod40)*2, method = "random", 
  rlayer = crp_prd#, calibarea = lbbufmcp#, 
  #calibarea = 
) %>% bind_rows()

bgp40pt <- part_random(data = bgpt40, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                  folds = 4))

bgpt40 <- bgp40pt

# Pseudo absences
psu1 <- sample_pseudoabs(data = prt_sod40, x = "x", y = "y", n = nrow(prt_sod40)*2, 
                         method = c("geo_const", width = "4000"),#, env = crp_prd), 
                         rlayer = crp_prd, sp_name = "Phytophthora_ramorum")

psu2 <- part_random(data = psu1[,2:4], pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds =4))

# Combine the presence and pseudo-absence
sod_all40 <- rbind(prt_sod40[,-1], psu2)

sod_all40 <- sod_all40[,1:4]
# Extract environmental data
sod_all40 <- sod_all40 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                       filter_na = T)

bgpt40 <- bgpt40[,1:4]
bgpt40 <- bgpt40 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd,
                                 filter_na = T)

sodnoa <- complete.cases(sod_all40)
bgnoa <- complete.cases(bgpt40)
sod_all40 <- sod_all40[sodnoa,]
bgpt40 <- bgpt40[bgnoa,]

tmax40 <- tune_max(data = sod_all40, response = "pr_ab", predictors = names(crp_prd),
                   #predictors_f = c("NLCD Land Cover Class")                 ,
                   background = bgpt40, partition = ".part", grid = expand.grid(regmult = seq(0.1, 
                                                                                              3, 0.5), classes = c("l", "lq", "lqhpt")), thr = c("max_sens_spec", 
                                                                                                                                                 "equal_sens_spec", "max_sorensen"), metric = "TSS", clamp = TRUE, 
                   pred_type = "cloglog")

glm40 <- fit_glm(data = sod_all40, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
                 select_pred = F, poly =0, inter_order = 0#,predictors_f = c("NLCD Land Cover Class"))
)

sod_ens40 <- fit_ensemble(models = list(tmax40, glm40), ens_method = "meanw", 
                          thr = c("equal_sens_spec", "max_sens_spec"),
                          thr_model = "max_sens_spec", metric = "TSS")

saveRDS(sod_ens40, "Z:/Late_blight/SDM/SOD/ensmodel40pctdata_filterthenreduc.RDS")

prd40 <- sdm_predict(models = sod_ens40, pred = crp_prd, thr = "max_sens_spec",
                     predict_area = crp_prd)

sodsdmmm <- project(rast(prd40), y = crs(rastemp))
sodsdmmm <- resample(sodsdmmm$meanw$meanw, rastemp)
sodsdmmm <- crop(sodsdmmm, y = rastemp, mask=T)
sodsdmq <- crop(sodsdmmm, y = qboun, mask=T)
writeRaster(sodsdmmm, "Z:/Late_blight/SDM/SOD/ensrast40pctdata_filterthenreduc.tif", overwrite=T)
writeRaster(sodsdmq, "Z:/Late_blight/SDM/SOD/ensrast40pctqurdata_filterthenreduc.tif", overwrite=T)

# 30%
sod30 <- sod5[sample(nrow(sod5), round(nrow(sod5)*0.3)),]
prt_sod30 <- part_random(data = sod30, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds = 4))

bgpt30 <- sample_background(
  data = prt_sod30, x= "x", y ="y", n = nrow(prt_sod30)*2, method = "random", 
  rlayer = crp_prd#, calibarea = lbbufmcp#, 
  #calibarea = 
) %>% bind_rows()

bgp30pt <- part_random(data = bgpt30, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                  folds = 4))

bgpt30 <- bgp30pt

psu1 <- sample_pseudoabs(data = prt_sod30, x = "x", y = "y", n = nrow(prt_sod30)*2, 
                         method = c("geo_const", width = "3000"),#, env = crp_prd), 
                         rlayer = crp_prd, sp_name = "Phytophthora_ramorum")

psu2 <- part_random(data = psu1[,2:4], pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds =4))

# Combine the presence and pseudo-absence
sod_all30 <- rbind(prt_sod30[,-1], psu2)

sod_all30 <- sod_all30[,1:4]
# Extract environmental data
sod_all30 <- sod_all30 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                       filter_na = T)

bgpt30 <- bgpt30[,1:4]
bgpt30 <- bgpt30 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd,
                                 filter_na = T)

sodnoa <- complete.cases(sod_all30)
bgnoa <- complete.cases(bgpt30)
sod_all30 <- sod_all30[sodnoa,]
bgpt30 <- bgpt30[bgnoa,]

tmax30 <- tune_max(data = sod_all30, response = "pr_ab", predictors = names(crp_prd),
                   #predictors_f = c("NLCD Land Cover Class")                 ,
                   background = bgpt30, partition = ".part", grid = expand.grid(regmult = seq(0.1, 
                                                                                              3, 0.5), classes = c("l", "lq", "lqhpt")), thr = c("max_sens_spec", 
                                                                                                                                                 "equal_sens_spec", "max_sorensen"), metric = "TSS", clamp = TRUE, 
                   pred_type = "cloglog")

glm30 <- fit_glm(data = sod_all30, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
                 select_pred = F, poly =0, inter_order = 0#,predictors_f = c("NLCD Land Cover Class"))
)

sod_ens30 <- fit_ensemble(models = list(tmax30, glm30), ens_method = "meanw", 
                          thr = c("equal_sens_spec", "max_sens_spec"),
                          thr_model = "max_sens_spec", metric = "TSS")

saveRDS(sod_ens30, "Z:/Late_blight/SDM/SOD/ensmodel30pctdata_filterthenreduc.RDS")

prd30 <- sdm_predict(models = sod_ens30, pred = crp_prd, thr = "max_sens_spec",
                     predict_area = crp_prd)

sodsdmmm <- project(rast(prd30), y = crs(rastemp))
sodsdmmm <- resample(sodsdmmm$meanw$meanw, rastemp)
sodsdmmm <- crop(sodsdmmm, y = rastemp, mask=T)
sodsdmq <- crop(sodsdmmm, y = qboun, mask=T)
writeRaster(sodsdmmm, "Z:/Late_blight/SDM/SOD/ensrast30pctdata_filterthenreduc.tif", overwrite=T)
writeRaster(sodsdmq, "Z:/Late_blight/SDM/SOD/ensrast30pctqurdata_filterthenreduc.tif", overwrite=T)

sod30binr <- sodsdmmm
sod30binr[sod30binr > 0.395] <- 1
sod30binr[sod30binr <= 0.395] <- 0
sod30brdf <- as.data.frame(sod30binr, xy =T)
ggplot(data = sod30brdf) + geom_tile(aes(x=x,y=y, fill = meanw)) + theme_void() +
  scale_fill_gradient(low = "antiquewhite", high = "red") + coord_equal()

sod30qbinr <- sodsdmq
sod30qbinr[sod30qbinr > 0.395] <- 1
sod30qbinr[sod30qbinr <= 0.395] <- 0
sod30qbrdf <- as.data.frame(sod30qbinr, xy =T)
ggplot(data = sod30qbrdf) + geom_tile(aes(x=x,y=y, fill = meanw)) + theme_void() +
  scale_fill_gradient(low = "antiquewhite", high = "red") + coord_equal()

# 20%
sod20 <- sod20 %>% occfilt_env(data = ., x = "x", y = "y", id = "id", 
                 nbins = c(4), env_layer = crp_prd) %>% left_join(sod4, by = c("id",
                 "x","y"))

# For filtering then doing the reduction. The prior line has to ignored/commentedout
#sod20 <- sod5[sample(nrow(sod5), round(nrow(sod5*0.2))),]
prt_sod20 <- part_random(data = sod20, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds = 4))

bgpt20 <- sample_background(
  data = prt_sod20, x= "x", y ="y", n = nrow(prt_sod20)*2, method = "random", 
  rlayer = crp_prd#, calibarea = lbbufmcp#, 
  #calibarea = 
) %>% bind_rows()


psu1 <- sample_pseudoabs(data = prt_sod20, x = "x", y = "y", n = nrow(prt_sod20)*2, 
                         method = c("geo_const", width = "2000"),#, env = crp_prd), 
                         rlayer = crp_prd, sp_name = "Phytophthora_ramorum")

psu2 <- part_random(data = psu1[,2:4], pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds =4))


bgp20pt <- part_random(data = bgpt20, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                  folds = 4))

bgpt20 <- bgp20pt

# Combine the presence and pseudo-absence
sod_all20 <- rbind(prt_sod20[,-1], psu2)

sod_all20 <- sod_all20[,1:4]
# Extract environmental data
sod_all20 <- sod_all20 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                       filter_na = T)

bgpt20 <- bgpt20[,1:4]
bgpt20 <- bgpt20 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd,
                                 filter_na = T)

sodnoa <- complete.cases(sod_all20)
bgnoa <- complete.cases(bgpt20)
sod_all20 <- sod_all20[sodnoa,]
bgpt20 <- bgpt20[bgnoa,]

tmax20 <- tune_max(data = sod_all20, response = "pr_ab", predictors = names(crp_prd),
                   #predictors_f = c("NLCD Land Cover Class")                 ,
                   background = bgpt20, partition = ".part", grid = expand.grid(regmult = seq(0.1, 
                                                                                              3, 0.5), classes = c("l", "lq", "lqhpt")), thr = c("max_sens_spec", 
                                                                                                                                                 "equal_sens_spec", "max_sorensen"), metric = "TSS", clamp = TRUE, 
                   pred_type = "cloglog")

#readr::write_rds(tmax2, "Z:/Late_blight/SDM/SOD/maxent.RDS")
#tmax2 <- readRDS("Z:/Late_blight/SDM/SOD/maxent.RDS")

glm20 <- fit_glm(data = sod_all20, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
                 select_pred = F, poly =0, inter_order = 0#,predictors_f = c("NLCD Land Cover Class"))
)

sod_ens20 <- fit_ensemble(models = list(tmax20, glm20), ens_method = "meanw", 
                          thr = c("equal_sens_spec", "max_sens_spec"),
                          thr_model = "max_sens_spec", metric = "TSS")

saveRDS(sod_ens20, "Z:/Late_blight/SDM/SOD/ensmodel20pctdata_filterthenreduc.RDS")

# Predict
# Threshold maximizes TSS
prd20 <- sdm_predict(models = sod_ens20, pred = crp_prd, thr = "max_sens_spec",
                     predict_area = crp_prd)

sodsdmmm <- project(rast(prd20), y = crs(rastemp))
sodsdmmm <- resample(sodsdmmm$meanw$meanw, rastemp)
sodsdmmm <- crop(sodsdmmm, y = rastemp, mask=T)
sodsdmq <- crop(sodsdmmm, y = qboun, mask=T)
writeRaster(sodsdmmm, "Z:/Late_blight/SDM/SOD/ensrast20pctdata_filterthenreduc.tif", overwrite=T)
writeRaster(sodsdmq, "Z:/Late_blight/SDM/SOD/ensrast20pctqurdata_filterthenreduc.tif", overwrite=T)

# 10%
sod10 <- sod5[sample(nrow(sod5), round(nrow(sod5)*0.1)),]
prt_sod10 <- part_random(data = sod10, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds = 4))

# Background
bgpt10 <- sample_background(
  data = prt_sod10, x= "x", y ="y", n = nrow(prt_sod10)*2, method = "random", 
  rlayer = crp_prd#, calibarea = lbbufmcp#, 
  #calibarea = 
) %>% bind_rows()

bgp10pt <- part_random(data = bgpt10, pr_ab = "pr_ab", method = c(method = "kfold",
                                                                  folds = 4))

bgpt10 <- bgp10pt

# Pseudo absence
psu1 <- sample_pseudoabs(data = prt_sod10, x = "x", y = "y", n = nrow(prt_sod10)*2, 
                         method = c("geo_const", width = "1000"),#, env = crp_prd), 
                         rlayer = crp_prd, sp_name = "Phytophthora_ramorum")

psu2 <- part_random(data = psu1[,2:4], pr_ab = "pr_ab", method = c(method = "kfold",
                                                                   folds =4))

# Combine the presence and pseudo-absence
sod_all10 <- rbind(prt_sod10[,-1], psu2)

sod_all10 <- sod_all10[,1:4]
# Extract environmental data
sod_all10 <- sod_all10 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                       filter_na = T)

bgpt10 <- bgpt10[,1:4]
bgpt10 <- bgpt10 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd,
                                 filter_na = T)

sodnoa <- complete.cases(sod_all10)
bgnoa <- complete.cases(bgpt10)
sod_all10 <- sod_all10[sodnoa,]
bgpt10 <- bgpt10[bgnoa,]

tmax10 <- tune_max(data = sod_all10, response = "pr_ab", predictors = names(crp_prd),
                   #predictors_f = c("NLCD Land Cover Class")                 ,
                   background = bgpt10, partition = ".part", grid = expand.grid(regmult = seq(0.1, 
                                                                                              3, 0.5), classes = c("l", "lq", "lqhpt")), thr = c("max_sens_spec", 
                                                                                                                                                 "equal_sens_spec", "max_sorensen"), metric = "TSS", clamp = TRUE, 
                   pred_type = "cloglog")

glm10 <- fit_glm(data = sod_all10, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
                 select_pred = F, poly =0, inter_order = 0#,predictors_f = c("NLCD Land Cover Class"))
)

sod_ens10 <- fit_ensemble(models = list(tmax10, glm10), ens_method = "meanw", 
                          thr = c("equal_sens_spec", "max_sens_spec"),
                          thr_model = "max_sens_spec", metric = "TSS")

saveRDS(sod_ens10, "Z:/Late_blight/SDM/SOD/ensmodel10pctdata_filterthenreduc.RDS")

prd10 <- sdm_predict(models = sod_ens10, pred = crp_prd, thr = "max_sens_spec",
                     predict_area = crp_prd)

sodsdmmm <- project(rast(prd10), y = crs(rastemp))
sodsdmmm <- resample(sodsdmmm$meanw$meanw, rastemp)
sodsdmmm <- crop(sodsdmmm, y = rastemp, mask=T)
sodsdmq <- crop(sodsdmmm, y = qboun, mask=T)
writeRaster(sodsdmmm, "Z:/Late_blight/SDM/SOD/ensrast10pctdata_filterthenreduc.tif", overwrite=T)
writeRaster(sodsdmq, "Z:/Late_blight/SDM/SOD/ensrast10pctqurdata_filterthenreduc.tif", overwrite=T)


# Change to binary after running the no filter version
# Change to binary after running the no filter version
# Change to binary after running the no filter version
nnacl <- 620530 
nofillis <- list.files("Z:/Late_blight/SDM/SOD/", full.names = T, pattern = "nofilter")
nofillis <- gtools::mixedsort(nofillis[c(10,seq(11,28, by = 2))])

jenlis <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/SOD/", full.names = T, 
                                       pattern = "sjenks"))

reclasf <- function(x) {
  # SDM rast
  rsdm <- rast(nofillis[x])
  # Jenks
  jenk1 <- read.csv(jenlis[x])
  # df to use for reclass
  jenk1$y <- jenk1$x[c(2:11,11)]
  jenk1$dec <- 1:11
  jenk1 <- round(jenk1[1:10,], 3)
  
  rclr1 <- classify(rsdm, jenk1, right =T)
  
  dfrq <- as.data.frame(freq(rclr1))
  dfrq$cusum <- cumsum(dfrq$count)
  dfrq$pct <- dfrq$cusum/nnacl
  # df for reference
  dfrq <- bind_cols(jenk1, dfrq)
  # value from 50 to get low&high
  lohithrsh <- dfrq$y[dfrq$dec == 5]
  
  rcl2 <- rsdm
  rcl2[rcl2 > lohithrsh] <- 1
  rcl2[rcl2 <= lohithrsh] <- 0
  frqdf <- as.data.frame(freq(rcl2))
  write.csv(dfrq, paste0("Z:/Late_blight/SDM/SOD/recldf_unf_",x, ".csv"), row.names=F)
  write.csv(frqdf, paste0("Z:/Late_blight/SDM/SOD/lohict_unf_",x, ".csv"), row.names=F)
  writeRaster(rcl2, paste0("Z:/Late_blight/SDM/SOD/relcls_unf_",x,".tif"), overwrite=T)
}

lapply(1:10, reclasf)

# Change to binary after running the filterthenreduce version
# Change to binary after running the filterthenreduce version
# Change to binary after running the filterthenreduce version
filredlis <- list.files("Z:/Late_blight/SDM/SOD/", full.names = T, pattern = "filterthenreduc")
filredlis <- gtools::mixedsort(filredlis[c(10,seq(11,28, by = 2))])

jenlis <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/SOD/", full.names = T, 
                                       pattern = "sjenkfr"))

reclasftd <- function(x) {
  # SDM rast
  rsdm <- rast(nofillis[x])
  # Jenks
  jenk1 <- read.csv(jenlis[x])
  # df to use for reclass
  jenk1$y <- jenk1$x[c(2:11,11)]
  jenk1$dec <- 1:11
  jenk1 <- round(jenk1[1:10,], 3)
  
  rclr1 <- classify(rsdm, jenk1, right =T)
  
  dfrq <- as.data.frame(freq(rclr1))
  dfrq$cusum <- cumsum(dfrq$count)
  dfrq$pct <- dfrq$cusum/nnacl
  # df for reference
  dfrq <- bind_cols(jenk1, dfrq)
  # value from 50 to get low&high
  lohithrsh <- dfrq$y[dfrq$dec == 5]
  
  rcl2 <- rsdm
  rcl2[rcl2 > lohithrsh] <- 1
  rcl2[rcl2 <= lohithrsh] <- 0
  frqdf <- as.data.frame(freq(rcl2))
  write.csv(dfrq, paste0("Z:/Late_blight/SDM/SOD/recldf_ftd_",x, ".csv"), row.names=F)
  write.csv(frqdf, paste0("Z:/Late_blight/SDM/SOD/lohict_ftd_",x, ".csv"), row.names=F)
  writeRaster(rcl2, paste0("Z:/Late_blight/SDM/SOD/relcls_ftd_",x,".tif"), overwrite=T)
}

lapply(1:10, reclasftd)



