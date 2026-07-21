#' Run the species distribution model: currently runs for Maxent only. Will add in other models later.
'The files in Z:/Late_blight/SDM/lateblight/flexsdm are all of eastern CONUS.
The files in Z:/Late_blight/SDM/lateblight/revise/flexsdm are a geographic subset.'

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

'theta_s 	saturated soil water content 	m3/m3
theta_r 	residual soil water content 	m3/m3'

#' 1. Load the required packages & set path
library(flexsdm)
library(terra)
library(data.table)
library(parallel)
library(spatialEco)
library(tidyverse)

set.seed(42)

path <- "Z:/pops_pesthostuse/pops.sdm/Data/"

setwd("Z:/Late_blight/SDM/lateblight/reduceddata/")
domain <- "USA"
extent <- "Z:/Late_blight/helpful_shapes_rasters/lb_extent.gpkg"

extent2 <- "Z:/Late_blight/helpful_shapes_rasters/lb_test_ext.gpkg"

# Area for validation (predict)
val_ext <- sf::st_difference(sf::st_as_sf(vect(extent)), sf::st_as_sf(crex1))
val_ext <- vect(val_ext)

rastemp <- rast("Z:/Late_blight/helpful_shapes_rasters/econus_mask_3857.tif")

res <- 1000

# Infections
lb2 <- vect("Z:/Late_blight/helpful_shapes_rasters/blightptsmovedforSDM.gpkg")
lbdf <- as.data.frame(lb2[,c(1,16:17)])
lbdf$id <- paste0(lbdf$x,"_",lbdf$y)
lbdf <- lbdf[!duplicated(lbdf$id),]
lb2sel <- lb2[lb2$Rec_ID %in% lbdf$Rec_ID,]


source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/sdm_helpers.R")
source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/get_envi_chunked.R") # get_topo_global # nolint
source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/raster_base.R") # base_raster at correct resolution and extent # nolint
source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/part_sblock.R") # nolint
source("C:/Users/japolo/Documents/code/blight_related/pops.sdm-main/sample_background.R") # nolint

# Predictors
#cropped_predictors <- c(cropped_predictors, rhm2)
cropped_predictors <- list.files("Z:/Late_blight/SDM/lateblight/flexsdm_results/1_Inputs/2_Predictors/1_Current/cropped/transformed/",
                                 pattern = "zscor", full.names = T)
cropped_predictors <- cropped_predictors[c(14,15,5,7,20,21,29,32,34,36:38)]
crp_prd <- rast(cropped_predictors)

# 5. Occurrence filtering
lb2sel <- as.data.frame(lb2sel[,c(1,16,17)])
lb2sel$pr_ab <- 1
lb2sel$id <- 1:nrow(lb2sel)
lb2sdf <- lb2sel

# 100
# The 100% should start here:
lb5 <- lb2sdf

# Spatial block partitioning
prt_lb2s <- part_sblock(env_layer = crp_prd, data = lb5, x = "x", 
                        y = "y", pr_ab = "pr_ab", n_part = 6, min_res_mult = 4, max_res_mult = 25)

# Get best block arrangement
bll100 <- get_block(crp_prd, prt_lb2s$grid)

lb5 <- prt_lb2s$part                        

# Background points
bgpt100 <- lapply(1:6, \(x) {
  sample_background(
    data = lb5, x= "x", y ="y", n = (sum(lb5$.part == x)) , 
    method = c("thickening", width = 6000), rlayer = crp_prd
  )
}) %>% bind_rows()

bgpt100 <- sdm_extract(bgpt100, x = "x", y = "y", bll100)

# Pseudo absences
lb100psa <- lapply(1:4, \(x) {
  sample_pseudoabs(
    data = lb5, x= "x", y ="y", n = sum(lb5$.part == x)*5, method = c("geo_const",
                            width = "3000"), rlayer = crp_prd[[1]]
  )
}) %>% bind_rows()

lb100psa <- sdm_extract(lb100psa, x="x", y="y", bll100)

#  Combine presences and pseudoabsences
lb100pa <- bind_rows(lb5, lb100psa)

## Extract environmental data
lb100pa <- lb100pa %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                 filter_na = T)

bgpt100 <- bgpt100 %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                 filter_na = T)

# Maxent
max100 <- tune_max(data = lb100pa, response = "pr_ab", predictors = names(crp_prd),
                  background = bgpt100, partition = ".part", thr = c("max_sens_spec",
                                                                    "equal_sens_spec"), metric = "TSS", grid = expand.grid(
                                                                      regmult = seq(0.1, 0.3, 0.5), classes=c("l","lq","lqhpt")))
# GLM
glm100 <- fit_glm(data = lb100pa, response = "pr_ab", predictors = names(crp_prd),
        #c("bio4","bio5","bio13","bio15","gdd","evi","sand_0_5","theta_s_0_5",
        #  "elevation","slope","dist2roads"),
                 partition = ".part", thr = c("max_sens_spec", "equal_sens_spec"), 
                 select_pred = FALSE, poly = 0, inter_order = 0)

# Ensemble
ens100 <- fit_ensemble(models = list(max100, glm100), ens_method = "meanw", 
                      thr = "max_sens_spec", thr_model = "max_sens_spec", 
                      metric = "TSS")

# Prediction
ensp100 <- sdm_predict(models = ens100, pred = crp_prd, thr = "max_sens_spec", 
                      predict_area = crp_prd)

saveRDS(ens100, "Z:/Late_blight/SDM/lateblight/ensmodellb100pct_filterthenreduc.RDS")

# Get same res and extent for prediction raster
lbsdmmm <- project(ensp100$meanw$meanw, y = crs(rastemp))
lbsdmmm <- resample(lbsdmmm, rastemp)
lbsdmmm <- crop(lbsdmmm, y = rastemp, mask=T)

writeRaster(lbsdmmm, "Z:/Late_blight/SDM/ensrastlb100_filterthenreduc.tif", overwrite=T)



# 90%
lb90t <- lb5[sample(nrow(lb5), round(nrow(lb5)*0.9)),]
lb90d <- lb90t

# Background
bgpt90 <- bgpt100[sample(nrow(bgpt100), round(nrow(bgpt100*0.9))),]

# Pseudo absences
lb90psa <- lb100psa[sample(nrow(lb100psa), round(nrow(lb100psa)*0.9)),1:4]

# Combine pseudo absences
lb90pa <- bind_rows(lb90d, lb90psa)

## Extract environmental data
lb90pa <- lb90pa %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
       filter_na = T)

# Max ent
max90 <- tune_max(data = lb90pa, response = "pr_ab", predictors = names(crp_prd),
            background = bgpt90, partition = ".part", thr = c("max_sens_spec",
            "equal_sens_spec"), metric = "TSS", grid = expand.grid(
              regmult = seq(0.1, 0.3, 0.5), classes=c("l","lq","lqhpt")))

# GLM
glm90 <- fit_glm(data = lb90pa, response = "pr_ab", predictors = names(crp_prd),
          partition = ".part", thr = c("max_sens_spec",
          "equal_sens_spec"), select_pred = FALSE, inter_order = 0)

ens90 <- fit_ensemble(models = list(max90, glm90), ens_method = "meanw", 
                      thr = "max_sens_spec", thr_model = "max_sens_spec", 
                      metric = "TSS")

ensp90 <- sdm_predict(models = ens90, pred = crp_prd, thr = "max_sens_spec", 
                      predict_area = crp_prd)

saveRDS(ens90, "Z:/Late_blight/SDM/lateblight/ensmodellb90_filterthenreduc.RDS")


lbsdmmm <- project(ensp90$meanw$meanw, y = crs(rastemp))
lbsdmmm <- resample(lbsdmmm, rastemp)
lbsdmmm <- crop(lbsdmmm, y = rastemp, mask=T)
writeRaster(lbsdmmm, "Z:/Late_blight/SDM/ensrastlb90_filterthenreduc.tif", overwrite=T)


# 80%
lb80t <- lb5[sample(nrow(lb5), round(nrow(lb5)*0.8)),]
lb80d <- lb80t

# Background
bgpt80 <- bgpt100[sample(nrow(bgpt100), round(nrow(bgpt100)*0.8)),]

# Pseudo absences
lb80psa <- lb100psa[sample(nrow(lb100psa), round(nrow(lb100psa)*0.8)),1:4]

lb80pa <- bind_rows(lb80d, lb80psa)

## Extract environmental data
lb80pa <- lb80pa %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                 filter_na = T)

# Maxent
max80 <- tune_max(data = lb80pa, response = "pr_ab", predictors = names(crp_prd),
                  background = bgpt80, partition = ".part", thr = c("max_sens_spec",
                                                                    "equal_sens_spec"), metric = "TSS", grid = expand.grid(
                                                                      regmult = seq(0.1, 0.3, 0.5), classes=c("l","lq","lqhpt")))

# GLM
glm80 <- fit_glm(data = lb90pa, response = "pr_ab", predictors = names(crp_prd),
       partition = ".part", thr = c("max_sens_spec",
          "equal_sens_spec"), select_pred = FALSE, inter_order = 0)

ens80 <- fit_ensemble(models = list(max80, glm90), ens_method = "meanw", 
                      thr = "max_sens_spec", thr_model = "max_sens_spec", 
                      metric = "TSS")

ensp80 <- sdm_predict(models = ens80, pred = crp_prd, thr = "max_sens_spec", 
                      predict_area = crp_prd)

saveRDS(ens80, "Z:/Late_blight/SDM/lateblight/ensmodellb80pctdata_filterthenreduc.RDS")


lbsdmmm <- project(ensp80$meanw$meanw, y = crs(rastemp))
lbsdmmm <- resample(lbsdmmm, rastemp)
lbsdmmm <- crop(lbsdmmm, y = rastemp, mask=T)
writeRaster(lbsdmmm, "Z:/Late_blight/SDM/ensrastlb80_filterthenreduc.tif", overwrite=T)
write.csv(lb80d, "Z:/Late_blight/SDM/lbftr80.csv", row.name =F)

# 70
lb70d <- lb5[sample(nrow(lb5), round(nrow(lb5)*0.7)),]

# Background
bgpt70 <- bgpt100[sample(nrow(bgpt100), round(nrow(bgpt100)*0.7)),]

lb70psa <- lb100pa[sample(nrow(lb100psa), round(nrow(lb100psa)*0.7)),1:4]

lb70pa <- bind_rows(lb70d, lb70psa)

## Extract environmental data
lb70pa <- lb70pa %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                 filter_na = T)

max70 <- tune_max(data = lb70pa, response = "pr_ab", predictors = names(crp_prd),
                  background = bgpt70, partition = ".part", thr = c("max_sens_spec",
                                                                    "equal_sens_spec"), metric = "TSS", grid = expand.grid(
                                                                      regmult = seq(0.1, 0.3, 0.5), classes=c("l","lq","lqhpt")))

glm70 <- fit_glm(data = lb90pa, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec",
                                              "equal_sens_spec"), select_pred = FALSE, inter_order = 0)

ens70 <- fit_ensemble(models = list(max70, glm90), ens_method = "meanw", 
                      thr = "max_sens_spec", thr_model = "max_sens_spec", 
                      metric = "TSS")

ensp70 <- sdm_predict(models = ens70, pred = crp_prd, thr = "max_sens_spec", 
                      predict_area = crp_prd)

saveRDS(ens70, "Z:/Late_blight/SDM/lateblight/ensmodellb70_filterthenreduc.RDS")

lbsdmmm <- project(ensp70$meanw, y = crs(rastemp))
lbsdmmm <- resample(lbsdmmm, rastemp)
lbsdmmm <- crop(lbsdmmm, y = rastemp, mask=T)
writeRaster(lbsdmmm, "Z:/Late_blight/SDM/lateblight/ensrastlb70_filterthenreduc.tif", overwrite=T)
write.csv(lb70d, "Z:/Late_blight/SDM/lbftr70.csv", row.names = F)


#60%
lb60d <- lb5[sample(nrow(lb5), round(nrow(lb5)*0.6)),]

# Background
bgpt60 <- bgpt100[sample(nrow(bgpt100), round(nrow(bgpt100)*0.6)),]

# Pseudo absences
lb60psa <- lb100psa[sample(nrow(lb100psa), round(nrow(lb100psa)*0.6)),1:4]

# Combine
lb60pa <- bind_rows(lb60d, lb60psa)

## Extract environmental data
lb60pa <- lb60pa %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                 filter_na = T)

# Maxent
max60 <- tune_max(data = lb60pa, response = "pr_ab", predictors = names(crp_prd),
                  background = bgpt60, partition = ".part", thr = c("max_sens_spec",
                                                                    "equal_sens_spec"), metric = "TSS", grid = expand.grid(
                                                                      regmult = seq(0.1, 0.3, 0.5), classes=c("l","lq","lqhpt")))

# GLM
glm60 <- fit_glm(data = lb90pa, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec",
                                              "equal_sens_spec"), select_pred = FALSE, inter_order = 0)

ens60 <- fit_ensemble(models = list(max60, glm90), ens_method = "meanw", 
                      thr = "max_sens_spec", thr_model = "max_sens_spec", 
                      metric = "TSS")

ensp60 <- sdm_predict(models = ens60, pred = crp_prd, thr = "max_sens_spec", 
                      predict_area = crp_prd)

saveRDS(ens60, "Z:/Late_blight/SDM/lateblight/ensmodellb60_filterthenreduc.RDS")

lbsdmmm <- project(ensp60$meanw, y = crs(rastemp))
lbsdmmm <- resample(lbsdmmm, rastemp)
lbsdmmm <- crop(lbsdmmm, y = rastemp, mask=T)
writeRaster(lbsdmmm, "Z:/Late_blight/SDM/ensrastlb60_filterthenreduc.tif", overwrite=T)
write.csv(lb60d,"Z:/Late_blight/SDM/lbftr60.csv", row.names = F)

# 50%
lb50d <- lb5[sample(nrow(lb5), round(nrow(lb5)*0.5)),]

# Background
bgpt50 <- bgpt100[sample(nrow(bgpt100), round(nrow(bgpt100)*0.5)),]

# Pseudo absences
lb50psa <- lb100psa[sample(nrow(lb100psa), round(nrow(lb100psa)*0.5)),]

lb50pa <- bind_rows(lb50d, lb50psa)

## Extract environmental data
lb50pa <- lb50pa %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                 filter_na = T)

max50 <- tune_max(data = lb50pa, response = "pr_ab", predictors = names(crp_prd),
                  background = bgpt50, partition = ".part", thr = c("max_sens_spec",
                                                                    "equal_sens_spec"), metric = "TSS", grid = expand.grid(
                                                                      regmult = seq(0.1, 0.3, 0.5), classes=c("l","lq","lqhpt")))

glm50 <- fit_glm(data = lb90pa, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec",
                                              "equal_sens_spec"), select_pred = FALSE, inter_order = 0)

ens50 <- fit_ensemble(models = list(max50, glm90), ens_method = "meanw", 
                      thr = "max_sens_spec", thr_model = "max_sens_spec", 
                      metric = "TSS")

ensp50 <- sdm_predict(models = ens50, pred = crp_prd, thr = "max_sens_spec", 
                      predict_area = crp_prd)

saveRDS(ens50, "Z:/Late_blight/SDM/lateblight/ensmodellb50_filterthanreduc.RDS")

lbsdmmm <- project(ensp50$meanw, y = crs(rastemp))
lbsdmmm <- resample(lbsdmmm, rastemp)
lbsdmmm <- crop(lbsdmmm, y = rastemp, mask=T)
writeRaster(lbsdmmm, "Z:/Late_blight/SDM/ensrastlb50_filterthenreduc.tif", overwrite=T)
write.csv(lb50d, "Z:/Late_blight/SDM/lbftr50.csv", row.names=F)

# 40
lb40d <- lb5[sample(nrow(lb5), round(nrow(lb5)*0.4)),]

# Background
bgpt40 <- bgpt100[sample(nrow(bgpt100), round(nrow(bgpt100)*0.4)),]

#Pseudo absences
lb40psa <- lb100psa[sample(nrow(lb100psa), round(nrow(lb100psa)*0.4)),]

lb40pa <- bind_rows(lb40d, lb40psa)

## Extract environmental data
lb40pa <- lb40pa %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                 filter_na = T)

max40 <- tune_max(data = lb40pa, response = "pr_ab", predictors = names(crp_prd),
                  background = bgpt40, partition = ".part", thr = c("max_sens_spec",
                                                                    "equal_sens_spec"), metric = "TSS", grid = expand.grid(
                                                                      regmult = seq(0.1, 0.3, 0.5), classes=c("l","lq","lqhpt")))

glm40 <- fit_glm(data = lb90pa, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec",
                                              "equal_sens_spec"), select_pred = FALSE, inter_order = 0)

ens40 <- fit_ensemble(models = list(max40, glm90), ens_method = "meanw", 
                      thr = "max_sens_spec", thr_model = "max_sens_spec", 
                      metric = "TSS")

ensp40 <- sdm_predict(models = ens40, pred = crp_prd, thr = "max_sens_spec", 
                      predict_area = crp_prd)

saveRDS(ens40, "Z:/Late_blight/SDM/lateblight/ensmodellb40_filterthenreduc.RDS")

lbsdmmm <- project(ensp40$meanw, y = crs(rastemp))
lbsdmmm <- resample(lbsdmmm, rastemp)
lbsdmmm <- crop(lbsdmmm, y = rastemp, mask=T)
writeRaster(lbsdmmm, "Z:/Late_blight/SDM/ensrastlb40p_filterthenreduc.tif", overwrite=T)
write.csv(lb40d, "Z:/Late_blight/SDM/lbftr40.csv", row.names=F)

# 30
lb30d <- lb5[sample(nrow(lb5), round(nrow(lb5)*0.3)),]

# Background
bgpt30 <- bgpt100[sample(nrow(bgpt100), round(nrow(bgpt100)*0.3)),]

# Pseudo absences
lb30psa <- lb100psa[sample(nrow(lb100psa), round(nrow(lb100psa)*0.3)),]

# Combine
lb30pa <- bind_rows(lb30d, lb30psa)

## Extract environmental data
lb30pa <- lb30pa %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                 filter_na = T)

max30 <- tune_max(data = lb30pa, response = "pr_ab", predictors = names(crp_prd),
                  background = bgpt30, partition = ".part", thr = c("max_sens_spec",
                                                                    "equal_sens_spec"), metric = "TSS", grid = expand.grid(
                                                                      regmult = seq(0.1, 0.3, 0.5), classes=c("l","lq","lqhpt")))

glm30 <- fit_glm(data = lb90pa, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec",
                                              "equal_sens_spec"), select_pred = FALSE, inter_order = 0)

ens30 <- fit_ensemble(models = list(max30, glm90), ens_method = "meanw", 
                      thr = "max_sens_spec", thr_model = "max_sens_spec", 
                      metric = "TSS")

ensp30 <- sdm_predict(models = ens30, pred = crp_prd, thr = "max_sens_spec", 
                      predict_area = crp_prd)

saveRDS(ens30, "Z:/Late_blight/SDM/lateblight/ensmodellb30_filterthenreduc.RDS")

lbsdmmm <- project(ensp30$meanw, y = crs(rastemp))
lbsdmmm <- resample(lbsdmmm, rastemp)
lbsdmmm <- crop(lbsdmmm, y = rastemp, mask=T)
writeRaster(lbsdmmm, "Z:/Late_blight/SDM/ensrastlb30_filterthenreduc.tif", overwrite=T)
write.csv(lb30d, "Z:/Late_blight/SDM/lbftr30.csv")

#20
lb20d <- lb5[sample(nrow(lb5), round(nrow(lb5)*0.2)),]

# Background
bgpt20 <- bgpt100[sample(nrow(bgpt100), round(nrow(bgpt100)*0.2)),]


# Pseudo absences
lb20psa <- lb100psa[sample(nrow(lb100psa), round(nrow(lb100psa)*0.2)),]

#Combine
lb20pa <- bind_rows(lb20d, lb20psa)

## Extract environmental data
lb20pa <- lb20pa %>% sdm_extract(data = ., x = "x", y = "y", env_layer = crp_prd, 
                                 filter_na = T)

max20 <- tune_max(data = lb20pa, response = "pr_ab", predictors = names(crp_prd),
                  background = bgpt20, partition = ".part", thr = c("max_sens_spec",
                                                                    "equal_sens_spec"), metric = "TSS", grid = expand.grid(
                                                                      regmult = seq(0.1, 0.3, 0.5), classes=c("l","lq","lqhpt")))

glm20 <- fit_glm(data = lb90pa, response = "pr_ab", predictors = names(crp_prd),
                 partition = ".part", thr = c("max_sens_spec",
                                              "equal_sens_spec"), select_pred = FALSE, inter_order = 0)

ens20 <- fit_ensemble(models = list(max20, glm90), ens_method = "meanw", 
                      thr = "max_sens_spec", thr_model = "max_sens_spec", 
                      metric = "TSS")

ensp20 <- sdm_predict(models = ens20, pred = crp_prd, thr = "max_sens_spec", 
                      predict_area = crp_prd)

saveRDS(ens20, "Z:/Late_blight/SDM/lateblight/ensmodellb20_filterthenreduc.RDS")

lbsdmmm <- project(ensp20$meanw, y = crs(rastemp))
lbsdmmm <- resample(lbsdmmm, rastemp)
lbsdmmm <- crop(lbsdmmm, y = rastemp, mask=T)
writeRaster(lbsdmmm, "Z:/Late_blight/SDM/ensrastlb20_filterthenreduc.tif", overwrite=T)
write.csv(lb20d, "Z:/Late_blight/SDM/lbftr20.csv", row.names=F)


# Change to binary after running the no filter version
# Change to binary after running the no filter version
# Change to binary after running the no filter version
nnacl <- 2628987 
nofillis <- list.files("Z:/Late_blight/SDM/lateblight", full.names = T, 
            pattern = "nof.*tif")
nofillis <- gtools::mixedsort(nofillis[c(10:18)])

jenlis <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/lateblight/", full.names = T, 
                                       pattern = "jenks.*nofil"))

reclasf <- function(x) {
  # SDM rast
  rsdm <- rast(nofillis[x])[[1]]
  print(nlyr(rsdm))
  # Jenks
  jenk1 <- read.csv(jenlis[x])
  # df to use for reclass
  jenk1$y <- jenk1$x[c(2:11,11)]
  jenk1$dec <- 1:11
  jenk1 <- round(jenk1[1:10,], 3)
  jenk1$x[1] = 0
  
  rclr1 <- classify(rsdm, jenk1, right =T)
  
  dfrq <- as.data.frame(freq(rclr1))
  dfrq$cusum <- cumsum(dfrq$count)
  dfrq$pct <- dfrq$cusum/nnacl
  print("jenk1")
  print(dim(jenk1))
  print("dfrq")
  print(dim(dfrq))
  # df for reference
  dfrq <- bind_cols(jenk1, dfrq)
  # value from 50 to get low&high
  lohithrsh <- dfrq$y[dfrq$dec == 5]
  
  rcl2 <- rsdm
  rcl2[rcl2 > lohithrsh] <- 1
  rcl2[rcl2 <= lohithrsh] <- 0
  plot(rcl2)
  frqdf <- as.data.frame(freq(rcl2))
  yr = as.numeric(as.character(str_match(jenlis[x], "[0-9]+")))
  write.csv(dfrq, paste0("Z:/Late_blight/SDM/lateblight/recldf_unf_",yr, ".csv"), row.names=F)
  write.csv(frqdf, paste0("Z:/Late_blight/SDM/lateblight/lohict_unf_",yr, ".csv"), row.names=F)
  writeRaster(rcl2, paste0("Z:/Late_blight/SDM/lateblight/relcls_unf_",yr,".tif"), overwrite=T)
}

lapply(1:9, reclasf)

# Change to binary after running the noe version
# Change to binary after running the noe version
# Change to binary after running the noe version
filredlis <- list.files("Z:/Late_blight/SDM/lateblight/", full.names = T, 
                        pattern = "filterth.*tif")
filredlis <- gtools::mixedsort(filredlis)

jenlis2 <- gtools::mixedsort(list.files("Z:/Late_blight/SDM/lateblight/", full.names = T, 
                                       pattern = "filtd"))
jenlis <- c(jenlis2, jenlis[9])

reclasftd <- function(x) {
  # SDM rast
  rsdm <- rast(filredlis[x])[[1]]
  # Jenks
  jenk1 <- read.csv(jenlis[x])
  # df to use for reclass
  jenk1$y <- jenk1$x[c(2:11,11)]
  jenk1$dec <- 1:11
  jenk1 <- round(jenk1[1:10,], 3)
  jenk1$x[1] = 0
  
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
  yr = as.numeric(as.character(str_match(jenlis[x], "[0-9]+")))
  write.csv(dfrq, paste0("Z:/Late_blight/SDM/lateblight/recldf_ftd_",yr, ".csv"), row.names=F)
  write.csv(frqdf, paste0("Z:/Late_blight/SDM/lateblight/lohict_ftd_",yr, ".csv"), row.names=F)
  writeRaster(rcl2, paste0("Z:/Late_blight/SDM/lateblight/relcls_ftd_",yr,".tif"), overwrite=T)
}

lapply(1:9, reclasftd)






