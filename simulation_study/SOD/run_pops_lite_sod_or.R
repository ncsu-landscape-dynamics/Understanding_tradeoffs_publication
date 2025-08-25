library(PoPS)
library(terra)
library(doParallel)
library(foreach)

config_file <- "/Volumes/cmjone25/Late_Blight/Manuscript_1_Data/simulation/SOD/inputs/config_SOD_OR.rds"

pops_lite(config_file = config_file, number_of_cores = NULL,
          number_of_iterations = NULL, new_dirs_path = NULL,
          random_seed = NULL)
