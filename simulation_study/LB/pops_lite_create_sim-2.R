library(PoPS)
library(terra)
library(doParallel)
library(foreach)

config_file <- "Z:/Late_Blight/Manuscript_1_Data/simulation/LB/inputs/config_LB_21_sim_locs10_up.rds"

pops_lite(config_file = config_file, number_of_cores = 4,
          number_of_iterations = NULL, new_dirs_path = NULL,
          random_seed = NULL)
