library(PoPS)
root <- "/Volumes/cmjone25/Late_Blight/Manuscript_1_Data/simulation/SOD"
inpath <- file.path(root, "inputs")

config <- list()

# Assign values to config
config$random_seed <- NULL
config$infected_file_list <- file.path(inpath, "infections", "mean_sd_end_inf_2021_eu1.tif")
config$host_file_list <- file.path(inpath, "hosts", "mean_sd_updated_tanoak.tif")
config$total_populations_file <- file.path(inpath, "hosts", "lemma_max100m.tif")
config$parameter_means <- t(read.table(file.path(inpath, "parameters", "sod_eu1_means.csv")))[1, ]
config$parameter_cov_matrix <- read.table(file.path(inpath, "parameters", "sod_eu1_cov_mat.csv"))
config$temp <- TRUE
config$temperature_coefficient_file <- file.path(inpath, "Weather", "weather_coefficients.tif")
config$precip <- FALSE
config$precipitation_coefficient_file <- ""
config$model_type <- "SI"
config$latency_period <- 0
config$time_step <- "week"
config$season_month_start <- 1
config$season_month_end <- 12
config$start_date <- "2022-01-01"
config$end_date <- "2024-12-31"
config$use_lethal_temperature <- FALSE
config$temperature_file <- ''
config$lethal_temperature <- -15
config$lethal_temperature_month <- 1
config$use_survival_rates <- FALSE
config$survival_rate_month <- 3
config$survival_rate_day <- 15
config$survival_rates_file <- ""
config$mortality_on <- FALSE
config$mortality_rate <- 0
config$mortality_time_lag <- 0
config$mortality_frequency <- "year"
config$mortality_frequency_n <- 1
config$management <- FALSE
config$treatment_dates <- c("2023-12-24")
config$treatments_file <- ""
config$treatment_method <- "ratio"
config$natural_kernel_type <- "cauchy"
config$anthropogenic_kernel_type <- "cauchy"
config$natural_dir <- "NONE"
config$anthropogenic_dir <- "NONE"
config$pesticide_duration <- c(0)
config$pesticide_efficacy <- 1.0
config$output_frequency <- "year"
config$output_frequency_n <- 1
config$movements_file <- ""
config$use_movements <- FALSE
config$start_exposed <- FALSE
config$generate_stochasticity <- TRUE
config$establishment_stochasticity <- TRUE
config$movement_stochasticity <- TRUE
config$dispersal_stochasticity <- TRUE
config$establishment_probability <- 0.5
config$dispersal_percentage <- 0.99
config$quarantine_areas_file <- ""
config$quarantine_directions <- ""
config$use_quarantine <- FALSE
config$use_spreadrates <- FALSE
config$use_overpopulation_movements <- FALSE
config$overpopulation_percentage <- 0.75
config$leaving_percentage <- 0.5
config$leaving_scale_coefficient <- 5
config$number_of_iterations <- 4000
config$number_of_cores <- 5
config$function_name <- "multirun"
config$failure <- NULL
config$exposed_file_list <- ""
config$mask <- NULL
config$write_outputs <- "all_simulations"
config$output_folder_path <- file.path(root, "outputs")
config$mortality_frequency <- "year"
config$mortality_frequency_n <- 1
config$network_filename <- ''
config$network_movement <- "walk"
config$use_initial_condition_uncertainty <- TRUE
config$use_host_uncertainty <- TRUE
config$weather_type <- "deterministic"
config$temperature_coefficient_sd_file <- ""
config$precipitation_coefficient_sd_file <- ""
config$dispersers_to_soils_percentage <- 0
config$multiple_random_seeds <- TRUE
config$file_random_seeds <- NULL
config$use_soils <- FALSE
config$soil_starting_pest_file <- ""
config$start_with_soil_populations <- FALSE
config$county_level_infection_data <- FALSE
config$pest_host_table <- file.path(inpath, "hosts", "pest_host_table_SOD_OR.csv")
config$competency_table <- file.path(inpath, "hosts", "competency_table_SOD_OR.csv")

config <- configuration(config)
saveRDS(config, file.path(inpath, "config_SOD_OR.rds"))
