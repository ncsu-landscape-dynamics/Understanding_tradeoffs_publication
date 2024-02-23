library(PoPS)
library(terra)
library(doParallel)
# library(plyr)

# Note that the google drive might have a different letter
ffIn <- "Z:/Late_blight/"
ffout <- "Z:/Late_blight/modelresults/"


# This has all the time steps for a calibration/validation duration
infected_years_file = paste0(ffIn, "infections/gpkg/years/2009/inf_years_file_2009.gpkg")#"infections/rasters/2009/usa_2009_infections.tif")
number_of_observations = 326
prior_number_of_observations = 0
prior_means = c(0, 0, 0, 0, 0, 0)
prior_cov_matrix = matrix(0, 6, 6)
params_to_estimate = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
number_of_generations = 7
generation_size = 1000
pest_host_table = paste0(ffIn, "hosts/pest_host_table_lb.csv")
competency_table = paste0(ffIn, "hosts/competency_table_original.csv")
# This is the starting infection for a calibration/validation perio
infected_file_list = paste0(ffIn, "infections/gpkg/filelist/2009/inf_filelist_2009.gpkg")
# Make the host files without NAs
host_file_list = paste0(ffIn, "hosts/host_2009_nona.tif")
# Make total populations files without NAs
total_populations_file = paste0(ffIn, "all_populations/tot_pop_2009_nona.tif")
temp = TRUE
temperature_coefficient_file = paste0(ffIn, "Weather_data/weather_us_2009.tif")
precip = FALSE
precipitation_coefficient_file = ""
model_type = "SI"
latency_period = 0
time_step = "month"
season_month_start = 1
season_month_end = 12
start_date = "2009-01-01"
end_date = "2009-12-31"
use_survival_rates = FALSE
survival_rate_month = 3
survival_rate_day = 15
survival_rates_file = ""
use_lethal_temperature = FALSE
temperature_file = ""
lethal_temperature = 5
lethal_temperature_month = 1
mortality_frequency = "year"
mortality_frequency_n = 1
management = FALSE
treatment_dates = c("")
treatments_file = ""
treatment_method = "ratio"
natural_kernel_type = "cauchy"
anthropogenic_kernel_type = "cauchy"
natural_dir = "NONE"
natural_kappa = 0
anthropogenic_dir = "NONE"
anthropogenic_kappa = 0
pesticide_duration = c(0)
pesticide_efficacy = 1
mask = NULL
output_frequency = "year"
output_frequency_n = 1
movements_file = ""
use_movements = FALSE
start_exposed = FALSE
generate_stochasticity = TRUE
establishment_stochasticity = TRUE
movement_stochasticity = TRUE
dispersal_stochasticity = TRUE
establishment_probability = 0.5
dispersal_percentage = 0.99
quarantine_areas_file = ""
use_quarantine = FALSE
use_spreadrates = FALSE
use_overpopulation_movements = FALSE
overpopulation_percentage = 0
leaving_percentage = 0
leaving_scale_coefficient = 1
calibration_method = "ABC"
number_of_iterations = 1e+05
exposed_file_list = ""
verbose = TRUE
write_outputs = "None"
output_folder_path = ffout
network_filename = ""
network_movement = "walk"
success_metric = "mcc"
use_initial_condition_uncertainty = FALSE
use_host_uncertainty = FALSE
weather_type = "deterministic"
temperature_coefficient_sd_file = ""
precipitation_coefficient_sd_file = ""
dispersers_to_soils_percentage = 0
quarantine_directions = ""
multiple_random_seeds = FALSE
file_random_seeds = NULL
use_soils = FALSE
soil_starting_pest_file = ""
start_with_soil_populations = FALSE
county_level_infection_data = TRUE

#started at 20:35 on 4/19
# Started at 9 on 2/13
lb_cal = calibrate(infected_years_file,
                               number_of_observations,
                               prior_number_of_observations,
                               prior_means,
                               prior_cov_matrix,
                               params_to_estimate,
                               number_of_generations,
                               generation_size,
                               pest_host_table,
                               competency_table,
                               infected_file_list,
                               host_file_list,
                               total_populations_file,
                               temp,
                               temperature_coefficient_file,
                               precip,
                               precipitation_coefficient_file,
                               model_type,
                               latency_period,
                               time_step,
                               season_month_start,
                               season_month_end,
                               start_date,
                               end_date,
                               use_survival_rates,
                               survival_rate_month,
                               survival_rate_day,
                               survival_rates_file,
                               use_lethal_temperature,
                               temperature_file,
                               lethal_temperature,
                               lethal_temperature_month,
                               mortality_frequency,
                               mortality_frequency_n,
                               management,
                               treatment_dates,
                               treatments_file,
                               treatment_method,
                               natural_kernel_type,
                               anthropogenic_kernel_type,
                               natural_dir,
                               natural_kappa,
                               anthropogenic_dir,
                               anthropogenic_kappa,
                               pesticide_duration,
                               pesticide_efficacy,
                               mask,
                               output_frequency,
                               output_frequency_n,
                               movements_file,
                               use_movements,
                               start_exposed,
                               generate_stochasticity,
                               establishment_stochasticity,
                               movement_stochasticity,
                               dispersal_stochasticity,
                               establishment_probability,
                               dispersal_percentage,
                               quarantine_areas_file,
                               use_quarantine,
                               use_spreadrates,
                               use_overpopulation_movements,
                               overpopulation_percentage,
                               leaving_percentage,
                               leaving_scale_coefficient,
                               calibration_method,
                               number_of_iterations,
                               exposed_file_list,
                               verbose,
                               write_outputs,
                               output_folder_path,
                               network_filename,
                               network_movement,
                               success_metric,
                               use_initial_condition_uncertainty,
                               use_host_uncertainty,
                               weather_type,
                               temperature_coefficient_sd_file,
                               precipitation_coefficient_sd_file,
                               dispersers_to_soils_percentage,
                               quarantine_directions,
                               multiple_random_seeds,
                               file_random_seeds,
                               use_soils,
                               soil_starting_pest_file,
                               start_with_soil_populations,
                               county_level_infection_data)

