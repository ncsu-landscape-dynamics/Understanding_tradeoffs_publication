library(terra)
root <- "/rs1/researchers/c/cmjone25/Late_blight/Manuscript_1_Data/simulation/LB"
output_path <- file.path(root, "outputs21/locs6x2/reg2")
raster_path <- file.path(root, "outputs21/locs6x2/reg2/rasts")
config_path <- file.path(root, "inputs")
inpath <- file.path(root, "input")


config <- readRDS(file.path(config_path, "config_LB_21_sim_locs6_lo_hpc.rds"))
#config$number_of_iterations <- 10
raster_template <- rast(config$host_file_list[[1]])[[1]]

filelist <- list.files(file.path(output_path), pattern = "pops_output*")

yearly_stats <- data.frame(number_infected = rep(0, config$number_of_iterations), infected_area = rep(0, config$number_of_iterations))
#summary_stats <- data.frame(number_infected_mean = rep(0, number_of_years),
#                            infected_area_mean = rep(0, number_of_years),
#                            number_infected_sd = rep(0, number_of_years),
#                            infected_area_sd = rep(0, number_of_years))

inf_indices <- seq(1, 2 * config$number_of_iterations, 2)
area_indices <- seq(2, 2 * config$number_of_iterations, 2)

rasts <- lapply(1:config$number_of_iterations, function(r){
  file <- readRDS(file.path(output_path, filelist[r]))
  values(raster_template) <- file$host_pools[[1]]$infected[[1]]
  raster_template
})
all_rasts <- rast(rasts)

x <- unlist(lapply(1:config$number_of_iterations, function(r){
  file <- readRDS(file.path(output_path, filelist[r]))
  yearly_stats[r, ] <- c(file$number_infected, file$area_infected)
}))

yearly_stats <- data.frame(number_infected = x[inf_indices], infected_area = x[area_indices])


all_means <- terra::mean(all_rasts)
sd_s <- terra::stdev(all_rasts)
which_median <- function(x) which.min(abs(x - terra::median(x, na.rm=T)))
median_run_index <- which_median(yearly_stats$number_infected)
min_run_index <- which.min(yearly_stats$number_infected)
max_run_index <- which.max(yearly_stats$number_infected)
median_run <- rasts[[median_run_index]]
min_run <- rasts[[min_run_index]]
max_run <- rasts[[max_run_index]]

writeRaster(median_run, file.path(raster_path, paste0("pops_median_Year_", 2021, ".tif")), overwrite = TRUE, gdal = c("COMPRESS=NONE"))
writeRaster(min_run, file.path(raster_path, paste0("pops_min_Year_", 2021, ".tif")), overwrite = TRUE, gdal = c("COMPRESS=NONE"))
writeRaster(max_run, file.path(raster_path, paste0("pops_max_Year_", 2021, ".tif")), overwrite = TRUE, gdal = c("COMPRESS=NONE"))
writeRaster(all_means, file.path(raster_path, paste0("pops_mean_Year_", 2021, ".tif")), overwrite = TRUE, gdal = c("COMPRESS=NONE"))
writeRaster(sd_s, file.path(raster_path, paste0("pops_sd_Year_", 2021, ".tif")), overwrite = TRUE, gdal = c("COMPRESS=NONE"))
write.csv(yearly_stats, (file.path(raster_path, paste0("yearly_stats", 2021, ".csv"))))
