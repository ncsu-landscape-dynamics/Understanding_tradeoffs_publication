#
# These come from https://github.com/cyborginhas/pops.sdm
#

#' Description: Checks and helpers for run_sdm.R

#' @description Function to check if the domain provided by the user is in the
#' correct format. If the domain is not in the correct format, the function
#' returns an error message.
#' @param domain The domain provided by the user. Must be "USA" or "Global".

check_domain <- function(domain) {
  if (domain != "USA" && domain != "Global") {
    stop("Domain must be 'USA' or 'Global'")
  }
}

#' @description Function to update the resolution based on the user input,
#' and the domain provided by the user. If the domain is "USA", the
#' resolution is set to can only be 30, 100, 250, 500, 1000; if the
#' domain is "Global", the resolution can only be 1000, 2500, 5000.
#' @param res The resolution provided by the user.
#' @param domain The domain provided by the user.

fix_resolution <- function(res, domain) {
  if (domain == "USA") {
    common_res_names <- c(30, 100, 250, 500, 1000)
    res <- common_res_names[which.min(abs(common_res_names - res))]
    print(paste0("Resolution set to: ", res, "m"))
  } else if (domain == "Global") {
    common_res_names <- c(1000, 2500, 5000)
    res <- common_res_names[which.min(abs(common_res_names - res))]
    print(paste0("Resolution set to: ", res, "m"))
  }
  return(res)
}


#' @description This function checks if the species name provided by the user is
#' in the correct format. If the species name is not in the correct format, the
#' function returns the species name in the correct format.
#' @param species The species name provided by the user.

format_species_name <- function(species) {
  # Replace " " with "_" in species name
  species <- gsub(" ", "_", species)
  # ensure the first letter is in uppercase and the rest are in lowercase
  species <- tolower(species)
  # make first letter uppercase
  species <- paste0(
    toupper(substr(species, 1, 1)),
    substr(species, 2, nchar(species))
  )
  return(species)
}

#' @description This function reads in predictors and crops them to the domain
#' of interest. The function returns a list of cropped predictors.
#' @param file The predictor file to be cropped.
#' @param pred_copies_path The directory path to copy the predictor files
#' will be copied to.
#' @param pred_cropped_path The directory path to the cropped predictor
#' files will be saved.
#' @param extent The base raster cropped to the study extent.

crop_predictor <- function(
    file,
    pred_copies_path,
    pred_cropped_path,
    extent) {
  # Check if the file exists, if it does not, copy, crop & zscore
  cropped_file <- paste0(pred_cropped_path, basename(file))
  cropped_file <- gsub(".tif", "_cropped.tif", cropped_file)
  zscore_file <- gsub("/cropped/", "/cropped/transformed/", cropped_file)
  zscore_file <- gsub(".tif", "_zscore.tif", zscore_file)

  if (!file.exists(cropped_file)) {
    print(paste0("Cropping ", basename(file)))
    s <- Sys.time()
    # Copy the file to the predictors directory
    file.copy(file, pred_copies_path)
    copy <- paste0(pred_copies_path, basename(file))
    r <- terra::rast(copy)
    # Replace the path with the new path
    new_filename <- paste0(pred_cropped_path, basename(copy))
    # Add cropped to the filename
    new_filename <- gsub(".tif", "_cropped.tif", new_filename)
    dtype <- terra::datatype(r)
    # Crop the raster to the study extent
    terra::crop(r, extent,
      filename = new_filename,
      wopt = list(
        gdal = "COMPRESS=ZSTD",
        datatype = dtype, overwrite = TRUE
      )
    )
    cropped_raster <- terra::rast(new_filename)
    # Throw an error if copy does not contain the correct path
    if (!grepl("1_Inputs/2_Predictors/1_Current/copies", copy)) {
      stop("Attemping to delete a file that is not in the correct directory.")
    }
    # Change file access permissions
    Sys.chmod(copy, mode = "0777")
    # Remove the file
    file.remove(copy)
    t <- Sys.time() - s
    print(paste0(
      "Time taken to crop ", basename(new_filename), ": ",
      t, attr(t, "units")
    ))
    if (dtype != "INT1U") {
      zscore_raster <- spatialEco::raster.transformation(
        cropped_raster,
        trans = "std"
      )
      terra::writeRaster(zscore_raster, zscore_file,
        overwrite = TRUE, gdal = "COMPRESS=ZSTD", datatype = "FLT4S"
      )
    } else {
      zscore_raster <- cropped_raster
    }
    terra::tmpFiles(orphan = TRUE, old = TRUE, remove = TRUE)
  } else {
    cropped_raster <- terra::rast(cropped_file)
    dtype <- terra::datatype(cropped_raster)
    if (dtype != "INT1U") {
      zscore_raster <- terra::rast(zscore_file)
    } else {
      zscore_raster <- cropped_raster
    }
  }
  return(zscore_raster)
}

#' @description Function to crop species occurrence data to a specified extent
#' and year.
#' @param occs A data table with the species occurrence data retrieved from
#' batch_get_pts.
#' @param path The path to the species occurrence data. Should point to the main Data folder.
#' @param extent The cropped base map raster.
#' @param year The year of interest. If the year is not provided, the function
#' will crop the species occurrence data to the extent of interest.
#' @return A cropped data table with the species occurrence data.
#' @export crop_species_occurrences

prep_occurrences <- function(occs, extent, year = NULL) {
  occs$lat <- as.numeric(occs$lat)
  occs$lon <- as.numeric(occs$lon)
  occs_pts <- terra::vect(occs,
    crs = terra::crs(extent),
    geom = c("lon", "lat")
  )
  occs_extent <- terra::crop(occs_pts, extent)
  occs_extent$lon <- terra::crds(occs_extent)[, 1]
  occs_extent$lat <- terra::crds(occs_extent)[, 2]
  occs_extent <- data.table::as.data.table(occs_extent)

  # Filter by year if provided
  if (!is.null(year)) {
    occs_pa$date <- as.numeric(format(as.Date(occs_pa$date,
      format = "%Y"
    ), "%Y"))
    occs_extent <- occs_extent[date >= year]
  }
  return(occs_extent)
}

#' @description Function to retrieve base raster for domain and resolution
#' and extent of interest.
#' @param domain The domain of interest. Options are "USA" or "Global".
#' @param res The resolution of the base raster to create.
#' @param path Path for the Data folder
#' @param extent Path to vector file of extent of interest to crop the base raster.
#' The cropped based raster will be saved locally in specified output path.
#' @note Will add option to specify resolutions beyond 30m, 100m, 250m, 500m, 1000m
#' in future versions.

crop_base_raster <- function(domain, res, path, extent) {
  cropped_base_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/3_Calibration_area/base_study_extent.tif") # no lint
  if (!file.exists(cropped_base_path)) {
    # Retrieve base raster
    if (domain == "USA") {
      common_res_names <- c(30, 100, 250, 500, 1000)
      base <- base_conus(path)
      # Pull the base raster at the specified resolution
      base <- base[[which(common_res_names == res)]]
    } else if (domain == "Global") {
      common_res_names <- c(1000, 2500, 5000)
      base <- base_global(path)
      # Pull the base raster at the specified resolution
      base <- base[[which(common_res_names == res)]]
    } else {
      stop("Domain must be 'USA' or 'Global'")
    }
    # Crop the base raster to the extent of interest
    extent <- terra::vect(extent)
    extent <- terra::project(extent, y = terra::crs(base))
    base <- terra::crop(base, extent)
    terra::writeRaster(base, cropped_base_path,
      overwrite = TRUE,
      gdal = "COMPRESS=ZSTD", datatype = "INT1U"
    )
  } else {
    base <- terra::rast(cropped_base_path)
  }
  return(base)
}

#' @description Function to export the cropped species occurrence
#' to the local output path. With option to filter NEON to use
#' as independent validation data. Data filtered so that no
#' more than one occurrence per lat/lon pair is included.
#' @param occs The cropped species occurrence data.

export_occurrences <- function(occs) {
  sciname <- unique(occs$sciname)
  sciname <- format_species_name(sciname)
  neon <- occs[db == "neon", ]
  occs <- occs[db != "neon", ]
  neon <- unique(neon[, .(x = lon, y = lat, pr_ab = as.integer(p_a))])
  neon$id <- 1:nrow(neon)
  occs <- unique(occs[, .(x = lon, y = lat, pr_ab = as.integer(p_a))])
  occs$id <- 1:nrow(occs)
  # Create outpath
  occs_dirpath <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/original/"
  )
  dir.create(occs_dirpath, showWarnings = FALSE, recursive = TRUE)
  write.csv(occs,
    file = paste0(occs_dirpath, sciname, "_train_occs.csv"),
    row.names = FALSE
  )
  write.csv(neon,
    file = paste0(occs_dirpath, sciname, "_test_occs.csv"),
    row.names = FALSE
  )
  return(occs)
}

#' @description Function to help determine the range of values to apply for
#' geographic filtering. The function returns the range of values to apply
#' for geographic filtering. with the flexsdm::occfilt_geo function.
#' @param domain The domain of interest. Options are "USA" or "Global".
#' If the domain is "USA", the range of values is set to 90m intervals;
#' if the domain is "Global", the range of values is set to 150m intervals.
#' @param res The resolution of the analysis.

get_geo_cellsizes <- function(domain, res) {
  if (domain == "USA") {
    range <- seq(res, 1000, by = 3 * res)
  } else {
    range <- seq(res, 5000, by = 150)
  }
  return(range)
}

#' @description Function to apply geographic filtering to the species occurrence
#' data. The function returns the filtered species occurrence data.
#' @param occs The species occurrence data.
#' @param extent The base raster cropped to the extent of interest.
#' @param d The distance to apply for geographic filtering.
#' @param species The species name.

geo_filter_occs <- function(occs, extent, d, species) {
  filt_geo <- flexsdm::occfilt_geo(
    data = occs,
    x = "x",
    y = "y",
    env_layer = extent,
    method = c("defined", d = d / 1000)
  )
  filt_geo$d <- d
  # Write out the filtered data
  occsfilt_path <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/filtered/"
  )
  dir.create(occsfilt_path, showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(filt_geo, paste0(
    occsfilt_path,
    species, "_occ_pa_geofilter_", d, "m2.csv"
  ),
  row.names = FALSE
  )
  return(filt_geo)
}

#' @description Function to remove the predictors that are categorical
#' and retain a list of the integer and continuous predictors.
#' @param predictor A predictor to be used in the analysis.

get_numerical_rasters <- function(predictor) {
  dt <- terra::datatype(predictor)
  if (dt != "INT1U") {
    predictor <- predictor
  } else {
    predictor <- NULL
  }
  return(predictor)
}


#' @description Function to keep the predictors that are categorical
#' @param predictor A predictor to be used in the analysis.

get_categorical_rasters <- function(predictor) {
  dt <- terra::datatype(predictor)
  if (dt == "INT1U") {
    predictor <- predictor
  } else {
    predictor <- NULL
  }
  return(predictor)
}

#' @description Function to partition the data using  environmenal
#' and spatial cross-validation, and the inputs from geo_filter_occs.
#' The function returns the partitioned data.
#' @param occs The filtered species occurrence data.
#' @param env_layer The list of predictors.
#' @param species The species name.
#' @param d The distance to apply for geographic filtering.

spatial_block_partition <- function(occs, env_layer, species, d) {
  s <- Sys.time()
  # Partition the data
  part_data <- part_sblock_upd(
    env_layer = env_layer,
    data = occs,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    n_part = 4,
    min_res_mult = 3,
    max_res_mult = 200,
    num_grids = 20,
    prop = 0.5,
    min_occ = 10
  )
  t <- Sys.time() - s
  # Write out the partitioned data
  # Create outpath
  part_path <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/partitioned/"
  )
  dir.create(part_path, showWarnings = FALSE, recursive = TRUE)

  # Write out the partitioned data
  names(part_data) <- paste0(names(part_data), "_filtgeo_", d, "m2")
  part_data[[2]]$d <- d
  # Best partition
  data.table::fwrite(part_data[[2]], paste0(
    part_path, species, "_best_part_filtgeo", d, "m2.csv"
  ))
  # Partition data
  part_data[[1]]$d <- d
  data.table::fwrite(part_data[[1]], paste0(
    part_path, species, "_part_data_filtgeo", d, "m2.csv"
  ))
  names(part_data[[3]]) <- paste0(names(part_data[[3]]), "_filtgeo_", d, "m2")
  # Partition grid
  terra::writeRaster(part_data[[3]], paste0(
    part_path, species, "_part_grid_filtgeo", d, "m2.tif"
  ), overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=ZSTD")

  t <- Sys.time() - s
  print(paste0(
    "Time taken to partition ", length(occs$x),
    " points: ", t, attr(t, "units")
  ))
  return(part_data)
}

#' @description A function to create randomly sampled background points data
#' and write the data to the local output path.
#' @param partitioned_data A list of partitioned data occurrences, and the
#' associated spatial block grid.
#' @param species The species name.
#' associated grid.
#' @param res The resolution of the base raster.

create_random_bg_pts <- function(partitioned_data, species, res) {
  occs_p <- partitioned_data[[1]]
  grid_env <- partitioned_data[[2]]
  d <- unique(occs_p$d)
  
  if (res == 30) {
    factor <- floor(terra::res(grid_env)[1] / terra::res(base)[1] / 5)
    grid_env <- terra::disagg(grid_env, factor)
  } else {
    grid_env <- terra::resample(grid_env, base, method = "near")
  }

  bg_pts <- flexsdm::sample_background(
    data = occs_p,
    x = "x",
    y = "y",
    n = length(occs_p$x),
    method = "random",
    rlayer = grid_env
  )

  # Extract grid values
  part <- terra::extract(grid_env, bg_pts[, c("x", "y")], ID = FALSE)
  part <- part[[1]]
  bg_pts$.part <- part

  # Combine the background points with the occurrences
  occs_p <- data.table::rbindlist(list(occs_p, bg_pts),
    use.names = TRUE,
    fill = TRUE
  )
  occs_p$d <- d

  # Write out the background points
  bg_path <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/background/random/"
  )
  dir.create(bg_path, showWarnings = FALSE, recursive = TRUE)

  data.table::fwrite(occs_p, paste0(
    bg_path, species, "_bg_pts_filtgeo_", d, "m2.csv"
  ), row.names = FALSE)
  return(occs_p)
}

#' @description A function to create thickened background points data
#' and write the data to the local output path.
#' @param partitioned_data A list of partitioned data occurrences, and the
#' associated spatial block grid.
#' @param base The base raster cropped to the extent of interest.
#' @param res The resolution of the base raster.
#' @param species The species name.
#' associated grid.

create_thickened_bg_pts <- function(
    partitioned_data,
    base, res, species) {
  s <- Sys.time()
  occs_p <- partitioned_data[[1]]
  d <- unique(occs_p$d)
  occs_p <- occs_p[, .(pr_ab, x, y, .part)]
  grid_env <- partitioned_data[[2]]

  # Disaggregate the spatial block grid

  if (res == 30) {
    factor <- floor(terra::res(grid_env)[1] / terra::res(base)[1] / 10)
    grid_env <- terra::disagg(grid_env, factor)
  } else {
    grid_env <- terra::resample(grid_env, base, method = "near")
  }

  no_parts <- length(unique(occs_p$.part))
  part_pts <- list()
  bg_pts <- list()
  grid_parts <- list()

  for (i in 1:no_parts) {
    s <- Sys.time()
    part_pts[[i]] <- occs_p[occs_p$.part == i, ]
    grid_parts[[i]] <- terra::ifel(grid_env == i, 1, NA)
    bg_pts[[i]] <- sample_bg_upd(
      data = dplyr::tibble(part_pts[[i]]),
      x = "x",
      y = "y",
      n = length(part_pts[[i]]$x),
      method = c("thickening", width = 10000),
      rlayer = grid_parts[[i]],
      maskval = 1
    )
    bg_pts[[i]]$.part <- i
    t <- Sys.time() - s
    print(paste0(
      "Time taken to create background points for part ", i, ": ",
      t, " ", attr(t, "units")
    ))
  }

  bg_pts <- data.table::rbindlist(bg_pts, use.names = TRUE, fill = TRUE)
  # Combine the background points with the occurrences
  occs_p <- data.table::rbindlist(list(occs_p, bg_pts),
    use.names = TRUE,
    fill = TRUE
  )
  occs_p$d <- d

  # Write out the background points
  bg_path <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/background/thicken/"
  )
  dir.create(bg_path, showWarnings = FALSE, recursive = TRUE)

  data.table::fwrite(occs_p, paste0(
    bg_path, species, "_bg_pts_filtgeo_", d, "m2.csv"
  ), row.names = FALSE)
  t <- Sys.time() - s
  print(paste0(
    "Time taken to create background points: ", t, " ", attr(t, "units"),
    " for ", length(occs_p$x), " points for filtgeo", d, "m2"
  ))
  return(occs_p)
}

#' @description Function to create a target group bias file from
#' iNaturalist tree occurrence data. The function returns the target
#' group bias file to be used with the sample_background function.
#' @param path A path to the Data folder.
#' @param extent The base raster cropped to the extent of interest.
#' @param domain The domain of interest. Options are "USA" or "Global".
#' @param res The resolution of the base raster.

target_group_rbias_lyr <- function(path, extent, domain, res) {
  # Copy the tree data to the target group directory
  path2 <- gsub(
    "/.*",
    "/auto_arborist_cvpr2022_v0.15/data/tree_classification/inat/metadata/", # nolint
    path
  )

  tg_outpath <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/background/target_group/"
  ) # nolint

  dir.create(tg_outpath, showWarnings = FALSE, recursive = TRUE)

  fn <- paste0(tg_outpath, "TargetGroup_biasfile.tif")
  if (!file.exists(fn)) {
    if (domain == "USA") {
      file.copy(
        paste0(path2, "tree_observations_uscanada.csv"),
        paste0(tg_outpath, "tree_observations_uscanada.csv")
      )
    } else {
      file.copy(
        paste0(path2, "tree_observations_global.csv"),
        paste0(tg_outpath, "tree_observations_global.csv")
      )
    }

    # Read in the tree data
    trees <- data.table::fread(paste0(
      tg_outpath,
      "tree_observations_uscanada.csv"
    ))

    base_extent <- terra::ext(extent)

    # Crop the tree data to the extent of interest
    trees <- trees[longitude >= base_extent[1] & longitude <= base_extent[2] &
      latitude >= base_extent[3] & latitude <= base_extent[4], ]

    trees <- unique(trees[, .(longitude, latitude,
      sciname = paste0(genus, " ", species)
    )])

    trees <- unique(trees[, .(count = .N), by = .(longitude, latitude)])

    # Write out the tree data
    data.table::fwrite(trees, paste0(
      tg_outpath,
      "tree_observations_cropped.csv"
    ),
    row.names = FALSE
    )

    # Delete the original tree data
    file.remove(paste0(tg_outpath, "tree_observations_uscanada.csv"))
    
    trees <- data.table::fread(paste0(
      tg_outpath,
      "tree_observations_cropped.csv"
    ))

    # Calculate scale for the tree data
    scale <- 1 / sum(trees$count)
    trees$scale <- trees$count * scale
    trees_crds <- as.matrix(trees[, .(longitude, latitude)])

    if (res == 30) {
      base_extent <- terra::aggregate(extent, fact = 6)
      base_extent <- terra::app(base_extent, function(x) ifelse(is.na(x), 0, x))
      ncols <- ncol(base_extent)
      nrows <- nrow(base_extent)
      # make tiles
    } else {
      base_extent <- extent
    }
    target_density <- ks::kde(as.matrix(trees[, c("longitude", "latitude")]), # no lint
          w = trees$scale, gridsize = c(
            nrow(base_extent),
            ncol(base_extent)
          )
        )

    target_raster <- base_extent
    terra::values(target_raster) <- target_density$estimate

    # Normalize the target bias file between 0 and 1
    min <- terra::global(target_raster, "min")
    min <- min$min

    f <- function(i) i - min
    target_raster <- terra::app(target_raster, f)
    target_raster <- spatialEco::raster.transformation(target_raster,
      trans = "norm"
    )

    if (res == 30) {
      target_raster <- terra::resample(target_raster, extent)
    } else {
      target_raster <- target_raster
    }
    terra::writeRaster(target_raster, paste0(
      tg_outpath,
      "TargetGroup_biasfile.tif"
    ),
    overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=ZSTD"
    )
  } else {
    target_raster <- terra::rast(fn)
  }
  return(target_raster)
}

#' @description Function to create a population density bias file from
#' the Population Density raster. The function returns the population
#' density bias file to be used with the sample_background function.
#' @param path A path to the Data folder.
#' @param extent The base raster cropped to the extent of interest.
#' @param rbias_type The type of bias file to use. Options are
#' "pop_density", "distance_roads", or "distance_rails".
#' @param domain The domain of interest. Options are "USA" or "Global".
#' @param res The resolution of the base raster.

human_factors_rbias_lyr <- function(path, extent, domain, rbias_type, res) {
  # Create the output path
  hf_outpath <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/background/", rbias_type, "/"
  )
  dir.create(hf_outpath, showWarnings = FALSE, recursive = TRUE)
  fn <- paste0(hf_outpath, paste0(rbias_type, "_biasfile.tif"))

  if (!file.exists(fn)) {
    # Pull in the population density data from the Data folder
    if (rbias_type == "pop_density") {
      hf <- subset_rasters(path, domain,
        all = FALSE, temp = FALSE,
        precip = FALSE, topo = FALSE, land = FALSE, soils = FALSE, pop = TRUE,
        gdd = FALSE, biotic = FALSE, light = FALSE, rails = FALSE, roads = FALSE,
        downscaled = FALSE
      )
    } else if (rbias_type == "distance_roads") {
      hf <- subset_rasters(path, domain,
        all = FALSE, temp = FALSE,
        precip = FALSE, topo = FALSE, land = FALSE, soils = FALSE, pop = FALSE,
        gdd = FALSE, biotic = FALSE, light = FALSE, rails = FALSE, roads = TRUE,
        downscaled = FALSE
      )
    } else if (rbias_type == "distance_rails") {
      hf <- subset_rasters(path, domain,
        all = FALSE, temp = FALSE,
        precip = FALSE, topo = FALSE, land = FALSE, soils = FALSE, pop = FALSE,
        gdd = FALSE, biotic = FALSE, light = FALSE, rails = TRUE, roads = FALSE,
        downscaled = FALSE
      )
    }
    hf_files <- as.vector(unlist(lapply(hf, function(x) {
      get_filename(x, "USA", path)
    })))
    subset_files <- hf_files[grep(res, hf_files)]
    hf <- terra::rast(subset_files)
    # If the resolution is 30m, aggregate the population density to 100m
    if (res == 30) {
      hf <- terra::aggregate(hf, fact = 4)
    } else {
      hf <- hf}
    # Crop the population density to the extent of interest and write it out
    hf <- terra::crop(hf[[1]], extent, filename = paste0(
      hf_outpath, rbias_type, "_original_cropped.tif"
    ), wopt = list(
      gdal = "COMPRESS=ZSTD", datatype = "FLT4S",
      overwrite = TRUE
    ))
    # plot all values <1000
    hf2 <- terra::ifel(hf[[1]] < 1000, NA, hf)

    hf <- log(hf)

    if (rbias_type == "distance_roads" | rbias_type == "distance_rails") {
      remove_infs <- function(x) {
        x[x == -Inf] <- NA
        return(x)
      }
      hf <- terra::app(hf, remove_infs)
      hf <- terra::app(hf, function(x) -x)
    } else {
      hf <- hf
    }
    # Normalize bias file to between 0 and 1.
    min <- terra::global(hf, "min", na.rm = TRUE)$min
    hf <- terra::app(hf, function(x) x - min)
    hf <- spatialEco::raster.transformation(hf, trans = "norm")
    # Add 0.1 to the bias file to avoid 0 values
    hf <- hf + 0.1
    # Convert NA values to 0
    hf <- terra::app(hf, function(x) ifelse(is.na(x), 0, x))
    # Normalize
    hf <- spatialEco::raster.transformation(hf, trans = "norm")
    terra::writeRaster(hf, paste0(
      hf_outpath, rbias_type, "_biasfile.tif"
    ), overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=ZSTD")
    unlink(paste0(hf_outpath, rbias_type, "_original_cropped.tif"))
  } else {
    hf <- terra::rast(fn)
  }

  return(hf)
}


#' @description Function to create sample background points using the
#' @description A function to create thickened background points data
#' and write the data to the local output path.
#' @param partitioned_data A list of partitioned data occurrences, and the
#' associated spatial block grid.
#' @param rbias_lyr The bias file: target group raster, population raster, or
#' distance rasters.
#' @param rbias_type The type of bias file to use. Options are "target_group",
#' "pop_density", "distance_roads", or "distance_rails".
#' @param extent The base raster cropped to the extent of interest.
#' @param species The species name.
#' associated grid.
#' @param res The resolution of the base raster.
#' @param extent


create_biased_bg_pts <- function(
    partitioned_data, rbias_lyr,
    rbias_type, species, res, extent) {
  s0 <- Sys.time()
  occs_p <- partitioned_data[[1]]
  grid_env <- partitioned_data[[2]]

  # Disaggregate the spatial block grid
  if (res == 30) {
    #factor <- floor(terra::res(grid_env)[1] / terra::res(extent)[1] / 10)
    grid_env <- terra::resample(grid_env, rbias_lyr, method = "near")
  } else {
    grid_env <- terra::resample(grid_env, rbias_lyr, method = "near")
  }
  d <- unique(occs_p$d)
  occs_p <- occs_p[, .(pr_ab, x, y, .part)]
  no_parts <- length(unique(occs_p$.part))

  part_pts <- list()
  bg_pts <- list()
  grid_parts <- list()

  for (i in 1:no_parts) {
    s <- Sys.time()
    part_pts[[i]] <- occs_p[occs_p$.part == i, ]
    grid_parts[[i]] <- terra::ifel(grid_env == i, 1, NA)
    bg_pts[[i]] <- flexsdm::sample_background(
      data = dplyr::tibble(part_pts[[i]]),
      x = "x",
      y = "y",
      n = length(part_pts[[i]]$x),
      method = "biased",
      rlayer = grid_parts[[i]],
      maskval = 1,
      rbias = rbias_lyr
    )
    bg_pts[[i]]$.part <- i
    t <- Sys.time() - s
    print(paste0(
      "Time taken to create background points for part ", i, ": ",
      t, " ", attr(t, "units")
    ))
  }

  bg_pts <- data.table::rbindlist(bg_pts, use.names = TRUE, fill = TRUE)
  # Combine the background points with the occurrences
  occs_p <- data.table::rbindlist(list(occs_p, bg_pts),
    use.names = TRUE,
    fill = TRUE
  )
  occs_p$d <- d

  # Write out the background points
  bg_path <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/background/", rbias_type, "/"
  )
  dir.create(bg_path, showWarnings = FALSE, recursive = TRUE)

  data.table::fwrite(occs_p, paste0(
    bg_path, species, "_", rbias_type, "_bg_pts_filtgeo_", d, "m2.csv"
  ), row.names = FALSE)
  t0 <- Sys.time() - s0
  print(paste0(
    "Time taken to create background points: ", t0, " ", attr(t0, "units"),
    " for ", length(occs_p$x), " points for filtgeo", d, "m2"
  ))
  return(occs_p)
}

#' @description Function to conduct cluster analysis on the predictors
#' based on the correlation threshold. A sample_size is used to reduce
#' the time taken to perform the cluster analysis.
#' @param data A list of  predictors that have been z-score
#' normalized.
#' @param sample_size The number of samples to use in the cluster analysis.
#' @param mincor The minimum correlation threshold to use in the cluster analysis.
#' @param categorical_vars A list of categorical variables to include in the

cluster_analysis <- function(data, sample_size, mincor, categorical_vars) {
  # Perform cluster analysis
  tryCatch(
    {
      sample <- terra::spatSample(data,
        size = sample_size, method = "regular",
        as.df = TRUE, na.rm = TRUE, values = TRUE,
        xy = TRUE
      )

      # Keep all columns except 1:2 and cat_cols
      ccres <- klaR::corclust(sample[, -c(1:2)])
    },
    error = function(e) {
      message("An error occurred: ", e$message)
      message("Please increase the sample_size.")
    }
  )
  # Create hierarchical cluster tree based correlation threshold
  cluster_dt <- klaR::cvtree(ccres, mincor = mincor)
  vars <- rownames(cluster_dt$correlations)
  cluster_dt <- as.data.table(cluster_dt$correlations)
  cluster_dt$var <- vars

  # Reassign clusters
  cluster_dt[, clustercount := .N, by = .(cluster)]
  singles <- cluster_dt[cluster_dt$clustercount == 1 &
    cluster_dt$av.cor2closest < mincor]
  multis <- cluster_dt[cluster_dt$clustercount > 1]
  multisfix <- multis[av.cor2closest > mincor, ]
  multisfix[, cluster := pmin(cluster, closest)]
  multis <- multis[av.cor2closest < mincor, ]
  cluster_dt <- rbind(singles, multis, multisfix)
  cluster_dt <- cluster_dt[, .(var, cluster)]
  cluster_dt <- cluster_dt[order(cluster_dt$cluster, cluster_dt$var)]
  # Add categorical variables to the cluster assignments
  if (!is.null(categorical_vars)) {
    cat_vars <- names(categorical_vars)
    cat_clusters <- max(cluster_dt$cluster) + 1:length(cat_vars)
    cat_vars <- data.table(var = cat_vars, cluster = cat_clusters)
    cluster_dt <- rbind(cluster_dt, cat_vars)
  }
  return(cluster_dt)
}


#' @description Generate combinations of variables based on clusters
#' @param vars A character vector containing the variable names
#' @param clusters A character vector containing the cluster assignments
#' @return A list containing data frames of variable combinations for each
#' cluster count

library(dplyr)
library(purrr)

# Sample data frame
df <- data.frame(
  var = c("A", "B", "C", "D", "E", "F"),
  cluster = c(1, 1, 2, 2, 3, 3)
)

# Function to generate combinations
generate_combinations <- function(df) {
  # Ensure the data frame is sorted by cluster
  df <- dplyr::arrange(df, cluster)

  # Split the data frame by cluster
  split_df <- split(df$var, df$cluster)

  # Get unique clusters
  unique_clusters <- unique(df$cluster)

  # Helper function to generate all valid combinations
  generate_valid_combinations <- function(vars_by_cluster) {
    combs <- list()
    for (i in 1:length(vars_by_cluster)) {
      cluster_combinations <- combn(seq_along(vars_by_cluster), i, simplify = FALSE)
      for (cluster_comb in cluster_combinations) {
        valid_comb <- do.call(expand.grid, vars_by_cluster[cluster_comb])
        combs <- append(combs, list(valid_comb))
      }
    }
    combs
  }

  # Generate all valid combinations of variables from different clusters
  all_combinations <- generate_valid_combinations(split_df)

  # Combine all combinations into a single data frame
  final_combinations <- purrr::map_dfr(all_combinations, dplyr::as_tibble)

  return(final_combinations)
}

#' @description Function to append enviromental data to the partitioned
#' presence-background data
#' @param part_data A list of partitioned data presence-background data.
#' @param env_layer The list of predictors to extract.
#' @param bg_method The method used to sample the background points: random,
#' thicken, or biased: target_group, pop_density, distance_roads, or
#' distance_rails.
#' @param species The species name.

append_env_data <- function(part_data, env_layer, bg_method, species) {
  s <- Sys.time()
  d <- unique(part_data$d)
  # Convert to spatvect
  part_data <- terra::vect(part_data, crs = terra::crs(env_layer), geom = c("x", "y"))
  # Extract the environmental data
  part_data <- terra::extract(env_layer, part_data, ID = FALSE, xy = TRUE, bind = TRUE, method = "simple")
  part_data <- data.table::as.data.table(part_data)
  fwrite(part_data, paste0(
    getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/background/", bg_method, "/",
    species, "_", bg_method, "_bg_pts_filtgeo_", d, "m2_wpreds.csv"
  ), row.names = FALSE)
  t <- Sys.time() - s
  print(paste0(
    "Time taken to append predictor data: ", t, " ", attr(t, "units")
  ))
  terra::tmpFiles(orphan = TRUE, old = TRUE, remove = TRUE)
  return(part_data)
}

#' @description Function to pull in systematically sampled occurrence data
#' from NEON and VMI NPS Inventory data and create testing sets for
#' independent validation.
#' @param path The path to the Data folder.
#' @param extent The base raster cropped to the extent of interest.
#' @param species The species name.


create_testing_set <- function(path, extent, species) {
  test_data_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/original/")
  neon <- data.table::fread(paste0(test_data_path, species, "_test_occs.csv"))
  vmi <- data.table::fread(paste0(path, "/Original/vmi_nps_inventory/data/cleaned/csv/aial.csv"))
  vmi <- unique(vmi[, .(x = lon, y = lat, pr_ab = p_a)])
  vmi <- vmi[, .(pr_ab = max(pr_ab)), by = .(x, y)]
  test_data <- data.table::rbindlist(list(neon, vmi), fill = TRUE, use.names = TRUE)
  test_data <- test_data[, .(x, y, pr_ab)]
  test_data <- terra::vect(test_data, crs = terra::crs(extent), geom = c("x", "y"))
  test_data <- terra::crop(test_data, extent)
  test_data$x <- terra::crds(test_data)[, 1]
  test_data$y <- terra::crds(test_data)[, 2]
  pr_ab_1 <- data.table::as.data.table(test_data[test_data$pr_ab == 1])
  pr_ab_0 <- data.table::as.data.table(test_data[test_data$pr_ab == 0])
  # Randomly select length(pr_ab_1) points from pr_ab_0 ten times
  testing_sets <- list()

  for (i in 1:5) {
    pr_ab_0_sample <- pr_ab_0[sample(1:nrow(pr_ab_0), nrow(pr_ab_1)), ]
    testing_set <- data.table::rbindlist(list(pr_ab_1, pr_ab_0_sample), fill = TRUE)
    testing_sets[[i]] <- testing_set
    testing_sets[[i]]$setid <- i
  }
  testing_sets <- data.table::rbindlist(testing_sets, fill = TRUE)
  fwrite(testing_sets, paste0(
    getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/original/",
    species, "_testing_sets.csv"
  ), row.names = FALSE)
  return(testing_sets)
}

#' Function to setup sdm prepped datasets for
#' analysis.
#' @param data The target group background points data
splitpts4sdm <- function(data) {
  #' Remove unnecessary columns: d, x, y using data.table
  data <- data[, c("d", "x", "y") := NULL]
  #' Rename NLCD Land Cover Class to landcoverrc
  setnames(data, "NLCD Land Cover Class", "landcoverrc")
  #' Convert to tibble
  data <- tibble::as_tibble(data)
  #' Train on .part == 1:4; test on .part == 5
  #no_parts <- length(unique(data$.part))
  #data <- data[data$.part %in% 1:(no_parts - 1), ]
  #' Pull out presence-absence data
  response <- data[data$pr_ab == 1, ]
  #' Add id column
  response$id <- 1:nrow(response)
  #' Pull out background data; pr_ab == 0
  bg <- data[data$pr_ab == 0, ]
  return(list(response = response, bg = bg))
}


#' @description Function to aggregate transformed predictors by a factor of 3
#' to speed up the analysis.
#' @param predictor A cropped and transformed predictor.
#' @param factor The factor to aggregate the predictor by.

aggregate_predictor <- function(predictor, factor) {
  s <- Sys.time()
  fn <- basename(terra::sources(predictor))
  fn <- gsub(".tif", "_agg.tif", fn)
  dt <- terra::datatype(predictor)
  agg_outpath <- paste0(
      getwd(),
      "/flexsdm_results/1_Inputs/2_Predictors/1_Current/cropped/transformed/aggregated/" # nolint
    )
    dir.create(agg_outpath, showWarnings = FALSE, recursive = TRUE)
  if (dt != "INT1U") {
    predictor <- terra::aggregate(predictor, fact = factor, fun = "mean")
    terra::writeRaster(predictor, paste0(
      agg_outpath, fn
    ), overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=ZSTD")
  } else {
    predictor <- NULL
  }
  t <- Sys.time() - s
  print(paste0(
    "Time taken to aggregate predictor: ", fn, " ", t, " ", attr(t, "units")
  ))
  return(predictor)
}

#####################################################################################
######
###### Formerly a separate script called get_envi_chunked.R
######
# Load packages
require(geodata)
require(terra)
require(sbtools)
require(rnaturalearth)
require(zen4R)
require(rgee)
require(foreach)

#' Function to get biovariables for the world (worldclim 30s)
#' @param path A character string specifying the path to write the biovars
#' @return A raster with 19 biovariable layers for the world.
#' @export

get_biovars_global <- function(path) {
  fn <- paste0(path, "Original/wc2.1_30s_bio/wc2.1_30s_bio_", 1:19, ".tif")
  if (!all(file.exists(fn))) {
    biovars <- geodata::worldclim_global(
      var = "bio", res = 0.5,
      path = tempdir()
    )
    dir.create(dirname(fn[1]), showWarnings = FALSE)
    # Write each layer to file
    lapply(1:terra::nlyr(biovars), function(x) {
      terra::writeRaster(biovars[[x]],
        overwrite = TRUE, gdal = "COMPRESS=ZSTD",
        filename = paste0(
          dirname(fn[1]),
          names(biovars)[x], ".tif"
        ),
        datatype = "FLT4S"
      )
    })
  } else {
    biovars <- lapply(fn, terra::rast)
  }
  return(biovars)
}

#' @description Function to get biovariables for CONUS (USGS biocomposite 15s)
#' @param path A character string specifying the path to write the biovars.
#' @return A raster with 19 biovariable layers for CONUS.
#' @export

get_biovars_conus <- function(path) {
  bio_names <- c(paste0("bio", 1:4), "bio4a", paste0("bio", 5:19))
  fn <- paste0(path, "Original/BioClimComposite/BioClimComposite_1971_2000_400m_", bio_names, ".tif") # nolint: line_length_linter.
  if (!all(file.exists(fn))) {
    sbtools::item_file_download(
      sb_id = "4ff32906e4b0e183ef5a2f16", dest_dir =
        tempdir(), names = "BioClimComposite_1971_2000_400m.tif", # nolint
      destinations = paste0(tempdir(), "BioClimComposite_1971_2000_400m.tif"),
      overwrite_file = TRUE
    )
    biovars <- terra::rast(paste0(
      tempdir(),
      "BioClimComposite_1971_2000_400m.tif" # nolint
    ))
    names(biovars) <- c(paste0("bio", 1:4), "bio4a", paste0("bio", 5:19))
    dir.create(dirname(fn[1]), showWarnings = FALSE)
    # Write each layer to file
    lapply(1:terra::nlyr(biovars), function(x) {
      terra::writeRaster(biovars[[x]],
        overwrite = TRUE, gdal = "COMPRESS=ZSTD",
        filename = paste0(
          dirname(fn[1]), "/BioClimComposite_1971_2000_400m_",
          names(biovars)[x], ".tif"
        ),
        datatype = "FLT4S"
      )
    })
  } else {
    biovars <- lapply(fn, terra::rast)
  }
  # Remove bio4a
  biovars <- biovars[-5]
  return(biovars)
}

#' @description Function to get biovariables for CONUS (daymet 1km downscaled
#' using elevation in regression)
#' @param path A character string specifying the path to write the biovars.
#' @return A raster with 19 biovariable layers for CONUS.
#' @export

get_biovars_conus_ds <- function(path) {
  bio_names <- c(paste0("bio", 1:19))
  fn <- paste0(path, "Raster/USA/bioclimatic/", bio_names, ".tif")
  biovars <- lapply(fn, terra::rast)
  return(biovars)
}


#' @description Function to obtain a vector of CONUS.
#' @return A vector of CONUS.
#' @export

get_l48_boundary <- function() {
  locs <- rnaturalearth::ne_states(
    country = "United States of America",
    returnclass = "sf"
  )
  locs <- unique(locs["iso_3166_2"])
  locs <- locs[locs$iso_3166_2 != "US-PR" & locs$iso_3166_2 != "US-VI" &
    locs$iso_3166_2 != "US-GU" & locs$iso_3166_2 != "US-MP" &
    locs$iso_3166_2 != "US-AS" & locs$iso_3166_2 != "US-AK" &
    locs$iso_3166_2 != "US-HI", ]
  locs <- sf::st_union(locs)
  return(locs)
}

#' @description Function to create terrain variables (slope, aspect,
#' & hillshade) from an elevation layer.
#' @param path A character string with the path to write the terrain variables.
#' @return A hillshade raster.
#' @export

get_terrain <- function(path, elevation) {
  fn <- paste0(c("slope", "aspect", "hillshade"), ".tif")
  # Calculate slope & export
  slp <- terra::terrain(elevation, v = "slope", unit = "radians", neighbors = 8)
  names(slp) <- "slope"
  terra::writeRaster(slp,
    filename = paste0(path, fn[1]), overwrite = TRUE, datatype = "FLT4S",
    gdal = "COMPRESS=ZSTD"
  )
  # Calculate aspect & export
  asp <- terra::terrain(elevation,
    v = "aspect", unit = "radians",
    neighbors = 4
  )
  names(asp) <- "aspect"
  terra::writeRaster(asp,
    filename = paste0(path, fn[2]), overwrite = TRUE, datatype = "FLT4S",
    gdal = "COMPRESS=ZSTD"
  )
  # Calculate hillshade & export
  shd <- terra::shade(slp, asp, direction = 180)
  names(shd) <- "hillshade"
  terra::writeRaster(shd,
    filename = paste0(path, fn[3]), overwrite = TRUE, datatype = "FLT4S",
    gdal = "COMPRESS=ZSTD"
  )
  return(list(slp, asp))
}

#' @description Function to obtain a global elevation raster (30s WorldClim).
#' and to calculate slope, aspect, and hillshade from the elevation layer.
#' @param path A character string specifying the path to write the elevation.
#' @return A list with elevation & hillshade rasters.
#' @export

get_topo_global <- function(path) {
  fn <- c(
    paste0(path, "Original/wc2.1_30s_elev/wc2.1_30s_elev.tif"),
    paste0(path, "Raster/Global/wc2.1_30s_elev/wc2.1_30s_slope.tif"),
    paste0(path, "Raster/Global/wc2.1_30s_elev/wc2.1_30s_aspect.tif")
  )
  if (!all(file.exists(fn))) {
    elevation <- geodata::elevation_global(res = 0.5, path = tempdir())
    dir.create(dirname(fn[1]), showWarnings = FALSE, recursive = TRUE)
    dir.create(dirname(fn[2]), showWarnings = FALSE, recursive = TRUE)
    terra::writeRaster(elevation,
      overwrite = TRUE, gdal = "COMPRESS=ZSTD",
      filename = fn, datatype = "INT2S"
    )
    asp_slope <- get_terrain(
      path = dirname(fn[2]),
      elevation
    )
    topo <- list(elevation, asp_slope[[1]], asp_slope[[2]])
  } else {
    topo <- lapply(fn, terra::rast)
  }
  return(topo)
}

#' @description Function to get a elevation layer for CONUS (NASA SRTM 1s),
#' and to calculate slope, aspect, and hillshade from the elevation layer.
#' @param key A character string specifying the OpenTopography API key.
#' @param path A character string specifying the path to write the elevation.
#' @export

get_topo_conus <- function(path, key) {
  fn <- c(
    paste0(path, "Original/SRTMGL1/SRTMGL1_1s_l48.tif"),
    paste0(path, "Raster/USA/SRTMGL1/SRTMGL1_1s_l48_slope.tif"),
    paste0(path, "Raster/USA/SRTMGL1/SRTMGL1_1s_l48_aspect.tif")
  )
  if (!all(file.exists(fn))) {
    # Get bounding box of CONUS
    loc <- get_l48_boundary()
    # Make grid of 450,000 sq km cells
    loc_bb <- sf::st_transform(loc, crs = 3857)
    loc_bb <- sf::st_make_grid(loc_bb, cellsize = 670820)
    loc_bb <- sf::st_transform(loc_bb, crs = 4326)
    loc_bb <- lapply(seq_along(loc_bb), function(x) {
      loc_bb[[x]] <- sf::st_bbox(loc_bb[[x]])
    })
    # Download tiles
    dir.create(paste0(dirname(fn[1]), "/tiles"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    lapply(seq_along(loc_bb), function(x) {
      url <- paste0(
        "https://portal.opentopography.org/API/",
        "globaldem?demtype=", "SRTMGL1", "&south=", loc_bb[[x]][2],
        "&north=", loc_bb[[x]][4], "&west=", loc_bb[[x]][1], "&east=",
        loc_bb[[x]][3], "&outputFormat=GTiff&API_Key=", key
      )
      # Use curl to request download url
      dl_url <- curl::curl_fetch_memory(url)$url
      if (.Platform$OS.type == "windows") {
        download.file(dl_url, destfile = paste0(
          path, "Original/SRTMGL1/tiles/srtm_",
          x, ".tif"
        ), mode = "wb")
      } else {
        download.file(dl_url, destfile = paste0(
          path, "Original/SRTMGL1/tiles/srtm_", x,
          ".tif"
        ))
      }
    })
    # Read in tiles
    srtm_files <- list.files(paste0(path, "Original/SRTMGL1/tiles"),
      pattern = ".tif",
      full.names = TRUE
    )
    srtm_files <- terra::sprc(srtm_files)
    srtm <- terra::mosaic(srtm_files, wopt = list(datatype = "INT2S"))
    names(srtm) <- "elevation"
    srtm <- terra::crop(srtm, terra::vect(loc),
      mask = TRUE, filename = fn[1],
      overwrite = TRUE, wopt = list(
        datatype = "INT2S",
        gdal = "COMPRESS=ZSTD"
      )
    )
    # Delete tiles folder
    unlink(paste0(path, "Original/SRTMGL1/tiles"), recursive = TRUE)
    # Calculate hillshade
    dir.create(dirname(fn[2]), showWarnings = FALSE, recursive = TRUE)
    asp_slope <- get_terrain(
      path = paste0(dirname(fn[2]), "/SRTMGL1_1s_l48_"),
      srtm
    )
    topo <- list(srtm, asp_slope[[1]], asp_slope[[2]])
  } else {
    topo <- lapply(fn, terra::rast)
  }
  return(topo)
}

#' @description Function to get landcover for the world (30s ESA World Cover)
#' @param path A character string specifying the path to write the landcover
#' @return A list with landcover rasters: built, cropland, grassland, shrubs,
#' trees, wetland.
#' @export

get_landcover_global <- function(path) {
  types <- c("built", "cropland", "grassland", "shrubs", "trees", "wetland")
  fn <- paste0(path, "Original/WorldCover/", types, "_30s.tif")
  if (!all(file.exists(fn))) {
    landcover <- lapply(seq_along(types), function(x) {
      geodata::landcover(var = types[[x]], path = tempdir())
    })
    dir.create(dirname(fn[1]), showWarnings = FALSE)
    lapply(seq_along(fn), function(x) {
      terra::writeRaster(landcover[[x]],
        overwrite = TRUE,
        gdal = "COMPRESS=ZSTD", filename = fn[x], datatype = "INT1U"
      )
    })
  } else {
    landcover <- lapply(seq_along(fn), function(x) {
      terra::rast(fn[x])
    })
  }
  return(landcover)
}
#' @description Function to obtain landcover data for CONUS (30m NLCD)
#' @param path A character string specifying the path to write the landcover
#' @return A list with landcover rasters: built, deciduous, evergreen, trees,
#' shrubs, grass, pasture, cropland, cultivated, wetland.
#' @export

get_landcover_conus <- function(path) {
  fn <- c(
    paste0(path, "Original/nlcd/nlcd_2021_land_cover_l48_20230630.tif"),
    paste0(path, "Raster/USA/nlcd/nlcd_2021_land_cover_l48.tif")
  )
  if (!all(file.exists(fn))) {
    url <- "https://s3-us-west-2.amazonaws.com/mrlc/nlcd_2021_land_cover_l48_20230630.zip" # nolint: line_length_linter.
    # Download zip file
    options(timeout = 60 * 5)
    download.file(url,
      destfile = paste0(tempdir(), "/", basename(url)),
      mode = "wb"
    )
    unzip(paste0(tempdir(), "/", basename(url)), exdir = paste0(
      tempdir(),
      "/nlcd"
    ))
    # Read in raster
    landcover <- terra::rast(paste0(tempdir(), "/nlcd/nlcd_2021_land_cover_l48_20230630.img")) # nolint: line_length_linter.
    dir.create(dirname(fn[1]),
      showWarnings = FALSE,
      recursive = TRUE
    )
    terra::writeRaster(landcover,
      overwrite = TRUE, gdal = "COMPRESS=ZSTD",
      filename = fn, datatype = "INT1U"
    )
    m <- c(
      0, 20, 0, 21, 24, 1, 25, 40, 0, 41, 41, 2, 42, 42, 3, 43, 43, 4, 44,
      51, 0, 52, 52, 5, 53, 70, 0, 71, 71, 6, 72, 80, 0, 81, 81, 7, 82,
      82, 8, 83, 89, 0, 90, 90, 9, 91, 94, 0, 95, 95, 9, 96, 255, 0
    )
    rclmat <- matrix(m, ncol = 3, byrow = TRUE)
    landcover <- terra::classify(landcover, rclmat, right = NA)
    dir.create(dirname(fn[2]), showWarnings = FALSE, recursive = TRUE)
    output_file <- paste0(dirname(fn[2]), "/nlcd_2021_land_cover_l48.tif")
    l48 <- get_l48_boundary()
    l48 <- sf::st_transform(l48, crs = terra::crs(landcover))
    terra::writeRaster(landcover,
      overwrite = TRUE, gdal = "COMPRESS=ZSTD",
      filename = output_file, datatype = "INT1U"
    )
  } else {
    landcover <- lapply(fn[-1], terra::rast)
  }
  return(landcover)
}

#' @description Function to obtain monthly tavg  for the world (worldclim 30s)
#' and then calculate growing degree days.
#' @param path A character string specifying the path to write the tavg normals
#' @export

get_gdd_global <- function(path) {
  fn <- c(
    paste0(path, "Original/wc2.1_30s_tavg/wc2.1_30s_tavg.tif"), # original
    paste0(path, "Raster/Global/wc2.1_30s_tavg/wc2.1_30s_gdd.tif")
  ) # gdd
  if (!all(file.exists(fn))) {
    tavg <- geodata::worldclim_global(var = "tavg", res = 0.5, path = tempdir())
    dir.create(dirname(fn[1]), showWarnings = FALSE, recursive = TRUE)
    terra::writeRaster(tavg,
      overwrite = TRUE, gdal = "COMPRESS=ZSTD",
      filename = fn[1], datatype = "FLT4S"
    )
    tbase <- 5
    gdd <- terra::tapp(x = tavg, tbase = tbase, fun = function(tavg, tbase) {
      days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      xbase <- tavg - tbase
      xbase[(xbase < 0)] <- 0
      xbase <- sum(xbase * days)
      return(xbase)
    })
    dir.create(dirname(fn[2]), showWarnings = FALSE, recursive = TRUE)
    names(gdd) <- "gdd_tbase5"
    terra::writeRaster(gdd,
      filename = fn[2], overwrite = TRUE, datatype = "FLT4S",
      gdal = "COMPRESS=ZSTD"
    )
  } else {
    gdd <- terra::rast(fn[2])
  }
  return(gdd)
}

#' @description Function to retrieve gdd for CONUS using daymet 1km downscaled
#' @param path A character string specifying the path to the gdd raster
#' @return A raster with growing degree days for CONUS.
#' @export

get_gdd_conus_ds <- function(path) {
  fn <- c(paste0(path, "Raster/USA/gdd/gdd.tif"))
  gdd <- lapply(fn, terra::rast)
}


#' @description Function to obtain monthly precip normals for the world
#' (Worldclim 30s)
#' Monthly precip normals are used to calculate precipitation timing.
#' @param path A character string specifying the path to write precip normals
#' @export
#' @note This function is not currently used in the package, because it isn't
#' computing the precipitation timing correctly. The function is left here for
#' future development.

# get_prectiming_global <- function(path) {
#   fn <- c(paste0(path, "Original/wc2.1_30s_prec/wc2.1_30s_prec.tif"),
#           paste0(path, "Raster/Global/wc2.1_30s_prec/wc2.1_30s_prec_timing.tif")) # nolint: line_length_linter.
#   if (!all(file.exists(fn))) {
#     precip <- geodata::worldclim_global(var = "prec", res = 0.5,
#                                         path = tempdir())
#     dir.create(dirname(fn[1]), showWarnings = FALSE, recursive = TRUE)
#     terra::writeRaster(precip,
#       overwrite = TRUE, gdal = "COMPRESS=ZSTD",
#       filename = fn[1], datatype = "INT2S"
#     )
#     precip <- terra::app(x = precip, fun = function(x) {
#       return(sum(x[[12]], x[[1]], x[[2]]) - sum(x[[6]], x[[7]], x[[8]]))
#     })
#     names(precip) <- "PrecipTiming"
#     precip <- terra::writeRaster(precip,
#       filename = paste0(dirname(fn[2]), "/wc2.1_30s_prec_timing.tif"),
#       overwrite = TRUE, datatype = "INT2S", gdal = "COMPRESS=ZSTD"
#     )
#   } else {
#     precip <- terra::rast(fn[2])
#   }
#   return(precip)
# }

#' @description Function to download human pop. density for the world
#' (30s GHS_POP_E2020_GLOBE)
#' @param path A character string specifying the path to write the pop data
#' @return A raster with human population density for the world
#' @export

get_pop_global <- function(path) {
  fn <- paste0(path, "Original/gpw_v4_population_density/gpw_v4_population_density_rev11_2020_30s.tif") # nolint: line_length_linter.
  if (!file.exists(fn)) {
    pop <- geodata::population(year = 2020, res = 0.5, path = tempdir())
    dir.create(dirname(fn), showWarnings = FALSE)
    terra::writeRaster(pop,
      overwrite = TRUE, gdal = "COMPRESS=ZSTD",
      filename = fn, datatype = "FLT4S"
    )
  } else {
    pop <- terra::rast(fn)
  }
  return(pop)
}

#' @description Function to download human pop. density for conus
#' @param path A character string specifying the path to write the pop data
#' @return A raster with human population density for the world
#' @export

get_pop_conus <- function(path) {
  fn <- paste0(path, "Original/meta_population/meta_population.tif")
  pop <- lapply(fn, terra::rast)
  return(pop)
}

#' @description Function to calc distance to roads/railroads for the world (30s)
#' @param lines A sf object of roads/railroads.
#' @param domain A character string specifying the domain (world or conus)
#' @param type A character string specifying the type of network (roads/rails)
#' @return A raster with distance to roads or railroads for the world
#' @import rgee
#' @export
dist_networks <- function(user, lines, domain, type) {
  ee_Initialize(user, drive = TRUE)
  ee_Authenticate()
  # Remove all columns except geometry
  lines <- lines[, 1]
  # Split into 1200 chunks if type == rails, else split into 3000 chunks
  if (type == "rails") {
    lines <- split(lines, ceiling(1:nrow(lines) / 1200)) # nolint
  } else {
    lines <- split(lines, ceiling(1:nrow(lines) / 3000)) # nolint
  }
  # Create rgee folder in assets
  asset <- ee_get_assethome()
  ee_manage_create(paste0(asset, "/rgee/"))
  # Merge the assets
  assetids <- paste0(asset, "/rgee/", type, "_", 1:length(lines)) # nolint
  # Convert sf to ee$FeatureCollection
  lapply(seq_along(lines), function(x) {
    sf_as_ee(
      x = lines[[x]],
      via = "getInfo_to_asset",
      assetId = assetids[x],
      overwrite = TRUE,
      monitoring = TRUE,
      proj = "EPSG:4326"
    )
  })

  assetids <- ee_manage_assetlist(paste0(ee_get_assethome(), "/rgee/"))$ID

  mergefeaturecollections <- function(assetids) {
    # Start with an empty FeatureCollection
    mergedcollection <- ee$FeatureCollection(list())

    # Merge each asset into the mergedCollection
    for (assetid in assetids) {
      collection <- ee$FeatureCollection(assetid)
      mergedcollection <- mergedcollection$merge(collection)
    }
    return(mergedcollection)
  }

  combinedassets <- mergefeaturecollections(assetids)

  # Load the domain boundaries
  if (domain == "world") {
    domainboundaries <- ee$FeatureCollection("USDOS/LSIB_SIMPLE/2017")
  } else if (domain == "conus") {
    # Load CONUS boundaries
    conus <- get_l48_boundary()
    sf_as_ee(
      x = conus,
      via = "getInfo_to_asset",
      assetId = paste0(ee, "/rgee/conus"),
      overwrite = TRUE,
      monitoring = TRUE,
      proj = "EPSG:4326"
    )
    # Replace with appropriate FeatureCollection for CONUS
    domainboundaries <- ee$FeatureCollection(paste0(ee, "/rgee/conus"))
  }

  # Calculate the distance to the nearest feature
  # If statement to determine distance value
  if (domain == "conus") {
    distancetoassets <- combinedassets$distance(ee$Number(2000252))
  } else if (type == "world") {
    distancetoassets <- combinedassets$distance(ee$Number(20002520))
  } else {
    stop("Invalid domain")
  }

  # Crop distanceToAssets with domainboundaries
  croppeddistance <- distancetoassets$clip(domainboundaries)

  # Convert the cropped image to UInt32
  croppeddistanceuint32 <- croppeddistance$toUint32()

  # Export to Google Drive
  exporttask <- ee$batch$Export$image$toDrive(
    image = croppeddistanceuint32,
    description = paste0("dist2", type, "_", domain),
    folder = "eu_distancee",
    fileNamePrefix = paste0("dist2", type, "_", domain),
    scale = 1000,
    fileFormat = "GeoTIFF",
    maxPixels = 1e13
  )
  exporttask$start()
  while (exporttask$status()$state != "COMPLETED") {
    Sys.sleep(300)
  }
}

#' @description Function to get roads data for the world (10m rnaturalearth)
#' @param path A character string specifying the path to write the roads
#' @param gee_path A character string specifying the path to write the roads
#' @param user A character string specifying the Google Earth Engine user name.
#' @return A euclidean distance raster with roads for the world.
#' @export

get_roads_global <- function(path, gee_path, user) {
  fn <- c(
    paste0(path, "Original/ne_roads/ne_10m_roads.gpkg"),
    paste0(path, "Raster/Global/ne_roads/ne_roads_distance_world_30s.tif")
  )
  if (!all(file.exists(fn))) {
    roads <- rnaturalearth::ne_download(type = "roads", scale = 10)
    dir.create(dirname(fn[1]), showWarnings = FALSE, recursive = TRUE)
    sf::st_write(roads, dsn = fn[1])
    dist_networks(user, lines = roads, domain = "world", type = "roads")
    roads <- list.files(gee_path,
      pattern = "dist2roads_world",
      full.names = TRUE
    )
    roads <- terra::sprc(roads)
    roads <- terra::mosaic(roads)
    names(roads) <- "dist2roads"
    dir.create(dirname(fn[2]), showWarnings = FALSE, recursive = TRUE)
    terra::writeRaster(roads,
      filename = fn[2], overwrite = TRUE, datatype = "INT4U",
      gdal = "COMPRESS=ZSTD"
    )
  } else {
    roads <- terra::rast(fn[2])
  }
  return(roads)
}

#' @description Function to get roads data for the conus
#' @param path A character string specifying the path to write the roads
#' @return A euclidean distance raster with roads for the conus.
#' @export

get_roads_conus <- function(path) {
  fn <- paste0(path, "Raster/USA/tiger_roads/tiger_roads_dist2roads.tif")
  roads <- lapply(fn, terra::rast)
  return(roads)
}

#' @description Function to download railroads for the world (10m rnaturalearth)
#' @param path A character string specifying the path to write the railroads
#' @param gee_path A character string specifying the path to write the railroads
#' @param user A character string specifying the Google Earth Engine user name.
#' @return A euclidean distance raster with railroads for the world.
#' @export

get_rails_global <- function(path, gee_path, user) {
  fn <- c(
    paste0(path, "Original/ne_railroads/ne_10m_railroads.gpkg"),
    paste0(path, "Raster/Global/ne_railroads/ne_railroads_dist_world_30s.tif")
  )
  if (!all(file.exists(fn))) {
    railroads <- rnaturalearth::ne_download(type = "railroads", scale = 10)
    dir.create(dirname(fn[1]), showWarnings = FALSE)
    sf::st_write(railroads, dsn = fn[1])
    dist_networks(user, lines = railroads, domain = "world", type = "rails")
    railroads <- list.files(gee_path,
      pattern = "dist2rails_world",
      full.names = TRUE
    )
    railroads <- terra::sprc(railroads)
    railroads <- terra::mosaic(railroads)
    names(railroads) <- "dist2rails"
    dir.create(dirname(fn[2]), showWarnings = FALSE, recursive = TRUE)
    terra::writeRaster(railroads,
      filename = fn[2], overwrite = TRUE, datatype = "INT4U",
      gdal = "COMPRESS=ZSTD"
    )
  } else {
    railroads <- terra::rast(fn[2])
  }
  return(railroads)
}

#' @description Function to get rails data for the conus
#' @param path A character string specifying the path to write the rails
#' @return A euclidean distance raster with rails for the conus.
#' @export
#'
get_rails_conus <- function(path) {
  fn <- paste0(path, "Raster/USA/tiger_railroads/tiger_railroads_dist2rails.tif")
  railroads <- lapply(fn, terra::rast)
  return(railroads)
}

#' @description Function to get soil pH & wc 33kpa & wc 1500kpa for the world
#' (250m OpenLandMap)
#' @param path A character string specifying the path to write soil vars
#' @param return A list pH & wc 33kpa & wc 1500kpa rasters for the world
#' @export

get_soilvars_global <- function(path) {
  # Create file names
  depths <- c(
    "b0..0cm", "b10..10cm", "b30..30cm", "b60..60cm", "b100..100cm",
    "b200..200cm"
  )
  orig_root <- "Original/OpenLandMap_soils/"
  fn <- c(
    paste0(
      path, orig_root, "sol_ph.h2o_usda.4c1a2a_m_250m_", depths,
      "_1950..2017_v0.2.tif"
    ),
    paste0(path, orig_root, "sol_watercontent.33kPa_usda.4b1c_m_250m_", depths, "_1950..2017_v0.1.tif"), # nolint: line_length_linter
    paste0(path, orig_root, "sol_watercontent.1500kPa_usda.3c2a1a_m_250m_", depths, "_1950..2017_v0.1.tif"), # nolint: line_length_linter
    paste0(path, "Raster/Global/OpenLandMap_soils/sol_ph.h2o_usda.4c1a2a_m_250m_mean_depth_1950..2017_v0.2.tif"), # nolint: line_length_linter
    paste0(path, "Raster/Global/OpenLandMap_soils/sol_watercontent.33kPa_usda.4b1c_m_250m_mean_depth_1950..2017_v0.1.tif"), # nolint: line_length_linter
    paste0(path, "Raster/Global/OpenLandMap_soils/sol_watercontent.1500kPa_usda.3c2a1a_m_250m_mean_depth_1950..2017_v0.1.tif")
  ) # nolint: line_length_linter
  # check if file exists
  if (!all(file.exists(fn))) {
    # Download soil pH and water content
    dl_urls <- c(
      "https://doi.org/10.5281/zenodo.2525664",
      "https://doi.org/10.5281/zenodo.2784001"
    )
    # Get list of files
    doi <- gsub("https://doi.org/", "", dl_urls)
    names(dl_urls) <- c("ph", "watercontent")
    lapply(seq_along(dl_urls), function(x) {
      files <- zen4R::get_zenodo(doi[x])
      files <- files$files
      files <- unlist(lapply(seq_along(files), function(y) {
        files[[y]]$filename[1]
      }))
      files <- files[grepl(".tif$", files)]
      files <- files[!grepl("_md_", files)]
      dir.create(paste0(tempdir(), "/OpenLandMap_soils/", names(dl_urls[x])),
        showWarnings = FALSE, recursive = TRUE
      )
      dir.create(paste0(path, "Original/OpenLandMap_soils/", names(dl_urls[x])),
        showWarnings = FALSE, recursive = TRUE
      )
      zen4R::download_zenodo(dl_urls[x],
        path = paste0(tempdir(), "/OpenLandMap_soils/", names(dl_urls)[x]),
        timeout = 60 * 30, files = list(files)
      )
    })
    # List files of .tifs in tempdir
    files <- list.files(paste0(tempdir(), "/OpenLandMap_soils/"),
      full.names = TRUE, recursive = TRUE
    )
    files <- files[grepl(".tif$", files)]
    files <- sort(files)
    dir.create(dirname(fn[1]), showWarnings = FALSE, recursive = TRUE)
    # For each file, check if pH or water content, & write to correct folder
    for (x in seq_along(files)) {
      filename <- paste0(dirname(fn[1]), "/", basename(files[x]))
      terra::writeRaster(terra::rast(files[x]),
        overwrite = TRUE,
        gdal = "COMPRESS=ZSTD", filename = filename,
        datatype = "INT1U"
      )
    }
    soilvars <- list.files(dirname(fn[1]), full.names = TRUE, recursive = TRUE)
    dir.create(dirname(fn[grep("Raster", fn)][1]),
      showWarnings = FALSE,
      recursive = TRUE
    )
    # Compute mean for ph and water content
    variables <- c("_ph", "33kPa", "1500kPa")
    output_filenames <- c(
      "sol_ph.h2o_usda.4c1a2a_m_250m_mean_depth_1950..2017_v0.2.tif", # nolint: line_length_linter
      "sol_watercontent.33kPa_usda.4b1c_m_250m_mean_depth_1950..2017_v0.1.tif", # nolint: line_length_linter
      "sol_watercontent.1500kPa_usda.3c2a1a_m_250m_mean_depth_1950..2017_v0.1.tif"
    ) # nolint: line_length_linter
    soilvars <- lapply(seq_along(variables), function(i) {
      var <- variables[i]
      raster <- terra::rast(soilvars[grepl(var, soilvars)])
      raster <- terra::mean(raster, na.rm = TRUE)
      names(raster) <- var
      terra::writeRaster(raster,
        overwrite = TRUE, gdal = "COMPRESS=ZSTD",
        filename = paste0(path, "Raster/Global/OpenLandMap_soils/", output_filenames[i]), # nolint: line_length_linter
        datatype = "INT1U"
      )
      return(raster)
    })
  } else {
    soilvars <- lapply(fn, terra::rast)
  }
  return(soilvars)
}

#' @description Function to get soil characteristic data for the conus
#' @param path A character string specifying the path to write soil chars
#' @return A euclidean distance raster with soil chars for the conus.
#' @export

get_soilvars_conus <- function(path) {
  chars <- c("clay", "om", "ph", "sand", "silt", "theta_r", "theta_s")
  fn <- paste0(path, "Original/polaris_soils/", chars, "_mean_0_5.tif")
  soilvars <- lapply(fn, terra::rast)
  return(soilvars)
}

#' @description Function to get PAR for conus (500m MODIS)
#' @param path A character string specifying the path to write PAR
#' @return A raster with PAR for the conus.
#' @export

get_par_conus <- function(path) {
  fn <- paste0(path, "Original/modis_par/modis_par.tif")
  par <- lapply(fn, terra::rast)
  return(par)
}

#' @description Function to get landsat ndvi data for the conus (30m)
#' @param path A character string specifying the path to write landsat evi
#' @return A list with landsat ndvi rasters for the conus.
#' @export

get_landsat_ndvi_conus <- function(path) {
  fn <- paste0(path, "Original/landsat_ndvi/landsat_ndvi.tif")
  ndvi <- lapply(fn, terra::rast)
  return(ndvi)
}

#' @description Function to get landsat evi data for the conus (30m)
#' @param path A character string specifying the path to write landsat evi
#' @return A list with landsat ndvi rasters for the conus.
#' @export

get_landsat_evi_conus <- function(path) {
  fn <- paste0(path, "Original/landsat_evi/landsat_evi.tif")
  evi <- lapply(fn, terra::rast)
  return(evi)
}

#' @description Function to create file names for resampled and reprojected predictors
#' @param pred A predictor raster
#' @param domain A character string specifying the domain (Global or USA)
#' @param path A character string specifying the path to write the reprojected

get_filename <- function(pred, domain, path) {
  # Get filename
  filename <- names(pred)
  dir <- gsub(".*\\/", "", dirname(terra::sources(pred)))
  filename <- paste0(dir, "/", dir, "_", filename)

  # Create filename based on domain and crop
  if (domain == "Global") {
    common_res <- c(1000, 2500, 5000)
    filename <- paste0(
      path, "Raster/Global/", filename, "_reproj_",
      common_res, "m.tif"
    )
    dir.create(dirname(filename[1]), showWarnings = FALSE, recursive = TRUE)
  } else if (domain == "USA") {
    common_res <- c(30, 100, 250, 500, 1000)
    filename <- paste0(
      path, "Raster/USA/", filename, "_reproj_",
      common_res, "m.tif"
    )
    dir.create(dirname(filename[1]), showWarnings = FALSE, recursive = TRUE)
  } else {
    stop("Invalid domain")
  }

  return(filename)
}


#' @description Function to reproject, resample, and extend a raster to a base
#' raster: 30m for CONUS and 1000m for the world.
#' @param pred The predictor raster to reproject, resample, and extend.
#' @param domain A character string of Global or USA
#' @param path A character string specifying the path to write the reprojected
#' and resampled rasters.
#' @param base The base raster to reproject and resample the predictor to.
#' @param crop A logical specifying whether to crop the raster to the L48
#' boundary BEFORE reprojecting and resampling. Only in cases where the
#' predictor is only available at the global resolution.
#' @return A print statement of the time taken to reproject and resample the
#' predictor raster.
#' @export

match_to_base <- function(pred, base, domain, path, crop) {
  # Create new filename for reprojected and resampled rasters
  filename <- get_filename(pred, domain, path)

  # Check if file exists
  if (!all(file.exists(filename))) {
    # Get method and cell size
    d_type <- terra::datatype(pred)
    method <- ifelse(d_type == "INT1U", "near", "bilinear")

    # If crop is TRUE, crop to l48 boundary
    if (crop == TRUE) {
      l48 <- terra::vect(get_l48_boundary())
      l48 <- terra::project(l48, terra::crs(pred))
      pred <- terra::crop(pred, l48, mask = TRUE)
    } else {
      pred <- pred
    }

    # For each base_raster, reproject and resample in a loop
    for (i in seq_along(base)) {
      start_time <- Sys.time() # Start time
      r <- terra::project(pred, base[[i]],
        method = method,
        wopt = list(datatype = d_type)
      )
      terra::crop(r, base[[i]],
        mask = TRUE, extend = TRUE,
        filename = filename[i], overwrite = TRUE,
        wopt = list(datatype = d_type, gdal = "COMPRESS=ZSTD")
      )
      terra::tmpFiles(orphan = TRUE, old = TRUE, remove = TRUE)
      end_time <- Sys.time() # End time
      time_taken <- (end_time - start_time)
      print(paste(
        "Time taken for", basename(filename[i]), ":", time_taken,
        attr(time_taken, "units")
      ))
    }
    msg <- print(paste0(basename(filename), " has been created"))
  } else {
    msg <- print(paste0(basename(filename), " has already been created"))
  }
  return(msg)
}

#' @description Function to get the rasters of interest based on domain
#' and type.
#' @param path A character string specifying the path to write the rasters.
#' @param domain A character string specifying the domain (Global or USA).

get_rasters <- function(path, domain = "Global") {
  if (domain == "Global") {
    rasters <- c(
      topo = get_topo_global(path),
      landcover = get_landcover_global(path),
      gdd = get_gdd_global(path),
      pop = get_pop_global(path),
      roads = get_roads_global(path),
      rails = get_rails_global(path),
      soilvars = get_soilvars_global(path),
      biovars_ = get_biovars_global(path),
      gdd = get_gdd_global(path)
    )
    names(rasters)[grep("biovars_1$", names(rasters))] <- "biovars_01"
  } else if (domain == "USA") {
    rasters <- c(
      topo = get_topo_conus(path),
      landcover = get_landcover_conus(path),
      biovars_400m_ = get_biovars_conus(path),
      biovars_ds_ = get_biovars_conus_ds(path),
      pop = get_pop_conus(path),
      biotic_evi = get_landsat_evi_conus(path),
      biotic_ndvi = get_landsat_ndvi_conus(path),
      par = get_par_conus(path),
      soilvars = get_soilvars_conus(path),
      roads = get_roads_conus(path),
      rails = get_rails_conus(path),
      gdd_ds = get_gdd_conus_ds(path),
      gdd_1km = get_gdd_global(path)
    )
    names(rasters)[grep("biovars_ds_1$", names(rasters))] <- "biovars_ds_01"
    names(rasters)[grep("biovars_400m_1$", names(rasters))] <- "biovars_400m_01"
  } else {
    stop("Invalid domain")
  }
  return(rasters)
}

#' @description Function to reproject, resample, and extend a raster to a base
#' raster: 30m for CONUS and 1000m for the world.
#' @param path A character string specifying the path to write the reprojected
#' and resampled rasters.
#' @return A print statement of the time taken to reproject and resample the
#' predictor raster.
#' @export

batch_match_to_base <- function(path) {
  # Get list of rasters
  global_rasters <- get_rasters(path, domain = "Global")

  # Get list of conus rasters
  conus_rasters <- get_rasters(path, domain = "USA")

  # Get base rasters
  base_global <- lapply(c(1000, 2500, 5000), function(x) {
    terra::rast(paste0(path, "Raster/Global/pops_sdm_base/base_", x, "m_global_epsg4326.tif"))
  })

  base_conus <- lapply(c(30, 100, 250, 500, 1000), function(x) {
    terra::rast(paste0(path, "Raster/USA/pops_sdm_base/base_", x, "m_conus_epsg4326.tif"))
  })

  # Use for loop to reproject and resample global rasters
  for (i in seq_along(global_rasters)) {
    pred <- global_rasters[[i]]
    match_to_base(pred, base = base_global, domain = "Global", path, crop = FALSE)
  }
  # Use for loop to reproject and resample conus rasters
  for (i in seq_along(conus_rasters)) {
    pred <- conus_rasters[[i]]
    match_to_base(pred, base = base_conus, domain = "USA", path, crop = FALSE)
    ## if (names(conus_rasters)[i] %in% c("gdd", "prectiming")) {
    ##   match_to_base(pred, base = base_conus, domain = "USA", path, crop = TRUE)
    ## } else {
    ##   match_to_base(pred, base = base_conus, domain = "USA", path, crop = FALSE)
    ## }
  }
}

#' @description Function to retrieve resampled rasters based on domain and res.
#' @param path A character string specifying the path to write the rasters.
#' @param domain A character string specifying the domain (Global or USA).
#' @param res A numeric vector specifying the resolutions to reproject and
#' resample the rasters to.
#' @return A list of resampled rasters.
#' @export

get_resampled_rasters <- function(path, domain, res) {
  if (domain == "Global") {
    rasters <- get_rasters(path, domain = "Global")
    files <- lapply(rasters, function(x) {
      get_filename(x, domain, path)
    })
    # Grep files based on res
    resampled_rasters <- lapply(files, function(x) {
      lapply(res, function(y) {
        terra::rast(x[grep(paste0(y, "m"), x)])
      })
    })
  } else if (domain == "USA") {
    rasters <- get_rasters(path, domain = "USA")
    files <- lapply(rasters, function(x) {
      get_filename(x, domain, path)
    })
    # Grep files based on res
    resampled_rasters <- lapply(files, function(x) {
      lapply(res, function(y) {
        terra::rast(x[grep(paste0(y, "m"), x)])
      })
    })
  } else {
    stop("Invalid domain")
  }
  return(resampled_rasters)
}

#' @description This function renames rasters returned by get_rasters
#' @param downscaled Logical. If TRUE, the function updates the names of the
#' downscaled rasters. If FALSE, the function updates the names of the
#' resampled rasters.
#' @param domain The domain provided by the user. Must be "USA" or "Global".
#' @param var_names A character vector of the variable names.

fix_raster_names <- function(downscaled, domain, var_names) {
  if ("precip" %in% var_names) {
    if (downscaled == TRUE && domain == "USA") {
      var_names <- var_names[var_names != "precip"]
      var_names <- c(var_names, paste0("biovars_ds_", 13:19))
      print("For precip, selecting downscaled rasters")
    } else if (downscaled == FALSE && domain == "USA") {
      var_names <- var_names[var_names != "precip"]
      var_names <- c(var_names, paste0("biovars_400m_", 13:19))
      print("For precip, selecting resampled rasters")
    } else {
      var_names <- var_names[var_names != "precip"]
      var_names <- c(var_names, paste0("biovars_", 13:19))
    }
  }
  if ("temp" %in% var_names) {
    if (downscaled == TRUE && domain == "USA") {
      var_names <- var_names[var_names != "temp"]
      var_names <- c(var_names, paste0("biovars_ds_", c("01", 2:12)))
      print("For temp, selecting downscaled rasters")
    } else if (downscaled == FALSE && domain == "USA") {
      var_names <- var_names[var_names != "temp"]
      var_names <- c(var_names, paste0("biovars_400m_", c("01", 2:12)))
      print("For temp, selecting resampled rasters")
    } else {
      var_names <- var_names[var_names != "temp"]
      var_names <- c(var_names, paste0("biovars_", c("01", 2:12)))
    }
  }
  return(var_names)
}

#' @description This function allows the users to reduce the set of predictors
#' returned by the function get_rasters. The function returns the reduced set.
#' @param domain The domain of interest. Options are "USA" or "Global".
#' @param path The path to the Data folder.
#' @param temp If TRUE, bioclimatic variables: 1-12 are included.
#' @param precip If TRUE, bioclimate variables: 13-19 are included.
#' @param topo If TRUE, topographic variables are included.
#' @param land If TRUE, land cover variables are included.
#' @param soils If TRUE, soil variables are included.
#' @param pop If TRUE, population variables are included.
#' @param gdd If TRUE, growing degree days variables are included.
#' @param biotic If TRUE, biotic variables (NDVI & EVI) are included. Not
#' applicable to Global domain.
#' @param light If TRUE, light variables (PAR) are included. Not applicable
#' to Global domain.
#' @param rails If TRUE, railroads variables are included.
#' @param roads If TRUE, roads variables are included.
#' @param downscaled If TRUE, only use downscaled bioclimatic
#' predictors returned. Only applies to USA domain.

subset_rasters <- function(
    path,
    domain,
    all = TRUE,
    temp = TRUE,
    precip = TRUE,
    topo = TRUE,
    land = TRUE,
    soils = TRUE,
    pop = TRUE,
    gdd = TRUE,
    biotic = TRUE,
    light = TRUE,
    rails = TRUE,
    roads = TRUE,
    downscaled = TRUE) {
  # Get rasters
  rasters <- get_rasters(path, domain)

  # If all is FALSE, only keep TRUE arguments
  if (all == FALSE) {
    # Subset rasters based on user input
    if (domain == "Global") {
      vars <- c(
        temp, precip, topo, land, soils, pop, gdd, rails, roads
      )
      var_names <- c(
        "temp", "precip", "topo", "landcover", "soilvars", "pop", "gdd",
        "rails", "roads"
      )
    } else if (domain == "USA") {
      vars <- c(
        temp, precip, topo, land, soils, pop, gdd, biotic,
        light, rails, roads
      )
      var_names <- c(
        "temp", "precip", "topo", "landcover", "soilvars", "pop", "gdd",
        "biotic", "par", "rails", "roads"
      )
    }
    # Only keep TRUE arguments
    var_names <- var_names[vars]
    print(paste("Selecting predictors:", paste(var_names, collapse = ", ")))
    var_names <- fix_raster_names(downscaled, domain, var_names)
    rasters <- rasters[grep(paste(var_names, collapse = "|"), names(rasters))]
  } else {
    print("Selecting all predictors")
    if (downscaled == FALSE) {
      rasters <- rasters[grep("_ds_", names(rasters), invert = TRUE)]
      print("For temp & precip, selecting downscaled rasters")
    } else {
      rasters <- rasters[grep("_400m_", names(rasters), invert = TRUE)]
      print("For temp & precip, selecting resampled rasters")
    }
  }
  return(rasters)
}

#####################################################################################
######
###### Formerly a separate script called get_part_sblock.R
######

#' pre_tr_te
#'
#' @noRd
pre_tr_te <- function(data, p_names, h) {
    train <- list()
    test <- list()

    if (any(c("train", "train-test", "test")
    %in%
        unique(data[, p_names[h]]))) {
        np2 <- 1

        filt <- grepl("train", data[, p_names[h]])
        train[[1]] <- data[filt, ] %>%
            dplyr::select(-p_names[!p_names == p_names[h]])

        filt <- grepl("test", data[, p_names[h]])
        test[[1]] <- data[filt, ] %>%
            dplyr::select(-p_names[!p_names == p_names[h]])
    } else {
        np2 <- max(data[p_names[h]])

        for (i in 1:np2) {
            train[[i]] <- data[data[p_names[h]] != i, ] %>%
                dplyr::select(-p_names[!p_names == p_names[h]])

            test[[i]] <- data[data[p_names[h]] == i, ] %>%
                dplyr::select(-p_names[!p_names == p_names[h]])
        }
    }
    return(list(train = train, test = test, np2 = np2))
}



# Inverse bioclim
#'
#' @noRd
#'
bio <- function(data, env_layer) {
    . <- NULL
    if (class(data)[1] != "data.frame") {
        data <- data.frame(data)
    }
    if (!methods::is(env_layer, "SpatRaster")) {
        env_layer <- terra::rast(env_layer)
    }

    data <- na.omit(data)

    result <- env_layer[[1]]
    result[] <- NA

    minv <- apply(data, 2, min)
    maxv <- apply(data, 2, max)
    vnames <- names(data)

    data_2 <- data %>%
        na.omit() %>%
        apply(., 2, sort) %>%
        data.frame()

    rnk <- function(x, y) {
        b <- apply(y, 1, FUN = function(z) sum(x < z))
        t <- apply(y, 1, FUN = function(z) sum(x == z))
        r <- (b + 0.5 * t) / length(x)
        i <- which(r > 0.5)
        r[i] <- 1 - r[i]
        r * 2
    }

    var_df <- terra::as.data.frame(env_layer)
    var_df <- na.omit(var_df)

    k <- (apply(t(var_df) >= minv, 2, all) &
        apply(t(var_df) <= maxv, 2, all))

    for (j in vnames) {
        var_df[k, j] <- rnk(
            data_2[, j],
            var_df[k, j, drop = FALSE]
        )
    }
    var_df[!k, ] <- 0
    res <- apply(var_df, 1, min)
    result[as.numeric(names(res))] <- res
    return(result)
}

inv_bio <- function(e, p) {
    if (!methods::is(e, "SpatRaster")) {
        e <- terra::rast(e)
    }
    r <- bio(data = terra::extract(e, p)[-1], env_layer = e)
    r <- (r - terra::minmax(r)[1]) /
        (terra::minmax(r)[2] - terra::minmax(r)[1])
    r <- r <= 0.01 # environmental constrain
    r[which(r[, ] == FALSE)] <- NA
    return(r)
}


#' Inverse geo
#'
#' @noRd
#'
inv_geo <- function(e, p, d) {
    colnames(p) <- c("x", "y")
    p <- terra::vect(p, geom = c("x", "y"), crs = terra::crs(e))
    b <- terra::buffer(p, width = d)
    b <- terra::rasterize(b, e, background = 0)
    e <- terra::mask(e, b, maskvalues = 1)
    return(e)
}

#' Boyce
#'
#' @description This function calculate Boyce index performance metric. Codes were adapted from
#' enmSdm package.
#'
#' @noRd
boyce <- function(pres,
                  contrast,
                  n_bins = 101,
                  n_width = 0.1) {
    lowest <- min(c(pres, contrast), na.rm = TRUE)
    highest <- max(c(pres, contrast), na.rm = TRUE) + .Machine$double.eps
    window_width <- n_width * (highest - lowest)

    lows <- seq(lowest, highest - window_width, length.out = n_bins)
    highs <- seq(lowest + window_width + .Machine$double.eps, highest, length.out = n_bins)

    ## initiate variables to store predicted/expected (P/E) values
    freq_pres <- NA
    freq_contrast <- NA

    # tally proportion of test presences/background in each class
    for (i in 1:n_bins) {
        # number of presence predictions in a class
        freq_pres[i] <-
            sum(pres >= lows[i] & pres < highs[i], na.rm = TRUE)

        # number of background predictions in this class
        freq_contrast[i] <-
            sum(contrast >= lows[i] & contrast < highs[i], na.rm = TRUE)
    }

    # mean bin prediction
    mean_pred <- rowMeans(cbind(lows, highs))

    # add small number to each bin that has 0 background frequency but does have a presence frequency > 0
    if (any(freq_pres > 0 & freq_contrast == 0)) {
        small_value <- 0.5
        freq_contrast[freq_pres > 0 & freq_contrast == 0] <- small_value
    }

    # remove classes with 0 presence frequency
    if (any(freq_pres == 0)) {
        zeros <- which(freq_pres == 0)
        mean_pred[zeros] <- NA
        freq_pres[zeros] <- NA
        freq_contrast[zeros] <- NA
    }

    # remove classes with 0 background frequency
    if (any(0 %in% freq_contrast)) {
        zeros <- which(freq_pres == 0)
        mean_pred[zeros] <- NA
        freq_pres[zeros] <- NA
        freq_contrast[zeros] <- NA
    }

    P <- freq_pres / length(pres)
    E <- freq_contrast / length(contrast)
    PE <- P / E

    # remove NAs
    rm_nas <- stats::complete.cases(data.frame(mean_pred, PE))
    # mean_pred <- mean_pred[rm_nas]
    # PE <- PE[rm_nas]

    # calculate Boyce index
    result <- stats::cor(
        x = ifelse(is.na(mean_pred), 0, mean_pred),
        y = ifelse(is.na(PE), 0, PE), method = "spearman"
    )
    return(result)
}

# maxnet:::predict.maxnet()
#' Predict maxnet
#' @importFrom stats model.matrix
#' @noRd
predict_maxnet <- function(object, newdata, clamp = TRUE, type = c("link", "exponential", "cloglog", "logistic"), ...) {
    categoricalval <- function(x, category) {
        ifelse(x == category, 1, 0)
    }
    thresholdval <- function(x, knot) {
        ifelse(x >= knot, 1, 0)
    }
    hingeval <- function(x, min, max) {
        pmin(1, pmax(0, (x - min) / (max - min)))
    }

    if (clamp) {
        for (v in intersect(names(object$varmax), names(newdata))) {
            newdata[, v] <- pmin(
                pmax(newdata[, v], object$varmin[v]),
                object$varmax[v]
            )
        }
    }
    terms <- sub(
        "hinge\\((.*)\\):(.*):(.*)$", "hingeval(\\1,\\2,\\3)",
        names(object$betas)
    )
    terms <- sub(
        "categorical\\((.*)\\):(.*)$", "categoricalval(\\1,\"\\2\")",
        terms
    )
    terms <- sub(
        "thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)",
        terms
    )
    f <- formula(paste("~", paste(terms, collapse = " + "), "-1"))
    mm <- model.matrix(f, data.frame(newdata))
    if (clamp) {
        mm <- t(pmin(
            pmax(t(mm), object$featuremins[names(object$betas)]),
            object$featuremaxs[names(object$betas)]
        ))
    }
    link <- (mm %*% object$betas) + object$alpha
    type <- match.arg(type)

    if (type == "link") {
        return(link)
    }
    if (type == "exponential") {
        return(exp(link))
    }
    if (type == "cloglog") {
        return(1 - exp(0 - exp(object$entropy + link)))
    }
    if (type == "logistic") {
        return(1 / (1 + exp(-object$entropy - link)))
    }
}

#' Outliers with Reverse Jackknife
#'
#' @noRd
#'
rev_jack <- function(v) {
    v2 <- v
    v <- unique(v)
    lgh <- length(v) - 1
    t1 <- (0.95 * sqrt(length(v))) + 0.2
    x <- sort(v)
    y <- rep(0, lgh)
    for (i in seq_len(lgh)) {
        x1 <- x[i + 1]
        if (x[i] < mean(v)) {
            y[i] <- (x1 - x[i]) * (mean(v) - x[i])
        } else {
            y[i] <- (x1 - x[i]) * (x1 - mean(v))
        }
    }
    my <- mean(y)
    z <- y / (sqrt(sum((y - my)^2) / lgh))
    out <- rep(0, length(v2))
    if (any(z > t1)) {
        f <- which(z > t1)
        v <- x[f]
        if (v < median(x)) {
            xa <- (v2 <= v) * 1
            out <- out + xa
        }
        if (v > median(x)) {
            xb <- (v2 >= v) * 1
            out <- out + xb
        }
    } else {
        out <- out
    }
    return(which(out == 1))
}

#' Calculate amount of data for each training dataset in a given partition
#'
#' @noRd
#'
n_training <- function(data, partition) {
    . <- partt <- NULL
    if (any(c("train", "train-test", "test")
    %in%
        (data %>%
            dplyr::select(dplyr::starts_with({{ partition }})) %>%
            dplyr::pull() %>%
            unique()))) {
        nn_part <- data %>%
            dplyr::select(dplyr::starts_with({{ partition }})) %>%
            apply(., 2, table) %>%
            data.frame()
        nn_part <- nn_part %>% dplyr::mutate(partt = rownames(nn_part))
        nn_part$partt[grepl("train", nn_part$partt)] <- "train"
        nn_part <- nn_part %>%
            dplyr::filter(partt == "train") %>%
            dplyr::select(-partt)
        nn_part <- colSums(nn_part)
    } else {
        data <- data %>%
            dplyr::select(dplyr::starts_with({{ partition }}))

        nn_part <- list()
        for (ppp in 1:ncol(data)) {
            nn_part[[ppp]] <- data %>%
                dplyr::pull(ppp) %>%
                table() %>%
                c()
            sm <- nn_part[[ppp]] %>% sum()
            nn_part[[ppp]] <- sapply(nn_part[[ppp]], function(x) sum(sm - x))
        }
        nn_part <- unlist(nn_part)
    }
    return(nn_part)
}

#' Calculate number of coefficient for gam models
#'
#' @noRd
#'
n_coefficients <- function(data, predictors, predictors_f = NULL, k = 10) {
    data <- data.frame(data)
    if (k < 0) {
        k <- 10
    }
    if (!is.null(predictors_f)) {
        n_levels <- rep(NA, length(predictors_f))
        for (fff in 1:length(predictors_f)) {
            n_levels[fff] <- unique(data[, predictors_f]) %>%
                na.omit() %>%
                length()
        }
        n_levels <- sum(n_levels)
    } else {
        n_levels <- 0
    }
    n <- (k - 1) * length(predictors) + n_levels
    return(n)
}

#' Euclidean distance for extrapolation
#'
#' @noRd
#'
euc_dist <- function(x, y) {
    if (!methods::is(x, "matrix")) {
        x <- as.matrix(x)
    }
    if (!methods::is(y, "matrix")) {
        y <- as.matrix(y)
    }
    result <- matrix(0, nrow = nrow(x), ncol = nrow(y))
    for (ii in 1:nrow(y)) {
        result[, ii] <- sqrt(colSums((t(x) - y[ii, ])^2))
    }
    rownames(result) <- rownames(x)
    colnames(result) <- rownames(y)
    return(result)
}

#' Moran I, based on ape package
#'
#' @noRd
#'
morani <- function(x, weight, na.rm = FALSE, scaled = TRUE) {
    if (dim(weight)[1] != dim(weight)[2]) {
        stop("'weight' must be a square matrix")
    }
    n <- length(x)
    if (dim(weight)[1] != n) {
        stop("'weight' must have as many rows as observations in 'x'")
    }
    ei <- -1 / (n - 1)
    nas <- is.na(x)
    if (any(nas)) {
        if (na.rm) {
            x <- x[!nas]
            n <- length(x)
            weight <- weight[!nas, !nas]
        } else {
            warning("'x' has missing values: maybe you wanted to set na.rm = TRUE?")
            return(NA)
        }
    }
    rs <- rowSums(weight)
    rs[rs == 0] <- 1
    weight <- weight / rs
    s <- sum(weight)
    m <- mean(x)
    y <- x - m
    cv <- sum(weight * y %o% y)
    v <- sum(y^2)
    res <- (n / s) * (cv / v)
    if (scaled) {
        imax <- (n / s) * (sd(rowSums(weight) * y) / sqrt(v / (n - 1)))
        res <- res / imax
    }

    return(res)
}

#' Mahalanobis distance
#'
#' @noRd
#'
mah_dist <- function(x, y, cov) {
    if (!methods::is(x, "matrix")) {
        x <- as.matrix(x)
    }
    if (!methods::is(y, "matrix")) {
        y <- as.matrix(y)
    }
    result <- matrix(0, nrow = nrow(x), ncol = nrow(y))
    for (ii in 1:nrow(y)) {
        # root square of squared Mahalanobis distance
        result[, ii] <- sqrt(mahalanobis(x = x, center = y[ii, ], cov = cov))
    }
    rownames(result) <- rownames(x)
    colnames(result) <- rownames(y)
    return(result)
}

#' @description Debugged version of the function 'spart_block' from the package
#' flexsdm. This function is used to partition the data into blocks based
#' on the spatial autocorrelation, and environmental similarity.

part_sblock_upd <- function(
    env_layer, data, x, y, pr_ab, n_part = 3, min_res_mult = 3,
    max_res_mult = 200, num_grids = 30, min_occ = 10, prop = 0.5) {
    data <- dplyr::tibble(data)
    data <- data[, c(pr_ab, x, y)]
    colnames(data) <- c("pr_ab", "x", "y")
    if (any(!unique(data[, "pr_ab"][[1]]) %in% c(0, 1))) {
        stop(
            "values in pr_ab column did not match with 0 and 1:\nunique list values in pr_ab column are: ",
            paste(unique(data[, "pr_ab"]), collapse = " ")
        )
    }
    data <- dplyr::tibble(data, terra::extract(env_layer, data[, 2:3])[-1])
    filt <- stats::complete.cases(data)
    if (sum(!filt) > 0) {
        data <- data[filt, ]
        message(sum(!filt), " rows were excluded from database because NAs were found")
    }
    rm(filt)
    pa <- data %>% dplyr::pull(pr_ab)
    cell_size <- seq(terra::res(env_layer[[1]])[1] * min_res_mult,
        terra::res(env_layer[[1]])[1] * max_res_mult,
        length.out = num_grids
    )
    message(
        "The following grid cell sizes will be tested:\n",
        paste(round(cell_size, 2), collapse = " | "), "\n"
    )
    message("Creating basic raster mask...\n")
    mask <- env_layer[[1]]
    names(mask) <- "group"
    mask[!is.na(mask)] <- 1
    e <- terra::ext(mask)
    message("Searching for the optimal grid size...\n")
    mask2 <- mask
    mask2[] <- 0
    presences2 <- data
    presences2 <- terra::vect(presences2,
        geom = c("x", "y"),
        crs = terra::crs(mask)
    )
    grid <- list()
    DIM <- matrix(0, length(cell_size), 2)
    colnames(DIM) <- c("R", "C")
    for (i in 1:length(cell_size)) {
        mask3 <- mask2
        terra::res(mask3) <- cell_size[i]
        mask3 <- terra::extend(mask3, y = c(1, 1))
        DIM[i, ] <- dim(mask3)[1:2]
        terra::values(mask3) <- 1
        NAS <- c(terra::extract(mask3, presences2)[-1])
        if (any(is.na(NAS))) {
            while (any(is.na(NAS))) {
                terra::ext(mask3) <- terra::ext(mask3) + cell_size[i]
                terra::res(mask3) <- cell_size[i]
                DIM[i, ] <- dim(mask3)[1:2]
                terra::values(mask3) <- 1
                NAS <- terra::extract(mask3, presences2)[-1]
            }
        }
        grid[[i]] <- mask3
    }
    rm(list = c("mask3", "mask2", "mask"))
    for (i in 1:length(grid)) {
        if (n_part %% 2 == 0) {
            group <- c(
                rep(1:n_part, DIM[i, 2])[1:DIM[i, 2]],
                rep(c((n_part / 2 + 1):n_part, 1:(n_part / 2)), DIM[
                    i,
                    2
                ])[1:DIM[i, 2]]
            )
        }
        if (n_part %% 2 == 1) {
            group <- c(
                rep(1:n_part, DIM[i, 2])[1:DIM[i, 2]],
                rep(
                    c(as.integer(n_part / 2 + 1):n_part, 1:(n_part / 2)),
                    DIM[i, 2]
                )[1:DIM[i, 2]]
            )
        }
        terra::values(grid[[i]]) <- rep(group, length.out = terra::ncell(grid[[i]]))
    }
    part <- data.frame(matrix(0, nrow(presences2), length(grid)))
    for (i in 1:length(grid)) {
        part[, i] <- terra::extract(grid[[i]], presences2)[
            ,
            2
        ]
    }
    part <- dplyr::tibble(part)
    pp <- sapply(part[pa == 1, ], function(x) {
        length(unique(x))
    })
    pp <- pp == n_part
    pf <- sapply(part[pa == 1, ], table)
    if (is.list(pf) == TRUE) {
        pf <- which(sapply(pf, min) < min_occ)
    } else {
        pf <- which(apply(pf, 2, min) < min_occ)
    }
    pp[pf] <- FALSE
    cell_size <- cell_size[pp]
    grid <- grid[pp]
    part <- part[, pp]
    names(part) <- names(which(pp == TRUE))
    if (any(unique(pa) == 0)) {
        pa <- presences2$pr_ab
        pp <- sapply(part[pa == 0, ], function(x) {
            length(unique(x))
        })
        pp <- pp == n_part
        pf <- sapply(part[pa == 0, ], table)
        if (is.list(pf) == TRUE) {
            pf <- which(sapply(pf, min) < min_occ)
        } else {
            pf <- which(apply(pf, 2, min) < min_occ)
        }
        pp[pf] <- FALSE
        cell_size <- cell_size[pp]
        grid <- grid[pp]
        part <- part[, pp]
        names(part) <- names(which(pp == TRUE))
    }
    if (ncol(part) == 0) {
        message("It was not possible to find a good partition. Try to change values in 'n_part', or in 'min_res_mult', 'max_res_mult', or 'num_grids'")
        return(NA)
    }
    ncell <- data.frame(matrix(0, nrow(presences2), length(grid)))
    for (i in 1:length(grid)) {
        ncell[, i] <- terra::cellFromXY(grid[[i]], terra::geom(presences2)[
            ,
            c("x", "y")
        ])
    }
    sd_p <- rep(NA, length(grid))
    if (any(unique(pa) == 0)) {
        sd_a <- rep(NA, length(grid))
    }
    for (i in 1:ncol(part)) {
        if (any(unique(pa) == 0)) {
            sd_a[i] <- stats::sd(table(part[pa == 0, i]))
        }
        sd_p[i] <- stats::sd(table(part[pa == 1, i]))
    }
    Env.P <- terra::extract(env_layer, presences2)[-1]
    env_sim <- rep(NA, length(grid))
    for (i in 1:ncol(part)) {
        cmb <- unique(part[, i][[1]]) %>% combn(2)
        Env.P1 <- cbind(part[i], Env.P)
        Env.P1 <- Env.P1[complete.cases(Env.P1), ]
        Env.P1 <- split(Env.P1[, -1], Env.P1[, 1])
        euq_c <- list()
        for (r in 1:ncol(cmb)) {
            euq_c[[r]] <- euc_dist(Env.P1[[cmb[1, r]]], Env.P1[[cmb[2, r]]]) %>% mean()
        }
        env_sim[i] <- euq_c %>%
            unlist() %>%
            mean()
        rm(list = c("Env.P1"))
    }
    spa_auto <- rep(NA, length(grid))
    presences2 <- terra::geom(presences2)[, c("x", "y")] %>%
        as.data.frame()
    dist <- euc_dist(presences2, presences2)
    dist <- 1 / dist
    diag(dist) <- 0
    dist[which(dist == Inf)] <- 0
    for (p in 1:ncol(part)) {
        cmb <- unique(part[, p][[1]]) %>% utils::combn(2)
        imoran_grid_c <- rep(NA, ncol(cmb))
        dff <- dplyr::tibble(
            nrow = 1:nrow(part), data["pr_ab"],
            group = part[p][[1]]
        )
        for (c in 1:ncol(cmb)) {
            filt <- dff %>%
                dplyr::group_by(group, pr_ab) %>%
                dplyr::slice_sample(prop = prop) %>%
                dplyr::pull(nrow) %>%
                sort()
            odd <- which((part[p][[1]] == cmb[1, c])[filt])
            even <- which((part[p][[1]] == cmb[2, c])[filt])
            dist2 <- dist[filt, filt]
            dist2[odd, odd] <- 0
            dist2[even, even] <- 0
            mins <- apply(dist2, 2, function(x) {
                max(x, na.rm = TRUE)
            })
            for (i in 1:length(mins)) {
                dist2[, i] <- ifelse(dist2[, i] == mins[i], mins[i],
                    0
                )
            }
            if (nrow(data) < 3) {
                imoran_grid_c[c] <- NA
            } else {
                im <- sapply(data[filt, names(env_layer)], function(x) {
                    suppressMessages(morani(x, dist2,
                        na.rm = TRUE,
                        scaled = TRUE
                    ))
                })
                imoran_grid_c[c] <- mean(abs(im))
            }
        }
        spa_auto[p] <- mean(imoran_grid_c)
    }
    Opt <- if (any(unique(pa) == 0)) {
        data.frame(
            n_grid = 1:length(cell_size), cell_size = cell_size,
            round(
                data.frame(spa_auto, env_sim, sd_p, sd_a),
                3
            )
        )
    } else {
        data.frame(
            n_grid = 1:length(cell_size), cell_size = cell_size,
            round(data.frame(spa_auto, env_sim, sd_p), 3)
        )
    }
    Opt2 <- Opt
    rownames(Opt2) <- colnames(part)
    Dup <- if (any(unique(pa) == 0)) {
        !duplicated(Opt2[c("spa_auto", "env_sim", "sd_p", "sd_a")])
    } else {
        !duplicated(Opt2[c("spa_auto", "env_sim", "sd_p")])
    }
    Opt2 <- Opt2[Dup, ]
    while (nrow(Opt2) > 1) {
        if (nrow(Opt2) == 1) {
            break
        }
        Opt2 <- Opt2[which(Opt2$spa_auto <= summary(Opt2$spa_auto)[2]), ]
        if (nrow(Opt2) == 1) {
            break
        }
        Opt2 <- Opt2[which(Opt2$env_sim >= summary(Opt2$env_sim)[5]), ]
        if (nrow(Opt2) == 1) {
            break
        }
        Opt2 <- Opt2[which(Opt2$sd_p <= summary(Opt2$sd_p)[2]), ]
        if (nrow(Opt2) == 2) {
            break
        }
        if (any(unique(pa) == 0)) {
            Opt2 <- Opt2[which(Opt2$sd_a <= summary(Opt2$sd_a)[2]), ]
            if (nrow(Opt2) == 2) {
                break
            }
        }
        if (unique(Opt2$spa_auto) && unique(Opt2$env_sim) &&
            unique(Opt2$sd_p)) {
            Opt2 <- Opt2[nrow(Opt2), ]
        }
    }
    if (nrow(Opt2) > 1) {
        Opt2 <- Opt2[nrow(Opt2), ]
    }
    result <- data.frame(data, .part = c(part[, rownames(Opt2)])[[1]])
    result <- result %>% dplyr::select(-names(env_layer))
    colnames(result) <- c("pr_ab", "x", "y", ".part")
    result <- result[c("x", "y", "pr_ab", ".part")]
    grid <- grid[[Opt2$n_grid]]
    names(grid) <- ".part"
    out <- list(
        part = dplyr::tibble(result), best_part_info = dplyr::tibble(Opt2),
        grid = grid
    )
    return(out)
}

#####################################################################################
######
###### Formerly a separate script called raster_base.R
######

#' Function to create a 1 arc second (30m) base raster for USA;
#' and aggregate to common resolutions: 30m, 100m, 250m, 500m, 1000m
#' @param path Path to write the data to
#' @return A print message indicating that the base rasters have been created
#' @export

base_conus <- function(path) {
  # Setup base raster filenames for check
  proj <- c(rep("epsg4326", 5))
  common_res_names <- c(30, 100, 250, 500, 1000)
  fn <- paste0(
    path, "Raster/USA/pops_sdm_base/base_", common_res_names, "m_conus_",
    proj, ".tif"
  )
  # Check if all fn do not exist, create base rasters
  if (!all(file.exists(fn))) {
    # Create 30m base raster for USA from 1s L48 boundary and 1s elevation
    l48 <- terra::vect(get_l48_boundary())
    elev <- get_topo_conus(path)[[1]]
    elev <- terra::clamp(elev, 1, 1, datatype = "INT1U")
    l48 <- terra::project(l48, y = terra::crs(elev))
    dir.create(paste0(path, "Raster/USA/pops_sdm_base/"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    elev <- terra::crop(elev, l48,
      mask = TRUE, filename =
        paste0(path, "Raster/USA/pops_sdm_base/base_30m_conus_epsg4326.tif"), # nolint line length
      overwrite = TRUE, wopt = list(
        gdal = "COMPRESS=ZSTD",
        datatype = "INT1U"
      )
    )

    # Base rasters at common resolutions: 30m, 100m, 250m, 500m, 1000m EPSG:4326
    base_res <- terra::res(elev)[1]
    common_res_names <- common_res_names[-1]
    common_res_vals <- c(
      base_res * 3, base_res * 7.5, base_res * 15,
      base_res * 30
    )
    lapply(seq_along(common_res_vals), function(x) {
      base_filename <- paste0(
        path, "Raster/USA/pops_sdm_base/base_",
        common_res_names[x], "m_conus_epsg4326.tif"
      )
      terra::project(elev,
        y = "EPSG:4326", threads = TRUE, method = "near",
        filename = base_filename, res = common_res_vals[x],
        overwrite = TRUE, wopt = list(
          gdal = "COMPRESS=ZSTD",
          datatype = "INT1U"
        )
      )
    })
    msg <- paste0(
      "Base rasters created for USA at resolutions: ",
      paste0(common_res_names, collapse = ", "), "m"
    )
  } else {
    msg <- paste0(
      "Base rasters already exist for USA at resolutions: ",
      paste0(common_res_names, collapse = ", "), "m"
    )
  }
  base <- lapply(fn, terra::rast)
  return(base)
}

#' Function to create a 30 arc second (1km) base raster for the world;
#' and aggregate to common resolutions: 2.5km, 5km, 10km.
#' @param path Path to write the data to
#' @return A print message indicating that the base rasters have been created.
#' @export

base_global <- function(path) {
  # Setup base raster filenames for check
  proj <- c(rep("epsg4326", 5))
  common_res_names <- c(1000, 2500, 5000)
  fn <- paste0(
    path, "Raster/Global/pops_sdm_base/base_", common_res_names, "m_conus_",
    proj, ".tif"
  )
  # Check if all fn do not exist, create base rasters
  if (!all(file.exists(fn))) {
    # Create 1000m base raster for globe from wc2.1_30s elevation
    elev <- get_topo_global(path)[[1]]
    elev <- terra::clamp(elev, 1, 1, datatype = "INT1U")
    dir.create(paste0(path, "Raster/Global/pops_sdm_base/"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    # Base rasters at common resolutions: 30m, 100m, 250m, 500m, 1000m EPSG:4326
    base_res <- terra::res(elev)[1]
    common_res_vals <- c(base_res, base_res * 2.5, base_res * 5)
    common_res_names <- common_res_names
    lapply(seq_along(common_res_vals), function(x) {
      base_filename <- paste0(
        path, "Raster/Global/pops_sdm_base/base_",
        common_res_names[x], "m_global_epsg4326.tif"
      )
      terra::project(elev,
        y = "EPSG:4326", threads = TRUE, method = "near",
        filename = base_filename, res = common_res_vals[x],
        overwrite = TRUE, wopt = list(
          gdal = "COMPRESS=ZSTD",
          datatype = "INT1U"
        )
      )
    })
    msg <- paste0(
      "Base rasters created for the Globe at resolutions: ",
      paste0(common_res_names, collapse = ", "), "m"
    )
  } else {
    msg <- paste0(
      "Base rasters already exist for the Globe at resolutions: ",
      paste0(common_res_names, collapse = ", "), "m"
    )
  }
  base <- lapply(fn, terra::rast)
  return(base)
}

#####################################################################################
######
###### Formerly a separate script called sample_background.R
######

sample_bg_upd <- function (data, x, y, n, method = "random", rlayer, maskval = NULL, 
    calibarea = NULL, rbias = NULL, sp_name = NULL) 
{
    . <- NULL
    if (!method[1] %in% c("random", "thickening", "biased")) {
        stop("argument 'method' was misused, available methods 'random', 'thickening'")
    }
    if (method[1] %in% c("biased") & is.null(rbias)) {
        stop("for using 'random' method a raster layer with biases data must be provided in 'rbias' argument")
    }
    if (class(rlayer)[1] != "SpatRaster") {
        rlayer <- terra::rast(rlayer)
    }
    if (!is.null(rbias)) {
        if (class(rbias)[1] != "SpatRaster") 
            rbias <- terra::rast(rbias)
    }
    rlayer <- rlayer[[1]]
    data <- data[, c(x, y)]
    rlayer[na.omit(terra::cellFromXY(rlayer, as.matrix(data)))] <- NA
    if (!is.null(calibarea)) {
        rlayer <- rlayer %>% terra::crop(., calibarea) %>% terra::mask(., 
            calibarea)
    }
    if (!is.null(maskval)) {
        if (is.factor(maskval)) {
            maskval <- which(levels(maskval) %in% as.character(maskval))
            rlayer <- rlayer * 1
        }
        filt <- terra::match(rlayer, maskval)
        rlayer <- terra::mask(rlayer, filt)
    }
    if (method[1] %in% c("biased")) {
        if (any(!(ext(rlayer)[1:4] %in% ext(rbias)[1:4])) | all(!res(rlayer) %in% 
            res(rbias))) {
            if (!all(res(rlayer) %in% res(rbias))) {
                rbias <- terra::resample(rbias, rlayer, method = "bilinear")
            }
            rbias2 <- rbias %>% terra::crop(., rlayer) %>% terra::mask(., 
                rlayer)
            rbias <- rbias2
            rm(rbias2)
        }
    }
    if (method[1] %in% c("biased")) {
        rlayer <- mask(rlayer, rbias)
    }
    ncellr <- terra::global(!is.na(rlayer), sum)
    if (any(method == "thickening")) {
        data2 <- terra::vect(data, geom = c(x, y), crs = crs(rlayer))
        if (is.na(method["width"])) {
            buf_with <- mean(terra::distance(data2))
        } else {
            buf_with <- as.numeric(method["width"])
        }
        buf <- terra::buffer(data2, buf_with, quadsegs = 10)
        buf_r <- terra::rasterize(buf, rlayer, background = 0, fun = "count")
        buf_r <- terra::mask(buf_r, rlayer)
    }
    if (ncellr < n) {
        message("Number of background-points exceeds number of cell will be returned ", 
            ncellr, " background-points")
        cell_samp <- terra::as.data.frame(rlayer, na.rm = TRUE, 
            cells = TRUE)[, "cell"]
        cell_samp <- terra::xyFromCell(rlayer, cell_samp) %>% 
            data.frame() %>% dplyr::tibble()
    } else {
        cell_samp <- terra::as.data.frame(rlayer, na.rm = TRUE, 
            cells = TRUE)[, "cell"]
        if (any(method == "random")) {
            cell_samp <- sample(cell_samp, size = n, replace = FALSE, 
                prob = NULL)
        } else if (any(method == "thickening")) {
            cell_samp <- sample(cell_samp, size = n, replace = FALSE, 
                prob = buf_r[cell_samp][, 1])
        } else if (any(method == "biased")) {
            cell_samp <- sample(cell_samp, size = n, replace = FALSE, 
                prob = rbias[cell_samp][, 1])
        }
        cell_samp <- terra::xyFromCell(rlayer, cell_samp) %>% 
            data.frame() %>% dplyr::tibble()
    }
    colnames(cell_samp) <- c(x, y)
    cell_samp$pr_ab <- 0
    if (!is.null(sp_name)) {
        cell_samp <- tibble(sp = sp_name, cell_samp)
    }
    return(cell_samp)
}


