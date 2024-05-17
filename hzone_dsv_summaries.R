library(data.table)
library(terra)
library(tidyverse)

# List all of the files for the years and hardiness zones
lh1 <- list.files("Z:/Late_blight/lbdsveconu/", pattern = "zone", full.names = T, recursive = T)
lh1 <- lh1[20:247]
lh1 <- gtools::mixedsort(lh1)

# Break up the big list to zone
# Number of names in each sublist
names_per_sublist <- 12

# Create a sequence of group numbers for each file name
groups <- ceiling(seq_along(lh1) / names_per_sublist)

# Split the file names into sublists based on the group numbers
sublists <- split(lh1, groups)

# Function to process a sublist's worth of rasters
yearlyval <- function(subind, dayind) {
  # Raster
  raszone = rast(sublists[[subind]])
  
  # Sequence for a given day of year, e.g 1, 366, etc.
  dayseq = seq(dayind, nlyr(raszone), by = 365)
  
  # Values for the layers that match 
  dayvalue = values(raszone[[dayseq]], na.rm=T, dataframe =T)

  # Add colnames
  colnames(dayvalue) = paste0("y",2009:2020)

  # Add day of year to track across the collected rasters
  dayvalue$day = dayind
  
  return(dayvalue)
}

# Apply the function to all rasters and days of years
yearcol_list <- lapply(1:length(sublists), function(subind) {
  yearcol_sub <- lapply(1:365, function(dayind) yearlyval(subind, dayind))
  collection <- do.call(rbind, yearcol_sub)
  write.csv(collection, paste0("Z:/Late_blight/lbdsveconu/appscores/csv/dsv/sub_",dayind,".csv), row.names=F)
})

