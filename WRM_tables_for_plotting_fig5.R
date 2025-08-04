# The tables created here are used to calculate the values that the bar charts
# in Fig. 5 rely on. They take daily WRM rasters and create matrices. Then the 
# lists of matrices are divided by growing season and region. 

library(data.table)
library(terra)
library(tidyverse)

# LB DSV raster file list
lh1 <- list.files("Z:/Late_blight/lbdsveconu/", pattern = "zone", full.names = T, recursive = T)
lh1 <- lh1[20:247]
lh1 <- gtools::mixedsort(lh1)

# Divide lh1 into lists based on hardiness zone and year
# Years
names_per_sublist <- 12

# Create a sequence of group numbers for each file name
groups <- ceiling(seq_along(lh1) / names_per_sublist)

# Split the file names into sublists based on the group numbers
sublists <- split(lh1, groups)

# Move through the layers and get the DSVs
yearlyval <- function(subind, dayind) {
  # Raster
  raszone = rast(sublists[[subind]])
  
  # Sequence for day of year
  dayseq = seq(dayind, nlyr(raszone), by = 365)
  
  # Values
  dayvalue = values(raszone[[dayseq]], na.rm=T, dataframe =T)
  colnames(dayvalue) = paste0("y",2009:2020)
  dayvalue$day = dayind
  
  return(dayvalue)
}

# Apply the function through the classes of lists
yearcol_list <- lapply(1:19, function(subind) {
  yearcol_sub <- lapply(1:365, function(dayind) yearlyval(subind, dayind))
  yearwrite = do.call(rbind, yearcol_sub)
  write.csv(yearwrite, paste0("Z:/Late_blight/lbdsveconu/appscores/csv/sub_",subind,".csv"),
            row.names = F)
})

# The CSVs written are huge... up to 6GB. Manipulating them on local desktop is difficult,
# even with data.table. Out of memory errors are returned. 



