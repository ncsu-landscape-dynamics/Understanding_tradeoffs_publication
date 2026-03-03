library(terra)
library(tidyverse)
library(lubridate)
library(tigris)
library(sf)

sodsv <- vect("Z:/Late_blight/helpful_shapes_rasters/SOD_all_3857.gpkg")

#table(sodsv$YEAR)
#2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 
#55   81   57   28   60  143  162  122  102  120  337  207  188  182   88  222 
#2017 2018 2019 2020 
#256  260  198   81 

#1   2     3    4    5     6   7    8    9   10   11 
#7  48    49   25   14   125  36   35    7    3    2
1    2    3    4    5    6    7    8    9   10   11 
0.02 0.14 0.14 0.07 0.04 0.36 0.10 0.10 0.02 0.01 0.01

odds1 <- c(0.02, 0.14, 0.14, 0.07, 0.04, 0.36, 0.1, 0.1, 0.02, 0.01, 0.01)

alldt <- seq(from = as.Date("2001-02-01"), to=as.Date("2020-12-31"), by=1)
# Create a new data frame with all dates
df_all_dates <- data.frame(Date = alldt, Max_14 = NA,
                           Min_14 = NA,
                           Mean_14 = NA,
                           Median_14 = NA,
                           SD_14 = NA,
                           Max_30 = NA,
                           Min_30 = NA,
                           Mean_30 = NA,
                           Median_30 = NA,
                           SD_30 = NA
)

l3 <- list.files("Z:/Late_blight/SODdsv/", full.names = T, recursive = T, pattern = 
                   "dsvdai.*tif$")

extracted_soddsv <- data.frame(ID = as.character(),
                            Date = as.character(),
                            #st_ct = as.character(),
                            Max_14 = as.numeric(),
                            Min_14 = as.numeric(),
                            Mean_14 = as.numeric(),
                            Median_14 = as.numeric(),
                            SD_14 = as.numeric(),
                            Max_30 = as.numeric(),
                            Min_30 = as.numeric(),
                            Mean_30 = as.numeric(),
                            Median_30 = as.numeric(),
                            SD_30 = as.numeric(),
                            #FIPS = as.character(),
                            stringsAsFactors = FALSE)


for (i in 1:nrow(sod2)){  #nrow(spinf)) {
  # Get the date from the current row
  print(paste0("Index: ",i))
  yearsel <- as_date(paste0(sodsv$YEAR[i],"-01-01"))
  
  # Extract the year and day from the date
  syear <- year(yearsel)
  
  # Sample for a month in the year. Weights are based on 2021 and 2022 data
  sampl_mo = 1:12
  sampl_wt = c(0.02, 0.14, 0.14, 0.07, 0.04, 0.34, 0.1, 0.1, 0.02, 0.01, 0.01, 0.01)
  pt_sampl_mo = sample(sampl_mo, 1, replace = T, sampl_wt)
  print(pt_sampl_mo)
  
  # Sample for a day to make the date complete.
  num_days <- days_in_month(as.Date(paste(syear, pt_sampl_mo, "01", sep="-")))
  print(num_days)
  # Randomly sample a day
  sampl_day <- sample(1:num_days, 1)
  
  # Create the date
  date <- as_date(paste(syear, pt_sampl_mo, sampl_day, sep="-"))
  print(date)
  # date as doy
  sday <- yday(date)
  
  # Find the raster file that matches the year and day
  dr1 <- rast(l3[grepl(syear, l3)])#paste0(syear,"masked"), l3)])
  
  # Get the dates for the preceding 14 days
  dates_14 <- seq(date - 14, date, by = "day")
  
  # Get the dates for the preceding 30 days
  dates_30 <- seq(date - 30, date, by = "day")
  
  # Create list of indices for the dates_14
  indx14 <- rev(seq(sday, length.out = 14, by = -1))
  
  # Create list of indices for the dates_30
  indx30 <- rev(seq(sday, length.out = 30, by = -1))
  
  if(sday < 30){
    # Prepend the earlier year to the selected year
    dr2 = c(rast(l3[grepl((syear-1), l3)]), dr1)
    
    # Figure out the number of days for prepended year
    dayext = nlyr(rast(l3[grepl( (syear-1), l3)]))
    
    # No. of days for selected year
    newday1 = seq(as_date(paste0(syear, "-01-01")), length.out = 365, by = 1)
    # No. of days for prepended year
    newday2 = rev(seq(as_date(paste0(syear, "-01-01")) - days(1), length.out = dayext, by = -1))
    # Add dates to the raster
    time(dr2) = c(newday2, newday1)
    dr1 = dr2
    # Adjust the fortnight for prepended year
    sdaynew = sday + dayext
    dates_14 = rev(seq(sdaynew, length.out=14, by = -1))
    #cr1 = c(rast(l1[grepl((syr1-1), l1)]), cr1)
    # Adjust the month for prepended year
    dates_30 = rev(seq(sdaynew, length.out=30, by = -1))
    indx14 = rev(seq(sdaynew, length.out = 14, by = -1))
    indx30 = rev(seq(sdaynew, length.out = 30, by = -1))
  } else {
    dr1 = dr1
    time(dr1) = seq(as_date(paste0(syear, "-01-01")), length.out = 365, by = 1)
    dates_14 = dates_14
    dates_30 = dates_30
    #cr1 = cr1
  } 
  

    # Test if the raster has values for the dates_14 in the location for extract with spinf. If not, use 0 for extract instead of NA.
  dsv_14 <- terra::extract(dr1[[indx14]], sodsv[i,])
  print(dsv_14)
  dsv_14 <- as.numeric(dsv_14)
  print(dsv_14)
  
  #if (any(is.na(terra::extract(dr1[[indx30]], spinf[i, ])))) {
  #  dsv_30 <- rep(0, length(indx30))
  #} else {
  dsv_30 <- terra::extract(dr1[[indx30]], sodsv[i,])
  #  #dsv_30[is.na(dsv_30)] <- 0
  #}  
  
  dsv_30 <- as.numeric(dsv_30)
  print(dsv_30)
  
  # Calculate the max, min, mean, median, and sd for each set
  sum_14 <- sum(dsv_14[2:15], na.rm = TRUE)
  mean_14 <- mean(dsv_14[2:15], na.rm = TRUE)
  median_14 <- median(dsv_14[2:15], na.rm = TRUE)
  sd_14 <- sd(dsv_14[2:15], na.rm = TRUE)
  
  sum_30 <- sum(dsv_30[2:31], na.rm = TRUE)
  mean_30 <- mean(dsv_30[2:31], na.rm = TRUE)
  median_30 <- median(dsv_30[2:31], na.rm = TRUE)
  sd_30 <- sd(dsv_30[2:31], na.rm = TRUE)
  
  # Add the extracted DSV values to the data frame
  extracted_dsv <- rbind(extracted_dsv, data.frame(ID = sodsv$SAMPLE__[i], #spinf$RptId[i],
                                                   Date = date,
                                                   Sum_14 = round(sum_14, 1),
                                                   #Min_14 = round(min_14, 1),
                                                   Mean_14 = round(mean_14, 1),
                                                   Median_14 = round(median_14, 1),
                                                   SD_14 = round(sd_14, 1),
                                                   Sum_30 = round(sum_30, 1),
                                                   #Min_30 = round(min_30, 1),
                                                   Mean_30 = round(mean_30, 1),
                                                   Median_30 = round(median_30, 1),
                                                   SD_30 = round(sd_30, 1),
                                                   stringsAsFactors = FALSE))
  
}

extracted_dsv$Year <- year(extracted_dsv$Date)
extracted_dsv$yday <- yday(extracted_dsv$Date)

write.csv(extracted_dsv, "Z:/Late_blight/SODdsv/sod_extr_dsvsum.csv", row.names = F)

### I've also called this dsv2roc.R
