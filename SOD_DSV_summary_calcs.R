library(terra)
library(tidyverse)
library(lubridate)
library(tigris)
library(sf)


# Load the shapefile containing information about blight in North America in 2023
spinf <- vect("Z:/Late_blight/helpful_shapes_rasters/blight_NA_2023.shp")
spinf <- spinf[spinf$Cnry == "United States",]
spinf$CTY <- gsub(" County", "", spinf$CTY) # remove the space in the county name

infx <-read.csv( "Q:/My Drive/P_infestans_data/USAB_NPDN_combined_USETHIS.csv")
infx$yday <- yday(infx$ObsDt)
#infx <- read.csv("Z:/Late_blight/relatedtables/infxco2.csv")
infxco <- infx %>% filter(!ST %in% c("AZ", "HI", "VI", "MP", "GU", "AS", "PR", 
                                     "IA", "KS","MN","ND","ID","NE","OR","SD",
                                     "AK","LA","MT","CA","CO","NB","UT","WA","AR",
                                     "MO","NM","OK","WY","NV","MS","DC","TX"))
infxco <- infxco %>% filter(!st_ct %in% c("IN-West Lafayette", "NJ-Agmart", "NY-Cutches", 
                                          "NY-Highland", "PA-PennState"))
#infxco <- infxco %>% group_by(year, State, week, County) %>% mutate(n_infx = n()) %>% slice_head()
infxco <- infxco %>% arrange(Year, yday)
infxco$week <- week(infxco$ObsDt)
stls <- unique(infxco$ST)
spinf <- spinf[spinf$ST %in% stls,]
spinf$ObsDt <- as_date(as_datetime(as.numeric(spinf$newdate)))
spinf$yday <- yday(spinf$ObsDt)
spinf$CTY <- gsub("St\\.", "St", spinf$CTY)
spinf$CTY <- gsub("Fond du", "Fond Du", spinf$CTY)
spinf$CTY <- gsub("LaGrange", "Lagrange", spinf$CTY)
spinf$CTY <- gsub("De [K|k]alb", "Dekalb", spinf$CTY)
spinf$CTY <- gsub("'s","s", spinf$CTY)
spinf$st_ct <- paste0(spinf$ST,"-",spinf$CTY)
spinf <- spinf[!spinf$st_ct %in% c("IN-West Lafayette", "NJ-Agmart", "NY-Cutches", 
                                          "NY-Highland", "PA-PennState"),]
#co1 <- counties(stls)
# This takes a long time
#co1 <- st_transform(co1, crs="epsg:3857")
# Saved from earlier
co1 <- vect("D:/data/ncsu/shapesothergis/counties3857.shp")
st1 <- states()
st1 <- st1[st1$STUSPS %in% stls,]
st2 <- as.data.frame(st1[,c(3,6)])
st2 <- st2[,1:2]

# Use letter codes instead of FIPS, easier to understand
co1 <- co1[co1$STATEFP %in% st2$STATEFP,]

# Make some names more easiiy interchangable with other data set
co1$NAME <- gsub("St\\.", "St", co1$NAME)
co1$NAME <- gsub("Fond du", "Fond Du", co1$NAME)
co1$NAME <- gsub("LaGrange", "Lagrange", co1$NAME)
co1$NAME <- gsub("De [K|k]alb", "Dekalb", co1$NAME)
co1$NAME <- gsub("'s","s", co1$NAME)

co1 <- co1[co1$GEOID != 24510,] # the REAL Baltimore is 24005

co3 <- as.data.frame(co1)
co3 <- co3[,c(1,2,4, 18, 19)]
co3 <- co3[co3$GEOID != 24510,]

infxco$FIPS <- as.character(co3[match(infxco$st_ct, co3$stco), "GEOID"])
spinf$FIPS <- as.character(co3[match(spinf$st_ct, co3$stco), "GEOID"])

# Or...
co2 <- vect("D:/data/ncsu/googledrop/newareas/infections/infections/county_infections_2009.shp")
co2$FIPS <- paste0(co2$STATEFP,co2$COUNTYFP)
# Take out Baltimore City (it's independent)
co2 <- co2[co2$FIPS != 24510,] # the REAL Baltimore is 24005

# A list of daily dsvs.
l2 <- list.files("Z:/Late_blight/lbdsveconu/", pattern="dsv[0-9]{4}masked.tif", recursive=T, full.names=T)
#l2 <- list.files("Z:/Late_blight/lbdsveconu/", pattern="maskcml.*tif", recursive=T, full.names=T)

l4 <- rast("D:/data/ncsu/googledrop/newareas/DSV/2009/dsv2009masked.tif")

# A list of daily cumulative dsvs.
l3 <- list.files("Q:/My Drive/popblightdata/newareas/DSV/2009/cumulthreshold/rasters/original_singles/",
                 full.names = T, pattern = "tif$")
l3 <- mixedsort(l3)
l3 <- rast("D:/data/ncsu/googledrop/popblightdata/cudsv2009.tif")
l5 <- rast("Q:/My Drive/temp/dsvs/2009/cumulative/cumuthresh_2009.tif")

rt1 = list()
rtw1 = list()

l1 <- list.files("D:/data/ncsu/lbdsveconu/", pattern = "cumld", full.names = T, 
                 recursive = T)
l1 <- c(l1[1], l1)

# Read in a list of year files with daily DSV values
l2 <- list.files("Z:/Late_blight/lbdsveconu/", pattern = "dsv_[0-9].*tif", full.names = T, 
                 recursive=T)
l2 <- l2[c(1:2,4:14)]
#l2 <- list.files("D:/data/ncsu/lbdsveconu/", pattern = "dsv_[0-9].*tif$", full.names = T, 

#                 recursive = T)
#l2 <- l2[2:length(l2)]




# Create a spatvect for all the points in infxco. Use coordinates from USAB and
# county center coordinates for NPDS.
spinfm <- spinf[, c(1,20,21)]
infxco2 <- merge(infxco, spinfm, by.x="Rec_ID", by.y="RptId", all.x=T)
delFIPS <- infxco2[which(is.na(infxco2$x)), "FIPS"]
delco <- co1[co1$GEOID %in% delFIPS, c(4, 16, 17)]
delco$INTPTLAT <- as.numeric(delco$INTPTLAT)
delco$INTPTLON <- as.numeric(delco$INTPTLON)
infxco3 <- merge(infxco2, delco, by.x="FIPS", by.y="GEOID", all.x=T)
infxco3$x <- if_else(is.na(infxco3$x), infxco3$INTPTLON, infxco3$x)
infxco3$y <- if_else(is.na(infxco3$y), infxco3$INTPTLAT, infxco3$y)
infxco3 <- infxco3 %>% select(2:14,1,15,16)
infpts <- vect(infxco3, geom=c("x","y"), crs="epsg:4326")
infpts <- infpts[order(infpts$Year, infpts$yday),] # Should've used ObsDt
#writeVector(infpts, "Z:/Late_blight/helpful_shapes_rasters/pointsUSABcrdsNPDNctycenters.gpkg")

spinf <- project(spinf, crs(rast(l2[1])))


dr1 <- NULL

infpts <- project(infpts, crs(rast(l2[1])))

extracted_dsv <- data.frame(ID = as.character(),
                            Date = as.character(),
                            st_ct = as.character(),
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
                            FIPS = as.character(),
                            stringsAsFactors = FALSE)


# Run the extracted_dsv chunk before this each time. 
#Loop through each row of spinf
for (i in 1:nrow(infpts)){  #nrow(spinf)) {
  # Get the date from the current row
  date <- as_date(infpts$ObsDt[i])#spinf$ObsDt[i]
  print(date)
  # Extract the year and day from the date
  syear <- year(date)
  sday <- yday(date)
  
  # Find the raster file that matches the year and day
  dr1 <- rast(l2[grepl(syear, l2)])#paste0(syear,"masked"), l2)])
  
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
    dr2 = c(rast(l2[grepl((syear-1), l2)]), dr1)
  
    # Figure out the number of days for prepended year
    dayext = nlyr(rast(l2[grepl( (syear-1), l2)]))
    
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
  dsv_14 <- terra::extract(dr1[[indx14]], spinf[i, ])
  
  dsv_14 <- as.numeric(dsv_14)
  print(dsv_14)
  
  #if (any(is.na(terra::extract(dr1[[indx30]], spinf[i, ])))) {
  #  dsv_30 <- rep(0, length(indx30))
  #} else {
  dsv_30 <- terra::extract(dr1[[indx30]], spinf[i, ])
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
  extracted_dsv <- rbind(extracted_dsv, data.frame(ID = infpts$Rec_ID[i], #spinf$RptId[i],
                                                   Date = date,
                                                   st_ct = infpts$st_ct[i],  #spinf$st_ct[i],
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
                                                   FIPS = infpts$FIPS[i], #spinf$FIPS[i],
                                                   stringsAsFactors = FALSE))
  
}

extracted_dsv$Year <- year(extracted_dsv$Date)
extracted_dsv$yday <- yday(extracted_dsv$Date)

write.csv(extracted_dsv, "Z:/Late_blight/lbdsveconu/lb20extracted_nomask_dsv.csv", row.names = F)

---:SOD


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
