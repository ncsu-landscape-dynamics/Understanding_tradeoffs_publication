library(terra)
library(tidyverse)

# List of daily DSVs
ldr1 <- list.files()
# List of cumulative DSVs
lcr1 <- list.files()

# Read in DSV 
dsvr = rast(ldr1[1])

# Read in cumulative DSV
cudsvr = rast(lcr1[1])

# Read in infections table
infectiontbl = read.csv()
# Modify table
infectiontbl = infectiontbl %>% filter(!State %in% c("AZ", "HI", "VI", "MP", "GU", "AS", "PR", "IA", "KS","MN","ND","ID","NE","OR","SD","AK","LA","MT","CA","CO","NB","UT","WA","AR","MO","NM","OK","WY","NV","MS","DC","TX"))
infectiontbl = infectiontbl %>% group_by(year, State, week, County) %>% mutate(n_infx = n()) %>% slice_head()
infectiontbl = infectiontbl %>% arrange(year, yday)

# Read in counties
counties = vect()

dsvinfext = function( indx, ...){
  # Select infection day
  infectionday = as.numeric(infectiontbl[indx, "yday"])
  infect_date = as.numeric(infectiontbl[indx, "Samp_date"])
  infection_year = as.numeric(infectiontbl[indx, "year"])
  
  # If sday is < 14 Jan., need to add in previous year's raster. This won't
  # work for before 2009.
  # dr1 = c(dr1, rast(year - 1))
  if(infectionday < 14){
    dsvr = c(dsvr, rast(ldr1[grepl((infection_year-1), ldr1)]))
    cudsvr = c(cudsvr, rast(lcr1[grepl((infection_year-1), lcr1)]))
  }
  cuml =   cudsvr[[1:infectionday]]#cd09[[1:sday]]
  cuml =   cuml[[nlyr(cuml)]]
  biweek = dsvr[[(infectionday-13):infectionday]]#cd09[[(sday-14):sday]]
  biweek = app(biweek, cumsum)
  biweek = biweek[[nlyr(biweek)]]
  
  # Identify location to match day
  infct_loc = infectiontbl[indx, "st_ct"]
  
  # Get the county to match infection location
  co_infct_loc = counties[counties$stco == infct_loc,]
  
  # Crop DSV 
  dsv_loc = crop(dsvr, co_infct_loc, mask=T)
  cudsv_loc = crop(cudsvr, co_infct_loc, mask=T)
  
  # Summary of cropped DSV
  dsvmean = mean(dsv_loc, na.rm=T)
  dsvmedn = median(dsv_loc, na.rm=T)
  
  cudsvmean = mean(cudsv_loc, na.rm=T)
  cudsvmedn = median(cudsv_loc, na.rm=T)
  
  # Combine summaries
  dsvdf = data.frame(infxn_date = infect_date, infct_loc = infct_loc, dsvlocmean=dsvmean, 
                     dsvlocmedian = dsvmedn, cmldsvmean = cudsvlocmean, cmldsvlocmedian = cudsvmedn)
  
  return(dsvdf)
}

