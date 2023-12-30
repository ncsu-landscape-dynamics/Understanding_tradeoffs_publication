# This is for counties that were actually infected. There is a similar file
# DSV_summaryval_calculation_random_uninfected_counties.R
library(terra)
library(tidyverse)

l1 <- list.files("...", pattern="dsv.*masked.tif", recursive=T, full.names=T)
l2 <- list.files("...", pattern="maskcml.*tif", recursive=T, full.names=T)

infx <- read.csv( ".../USAB_NPDN_combined_USETHIS.csv")
infxco <- infx %>% filter(!State %in% c("AZ", "HI", "VI", "MP", "GU", "AS", "PR", "IA", "KS","MN","ND","ID","NE","OR","SD","AK","LA","MT","CA","CO","NB","UT","WA","AR","MO","NM","OK","WY","NV","MS","DC","TX"))
infx_df <- infxco %>% arrange(year, yday)

co1 <- vect(".../counties3857.shp")

dsvext <- function(infdf, indx, cr1, dr1, ...){ #cr1 = cumul, dr1 = daily
  print(indx)
  # Get the days/durations relevant to each period
  sday = as.numeric(infdf$yday)
  print(sday)
  syr1 = as.numeric(infdf$year)
  # If sday is < 30 Jan., need to add in previous year's raster. This won't
  # work for before 2009.
  if(sday < 30){
    dr1 = c(rast(l2[grepl((syr1-1), l2)]), dr1)
    cr1 = c(rast(l1[grepl((syr1-1), l1)]), cr1)
    sday = sday + 365
  } else {
    dr1 = dr1
    cr1 = cr1
  }
  # The annual cumulative values
  cuml    = cr1[[1:sday]]
  cuml    = cuml[[nlyr(cuml)]]
  # The cumulative fortnight values
  fortnt  = dr1[[(sday-13):sday]]
  fortnt  = app(fortnt, cumsum)
  fortnt  = fortnt[[nlyr(fortnt)]]
  # The cumulative monthly values
  monthly = dr1[[ (sday-29):sday]]
  monthly = app(monthly, cumsum)
  monthly = monthly[[nlyr(monthly)]]

  # Get the state and county for selected row
  selecstcty = as.character(infdf[indx, "st_ct"])
  print(paste0("selectct ", selecstcty))
  # County matching the selected location
  stecty = co1[co1$stco == selecstcty,]
  print(paste0("sel cty ", as.character(as.data.frame(stecty)[1,6])))
  # Crop DSVs to selected county
  dsvcumlcrp = crop(cuml, stecty, mask=T)
  dsvftntcrp = crop(fortnt, stecty, mask=T)
  dsvmnthcrp = crop(monthly, stecty, mask=T)
  
  # Calcuate the summary stats for each accumulated period
  cmlmin = min(freq(dsvcumlcrp)$value)
  cmlmax = max(freq(dsvcumlcrp)$value)
  cmlavg = weighted.mean(freq(dsvcumlcrp)$value, freq(dsvcumlcrp)$count)
  cmlsd = sd(freq(dsvcumlcrp)$value)
  cmlmed = median(freq(dsvcumlcrp)$value)

  fntmin =  min(freq(dsvftntcrp)$value)
  fntmax =  max(freq(dsvftntcrp)$value)
  fntavg = weighted.mean(freq(dsvftntcrp)$value, freq(dsvftntcrp)$count)
  fntsd =  sd(freq(dsvftntcrp)$value)
  fntmed = median(freq(dsvftntcrp)$value)

  mnthmin =  min(freq(dsvmnthcrp)$value)
  mnthavg = weighted.mean(freq(dsvmnthcrp)$value, freq(dsvmnthcrp)$count)
  mnthsd =  sd(freq(dsvmnthcrp)$value)
  mnthmed = median(freq(dsvmnthcrp)$value)
  
  # combine
  ddf1 = data.frame(loc=as.character(infdf[indx,"st_ct"]), date=as.character(infdf[indx,"Samp_date"]),
                    cmlmin=cmlmin, cmlmax=cmlmax, cmlavg=cmlavg, cmlsd=cmlsd,
                    cmlmedian = cmlmed, fortntmin=fntmin, fortntmax=fntmax,
                    fortntavg=fntavg, fortntsd=fntsd, fortntmedian = fntmed,
                    mnthmin=mnthmin, mnthmax=mnthmax, mnthavg=mnthavg,
                    mnthsd=mnthsd, mthmedian = mnthmed)
  return(ddf1)

}

infx9 <- infx_df %>% filter(year == 2009) %>% arrange(yday)
infx10 <- infx_df %>% filter(year == 2010) %>% arrange(yday)
infx11 <- infx_df %>% filter(year == 2011) %>% arrange(yday)
infx12 <- infx_df %>% filter(year == 2012) %>% arrange(yday)
infx13 <- infx_df %>% filter(year == 2013) %>% arrange(yday)
infx14 <- infx_df %>% filter(year == 2014) %>% arrange(yday)
infx15 <- infx_df %>% filter(year == 2015) %>% arrange(yday)
infx16 <- infx_df %>% filter(year == 2016) %>% arrange(yday)
infx17 <- infx_df %>% filter(year == 2017) %>% arrange(yday)
infx18 <- infx_df %>% filter(year == 2018) %>% arrange(yday)
infx19 <- infx_df %>% filter(year == 2019) %>% arrange(yday)
infx20 <- infx_df %>% filter(year == 2020) %>% arrange(yday)

# Note that running this takes a couple of days on HPC for each year of first 6 years.  
# The later years have fewer reports and can run in less than a day for each year.

cr1 <- rast(l1[1])
dr1 <- rast(l2[1])
del09 <- lapply(12:nrow(infx9), function(x) dsvext(infx9, x, cr1, dr1))
del09 <- do.call(rbind, del09)
write.csv(del09, ".../dsvstats2009.csv", row.names = F)

cr1 <- rast(l1[2])
dr1 <- rast(l2[2])
del10 <- lapply(1:nrow(infx10), function(x) dsvext(infx10, x, cr1, dr1))
del10 <- do.call(rbind, del10)
write.csv(del10, ".../dsvstats2010.csv", row.names = F)

#Sys.time()
cr1 <- rast(l1[3])
dr1 <- rast(l2[3])
del11 <- lapply(1:nrow(infx11), function(x) dsvext(infx11, x, cr1, dr1))
del11 <- do.call(rbind, del11)
write.csv(del11, ".../dsvstats2011.csv", row.names = F)

cr1 <- rast(l1[4])
dr1 <- rast(l2[4])
del12 <- lapply(1:nrow(infx12), function(x) dsvext(infx12, x, cr1, dr1))
del12 <- do.call(rbind, del12)
write.csv(del12, ".../dsvstats2012.csv", row.names = F)

cr1 <- rast(l1[5])
dr1 <- rast(l2[5])
del13 <- lapply(1:nrow(infx13), function(x) dsvext(infx13, x, cr1, dr1))
del13 <- do.call(rbind, del13)
write.csv(del13, ".../dsvstats2013.csv", row.names = F)

cr1 <- rast(l1[6])
dr1 <- rast(l2[6])
del14 <- lapply(1:infx14, function(x) dsvext(infx14, x, cr1, dr1))
del14 <- do.call(rbind, del14)
write.csv(del14, ".../dsvstats2014.csv", row.names = F)

cr1 <- rast(l1[7])
dr1 <- rast(l2[7])
del15 <- lapply(1:nrow(infx15), function(x) dsvext(infx15, x, cr1, dr1))
del15 <- do.call(rbind, del15)
write.csv(del15, ".../dsvstats2015.csv", row.names = F)

cr1 <- rast(l1[8])
dr1 <- rast(l2[8])
del16 <- lapply(1:nrow(infx16), function(x) dsvext(infx16, x, cr1, dr1))
del16 <- do.call(rbind, del16)
write.csv(del16, ".../dsvstats2016.csv", row.names = F)

cr1 <- rast(l1[9])
dr1 <- rast(l2[9])
del17 <- lapply(1:nrow(infx17), function(x) dsvext(infx17, x, cr1, dr1))
del17 <- do.call(rbind, del17)
write.csv(del17, ".../dsvstats2017.csv", row.names = F)

cr1 <- rast(l1[10])
dr1 <- rast(l2[10])
del18 <- lapply(1:nrow(infx18), function(x) dsvext(infx18, x, cr1, dr1))
del18 <- do.call(rbind, del18)
write.csv(del18, ".../dsvstats2018.csv", row.names = F)

cr1 <- rast(l1[11])
dr1 <- rast(l2[11])
del19 <- lapply(1:nrow(infx19), function(x) dsvext(infx19, x, cr1, dr1))
del19 <- do.call(rbind, del19)
write.csv(del19, ".../dsvstats2019.csv", row.names = F)

cr1 <- rast(l1[12])
dr1 <- rast(l2[12])
del20 <- lapply(1:nrow(infx20), function(x) dsvext(infx20, x, cr1, dr1))
del20 <- do.call(rbind, del20)
write.csv(del20, ".../dsvstats2020.csv", row.names = F)

