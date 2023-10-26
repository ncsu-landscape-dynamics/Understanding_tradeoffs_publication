library(terra)
library(tidyverse)
library(lubridate)
library(tigris)
library(sf)

infx <-read.csv( "Q:/My Drive/P_infestans_data/USAB_NPDN_combined_USETHIS.csv")
infx <- infx[,-1]
infx <- infx %>% filter(!State %in% c("KS","MN","ND","ID","NE","OR","SD","AK","LA","MT","CA","CO","NB","UT","WA","AR","MO","NM","TX"))
infxco <- infx %>% group_by(year, State, week, County) %>% mutate(n_infx = n()) %>% slice_head()
stls <- unique(infxco$State)

#co1 <- counties(stls)
# This takes a long time
#co1 <- st_transform(co1, crs="epsg:3857")
# Saved from earlier
co1 <- vect("Q:/My Drive/temp/counties3857.shp")
st1 <- states()
st2 <- as.data.frame(st1[,c(3,6)])
# Use letter codes instead of FIPS, easier to understand
co1$stcd <- st2[match(co1$STATEFP, st2$STATEFP), "STUSPS"]
# Make some names more easiiy interchangable with other data set
co1$NAME <- gsub("St\\.", "St", co1$NAME)
co1$NAME <- gsub("Fond du", "Fond Du", co1$NAME)
co1$NAME <- gsub("LaGrange", "Lagrange", co1$NAME)
co1$NAME <- gsub("De [K|k]alb", "Dekalb", co1$NAME)
co1$NAME <- gsub("'s","s", co1$NAME)

# State templates 
statemp <- list.files("Q:/My Drive/popblightdata/templaterastershapes1k/", pattern="tif$", full.names = T)
crln <- str_sub(statemp, 50,51)

# Hosts
hrls1 <- list.files("Q:/My Drive/popblightdata/newareas/binaryfrompottom/", pattern="tif$", full.names = T)
hrln <- str_sub(hrls1, 48,51)

# Select the state raster based on row in infections data set
stsel_f <- function(x, y){
  stsel = as.character(x[y, "State"])
  styse = x[y, "year"]
  strsel = statemp[which(crln == stsel)]
  str1 = rast(strsel)
  str1[str1 < 0] = NA
  str1[str1 > 0] = 0
  names(str1) = "infection"
  return(str1)
}

# Build a raster with each week
bd_rast <- function(x) {
  str2 = c(x, x)
  for(i in 3:52) str2 = c(str2, x)
  return(str2)
}

# Subset to county and set values to 0, with cells outside boundary NA
ctysub <- function(x, y, selst) {
  #wk1 = x[y, "week"]
  cosel = as.character(x[y, "County"])
  print(cosel)
  stasel = as.character(x[y, "State"])
  print(stasel)
  co2 = co1[co1$NAME == cosel & co1$stcd == stasel,]
  print(co2)
  selco = crop(selst, co2, mask=T)
  selco[selco < 0] = NA
  selco[selco > 0] = 0
  return(selco)
}

# Add values to the county. Whole county is assigned value.
# Write out the county raster. It will be used later.
addvl <- function(x, y, ctysub, stsel_f){#, bd_rast) {
  selco = ctysub
  selco[selco == 0] = x[y, "n_infx"]
  tesl = mosaic(stsel_f, selco, fun=sum)
  st1 = x[y, "State"]
  co1 = x[y, "County"]
  wk1 = x[y, "week"]
  if(str_length(wk1) == 1){
    wk1 = paste0(0,wk1)
  } else {
    wk1
  }
  writeRaster(tesl, paste0("Q:/My Drive/popblightdata/newareas/infections/", st1,"_wk",wk1,"_",co1,".tif"))
}

# Subset the infection collection to a year
infyr <- infxco[infxco$year == 2010 ,]

# Loop through the year's infections and create files for each row (county & year)
for(i in 1:nrow(infyr)){
  y = infyr
  # Get state
  selst = stsel_f(y, i)
  # Build year's worth
  yrrst = bd_rast(selst)
  # County and week
  countysub = ctysub(y, i, selst)
  # Add values
  add_val = addvl(y, i, countysub, selst)
  print(add_val)
}

# Read in list of names of raster files for a state & year that were saved
# in earlier step
infrls <- list.files("Q:/My Drive/popblightdata/newareas/infections/", pattern = "tif$", 
                     full.names = T)

# Get the state and week of each file from the raster name
stnl <- str_sub(infrls, 47, 48)
stnl <- unique(stnl)
wkl1 <- str_sub(infrls, 52, 53)

# Combine the counties, via addition, for a given week for a given state
combr <- function(x,i){
  for(j in x$lengths[x$values == uwk1[i]]){
    print(paste0("j:",j))
    smrstl = infrst[names(infrst) == uwk1[i]]
    print(smrstl)
    #temr = 
    smrs = rast(smrstl[1])+rast(smrstl[2])
    if(length(smrstl) > 2){
      for(m in 3:length(smrstl)){
        print(paste0("m:",m))
        smrs = smrs+rast(smrstl[m])
        return(smrs)
      }
    }
  } 
  return(smrs)
} 

assign_val <- function(uwk1){
  for(k in 1:length(uwk1)){
    print(paste0("k:",k))
    if(wrle1$lengths[wrle1$values==uwk1[k]] == 1){
      yrrst[[as.numeric(wrle1$values[which(wrle1$values==uwk1[k])])]] = rast(infrst[k])
    } else {
      yrrst[[as.numeric(wrle1$values[which(wrle1$values==uwk1[k])])]] = combr(wrle1, k)
    }
  }
  return(yrrst)
} 

# This doesn't work as a loop. Have to manually change the indexe for stnl
# Probably need to make this a function and loop that way.
for(i in stnl[18]){
  inrsl = infrls[grepl(i, infrls)]
  infrst = lapply(infrls[grepl(i, infrls)], rast
  # State code to add to file name
  stlet = unique(str_sub(inrsl, 47, 48))
  strst = statemp[grepl(i, statemp)]
  strst = rast(strst)
  strst[strst <= 0] = NA
  strst[strst > 0] = 0
  names(strst) <- "no_infection"
  wkl1 = str_sub(inrsl, 52, 53)
  uwk1 = unique(wkl1)
  print(paste0("uwk1:",uwk1))
  names(infrst) <- wkl1
  wrle1 = rle(wkl1)
  wrle2 = wrle1$values[which(wrle1$lengths > 1)]
  yrrst = bd_rast(strst)
  temout = assign_val(uwk1)
  #return(temout)
  writeRaster(temout, paste0( outpath, stlet, "_2010_infections.tif"))
}

