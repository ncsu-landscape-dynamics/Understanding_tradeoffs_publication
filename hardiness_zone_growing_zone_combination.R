# Seasons by area are based on the dates and areas found here:
'https://www.tomatofest.com/Tomato_Growing_Zone_Maps_s/164.htm '

library(terra)
library(tidyverse)
library(lubridate)

l1 <- list.files("Z:/Late_blight/lbdsveconu/hzone", full.names = T")
hz <- gtools::mixedsort(unique(str_split_i(str_split_i(l1, "/", 5), "_", 2)))

# Use 4, (4b&5a), (5b&6a), 
4 - darkblue    - late May:mid-Sept
5 - blue        - late May:mid-Sept

6 - lightblue   - early May:early Oct

7 - brown       - mid Apr: Oct
8 - orange      - 
9 - yellow
10 - green
11 - lightgreen

"Jan: Dec"
"early Feb: late Dec"
"mid Mar: early Nov"
"mid Mar: early Dec"
"late Mar: early Nov"
"mid Apr: Oct"
"late Apr: late Oct"
"early May: early Oct"
"May: early Oct"      
"May: mid Oct"
"May: late Oct"
"late May: Sept"
"late May:mid Sept"
"late May: late Sept"
"Jun: mid Sept"

grwson$NOTE <- factor(grwson$NOTE, levels = c("Jan: Dec",
                                              "early Feb: late Dec",
                                              "mid Mar: early Nov",
                                              "mid Mar: early Dec",
                                              "late Mar: early Nov",
                                              "mid Apr: Oct",
                                              "late Apr: late Oct",
                                              "early May: early Oct",
                                              "May: early Oct",      
                                              "May: mid Oct",
                                              "May: late Oct",
                                              "late May: Sept",
                                              "late May:mid Sept",
                                              "late May: late Sept",
                                              "Jun: mid Sept"
))

# Several of these dates are modified. 
grwson$season_st <- case_when(str_detect(grwson$NOTE, "Jan") == T ~ yday(as.Date("2009/01/01")),
                              str_detect(grwson$NOTE, "early Feb") == T ~ yday(as.Date("2009/02/03")),
                              str_detect(grwson$NOTE, "mid Mar") == T ~ yday(as.Date("2009/03/15")),                              
                              str_detect(grwson$NOTE, "late Mar") == T ~ yday(as.Date("2009/03/24")),
                              str_detect(grwson$NOTE, "mid Apr") == T ~ yday(as.Date("2009/04/15")),
                              str_detect(grwson$NOTE, "late Apr") == T ~ yday(as.Date("2009/04/24")),
                              str_detect(grwson$NOTE, "early May") == T ~ yday(as.Date("2009/02/03")),
                              str_detect(grwson$NOTE, "May") == T ~ yday(as.Date("2009/05/01")),
                              str_detect(grwson$NOTE, "late May") == T ~ yday(as.Date("2009/05/24")),
                              str_detect(grwson$NOTE, "Jun") == T ~ yday(as.Date("2009/06/01"))
                              )

grwson$season_end <- case_when(str_detect(grwson$NOTE, "Dec") == T ~ yday(as.Date("2009/12/31")),
                              str_detect(grwson$NOTE, "late Dec") == T ~ yday(as.Date("2009/12/24")),
                              str_detect(grwson$NOTE, "early Nov") == T ~ yday(as.Date("2009/11/05")),                              
                              str_detect(grwson$NOTE, "Oct") == T ~ yday(as.Date("2009/10/01")),
                              str_detect(grwson$NOTE, "late Oct") == T ~ yday(as.Date("2009/10/24")),
                              str_detect(grwson$NOTE, "early Oct") == T ~ yday(as.Date("2009/10/05")),
                              str_detect(grwson$NOTE, "mid Oct") == T ~ yday(as.Date("2009/10/14")),
                              str_detect(grwson$NOTE, "Sept") == T ~ yday(as.Date("2009/09/01")),
                              str_detect(grwson$NOTE, "mid Sept") == T ~ yday(as.Date("2009/09/14")),
                              str_detect(grwson$NOTE, "late Sept") == T ~ yday(as.Date("2009/09/24"))
)

# The file Z:/Late_blight/relatedtables/seasondatesdirs.csv has the dates used to delineate the season start and end.

# List of DSV
dl1 <- list.files("Z:/Late_blight/lbdsveconu/", recursive = T, pattern = "mask.*new.*tif$",
                  full.names = T)
# List of seasons for subsetting
seasn_cat <- unique(grwson$NOTE)

grwson <- project(grwson, y=crs(rast(dl1)))

thinoutf <- function(x, y) {
  # Year rast
  ras1 = rast(y)
  print(ras1)
  # Subset the shapefile
  seas_shp = grwson[grwson$NOTE == x,]
  seas_shp = project(seas_shp, y=crs(ras1))
  # String for name
  nm1 = paste0("seas_",gsub(":","",gsub(" ","_", x)),"_", nrow(seas_shp))
  # Year for name
  yr1 = str_sub(y, 27, 30)
  # Crop & mask
  seas_ras = crop(ras1, seas_shp)
  seas_ras = mask(seas_ras, seas_shp)
  # Days of season
  seas_st = unique(seas_shp$season_st)
  seas_en = unique(seas_shp$season_end)
  # Dates of season, possibly assign as time() later
  day_sst = (as.Date(paste0(yr1,"/01/01"))-1)+days(seas_st)
  day_sen = (as.Date(paste0(yr1,"/01/01"))-1)+days(seas_en)
  # Subset 
  values(seas_ras[[c(1:seas_st,seas_en:365)]]) = NA
  # Describe layers used in subset by assigning as name
  names(seas_ras) = paste0("l_",seq(seas_st, seas_en, by = 1))
  writeRaster(seas_ras, paste0("Z:/Late_blight/lbdsveconu/seasonzone/",nm1,
                               "_",yr1,".tif"))
  #return(seas_ras)
}


# This uses the two files, allmax & allmin.
# They are the 10 year means of MAXT and MINT for 2010-2019.
dl1 <- list.files("Z:/Late_blight/temp/", pattern = "allm", full.names = T)
dl1 <- dl1[c(1,3)]
grwson <- project(grwson, y=crs(rast(dl1[1])))

for (i in seasn_cat) {
  for (j in 1:length(dl1)){
    nm1 = paste0("seas_",gsub(":","",gsub(" ","_", i)),"_")
    # Year for name
    yr1 = str_sub(dl1[j], 24, 26)
    r1 = rast(dl1[j])
    r1 = crop(r1, grwson[grwson$NOTE == i,])
    r1 = mask(r1, grwson[grwson$NOTE == i,])
    writeRaster(r1, paste0("Z:/Late_blight/temp/",nm1,"_",yr1,".tif"))
  }
}

sl1 <- list.files("Z:/Late_blight/temp/", pattern = "seas", full.names = T)
sl1 <- sl1[-1]

seasns <- unique(gsub("__m.*","",gsub("seas_","",str_split_i(sl1, "/", 4))))

tempselcf <- function(x) {
  minsel = sl1[str_detect(sl1, x)][str_detect(sl1[str_detect(sl1, x)], "min")]
  print(minsel)
  maxsel = sl1[str_detect(sl1, x)][str_detect(sl1[str_detect(sl1, x)], "max")]
  print(maxsel)
  rmin = rast(minsel)
  rmax = rast(maxsel)
  minthrs = ifel(rmin > 10, 1, 0)
  maxthrs = ifel(rmax > 29, 0, 1)
  grwthrs = minthrs * maxthrs
  writeRaster(grwthrs, paste0("Z:/Late_blight/lbdsveconu/seasonzone/",x,"_thrsold.tif"))
}

for (i in seasns[10:15]) {
  tempselcf(i)
}

###
l1 <- list.files("Z:/Late_blight/temp/seas_dsv/", pattern="_a", full.names = T)

r1 <- rast(l1)

delavg <- list()

for (i in 1:365) {
  delavg <- app(del1[[seq(i, nlyr(del1), by=365)]], median)
}
