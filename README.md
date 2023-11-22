# Late blight and sudden oak death model comparison paper workflow repository

## Data Preparation
PoPS uses multiple inputs for infections, weather, and settings.

### Late blight
To forecast late blight, infection data came from a data set that combined infections in the USA Blight and National Plant Diagnostic Network. Those data were in a tabular format and were assigned to a raster based on county and state. Daymet's min. and max. temperatures, vapor pressure, and daylength rasters were processed in R to calculate an hourly mean temperature and relative humidity. The mean temperature and relative humidity were then used to calculate a weather coefficient for use in PoPS.

Late blight was also modeled with a weather risk model called a degree severity value. This model used the same mean temperature and relative humdity estimates used in PoPS. The DSV uses if-then logic to assign a score from 0-4, 0 equating to no risk and 4, infection imminent. 

### SOD
To forecast SOD with PoPS, we used infection data collected by the Oregon Department of Forestry and U.S. Forest Service in Coos County, Oregon and the 3 surrounding Oregon counties. Weather inputs were Daymet temperatures and precipition. Host and total populations data came from the ODF and USFS.

SOD was modeled with a DSV that used the same range of scores as that of late blight. The same mean temperature and precipitation data as those for late blight were used in the if-then logic for the DSV for SOD.
