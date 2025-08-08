# Late blight and sudden oak death model comparison paper workflow repository

## Data Preparation
PoPS uses multiple inputs for infections, weather, and settings.

### Late blight
To forecast late blight, infection data came from a data set that combined infections in the USA Blight and National Plant Diagnostic Network. Those data were in a tabular format and were assigned to a raster based on county and state. Daymet's min. and max. temperatures, vapor pressure, and daylength rasters were processed in R to calculate an hourly mean temperature and relative humidity. The mean temperature and relative humidity were then used to calculate a weather coefficient for use in PoPS.

Late blight was also modeled with a weather risk model, "WRM". This model used the same mean temperature and relative humdity estimates used in PoPS. The WRM uses if-then logic to assign a score from 0-4, 0 equating to no risk and 4, infection imminent. 

### SOD
To forecast SOD with PoPS, we used infection data collected by the Oregon Department of Forestry and U.S. Forest Service in Coos County, Oregon and the 3 surrounding Oregon counties. Weather inputs were Daymet temperatures and precipitation. Host and total populations data came from the ODF and USFS.

SOD was modeled with a DSV that used the same range of scores as that of late blight, however with the thresholds at different precipitation and temperature values than those of late blight. Mean temperature, calculated from Daymet temperature data, and Daymet precipitation were used in the if-then logic for the WRM for SOD.

### Guide
This is a sort of map to indicate how the scripts above are used to create the results for the three models. If a script above doesn't appear in the guide, then it was probably used to create a figure or do some ancillary calculations.

WRM
	late blight:
		late_blight_weather_time_prep.R: calculates the required inputs for 2
		
		late_blight_temp_hour.R; late_blight_WRM_calc.R >>> late_blight_wc_and_WRM.R 
			
	SOD:
		SOD_wrm_creation.R
		

SDM
	late blight:
		run_sdmLB.R
		
		SDM_pred_reduction_loop.R
		
		rerun run_sdmLB with the reduced set of predictors
		
	SOD:
		run_sdmSOD.R
		
		SDM_pred_reduction_loop.R
		
		rerun run_sdmSOD with the reduced set of predictors
		
		
PoPS
	late blight:
		late_blight_infection_creation.R
		late_blight_host_nonhost_creation.R
		[weather coeffient from late_blight_wc_and_WRM.R]
		config
		calibration
		
	SOD:
		
