# Late blight and sudden oak death model comparison paper workflow repository

## Data requirements
PoPS uses multiple inputs for infections, weather, and settings. We relied on Daymet for weather, Cropland Data Layer for host and all populations, and our infection data for late blight came usablight.org, The National Plant Diagnostic Network, and Twitter. 

### Late blight
The data were in a tabular format and were assigned to a raster based on county and state. Daymet's min. and max. temperatures, vapor pressure, and daylength rasters were processed in R to calculate an hourly mean temperature and relative humidity. The mean temperature and relative humidity were then used to calculate a weather coefficient for use in PoPS.

Late blight was also modeled with a weather risk model, "WRM". We started this project calling it "DSV", so if that acronym turns up, it's referring to the same model/values. This model used the same mean temperature and relative humdity estimates used in PoPS. The WRM uses if-then logic to assign a score from 0-4, 0 equating to no risk and 4, infection imminent. 

The SDM used several predictor inputs. The flexsdm library will do almost all of the work in the script. It will take some time to get the predictors' original files in the proper directories before running the script. See the Methods section for more specifics on which predictors were used. The selection is up to the author.

### SOD
To forecast SOD with PoPS, we used infection data collected by the Oregon Department of Forestry and U.S. Forest Service in Coos County, Oregon and the 3 surrounding Oregon counties. Weather inputs were Daymet temperatures and precipitation. Host and total populations data came from the LEMMA product from Oregon State.

SOD was modeled with a WRM that used the same range of scores as that of late blight, however with the thresholds at different precipitation and temperature values than those of late blight. Mean temperature, calculated from Daymet temperature data, and Daymet precipitation were used in the if-then logic for the WRM for SOD. 

The SDM for SOD is the same as PoPS.

### Guide
This is a sort of map to indicate how the scripts above are used to create the results for the three models. If a script above doesn't appear in the guide, then it was probably used to create a figure or do some ancillary calculations, like `SOD_WRM_temp_curve_threshold_estimation.R`.
```
*WRM*
	_late blight_:
		late_blight_weather_time_prep.R: calculates the required inputs for ...
		
		late_blight_temp_hour.R; late_blight_WRM_calc.R >>> late_blight_wc_and_WRM.R 
			
	_SOD_:
		SOD_wrm_creation.R
		

*SDM*
	_late blight_:
		run_sdmLB.R
		
		SDM_pred_reduction_loop.R
		
		rerun run_sdmLB with the reduced set of predictors
		
	_SOD_:
		run_sdmSOD.R
		
		SDM_pred_reduction_loop.R
		
		rerun run_sdmSOD with the reduced set of predictors
		
		
*PoPS*
	_late blight_:
		late_blight_infection_creation.R
		late_blight_host_nonhost_creation.R
		[weather coeffient from late_blight_wc_and_WRM.R]
		config
		calibration
		
	_SOD_:
		
```
