library(terra)

rootdir <- "z:/"
# I use Z:/

# 2021, 2022 had very small counts.
# These are the rasters for the data as they come from Plant Aid DB
padb21 <- rast(file.path(rootdir, "Late_blight/Manuscript_1_Data/simulation/LB/inputs/infection/infection_2021.tif"))
# 6 locations x 2 infections
padb21 <- rast(file.path(rootdir, "Late_blight/Manuscript_1_Data/simulation/LB/inputs/infection/infection_2022.tif"))
# 1 loc x 1 infection

# Created new rasters with infections at random in place of those.
# These are the files that had the augmented counts.
y21 <- rast(file.path(rootdir, "Late_blight/Manuscript_1_Data/simulation/LB/inputs/infection/infection_2021_rev_10locs.tif"))
y22 <- rast(file.path(rootdir, "Late_blight/Manuscript_1_Data/simulation/LB/inputs/infection/infection_2022_rev_10locs.tif"))

# These each had 10 locs with 2 infections each: 20 infections

# Each of these rasters had just one layer, not multiple layers for different time steps.

# Then they were run through steps 1 - 3 from the 
https://github.com/ncsu-landscape-dynamics/Understanding_tradeoffs_publication/tree/main/simulation_study/LB

# Step 3 produces mean, median, etc. and they each have one layer. 

# I used the means in the calibrations. Those were
new21 <- rast(file.path(rootdir, "Late_blight/Manuscript_1_Data/simulation/LB/outputs21/locs10x2/upwdrev/rasts/pops_mean_Year_2021.tif"))
new22 <- rast(file.path(rootdir, "Late_blight/Manuscript_1_Data/simulation/LB/outputs22/locs10x2/upwdrev/rasts/pops_mean_Year_2022.tif"))
              
# No time steps

# Their histograms
2021:
  layer value   count
1      1     0 2274225
2      1     1      47
3      1     2      19
4      1     3       6
5      1     4       5
6      1     5       2
7      1     7       2
8      1     8       1
9      1     9       5
10     1    11       1
11     1    14       1
12     1    25       1
13     1    29       1
# total infections: 279

2022:
  layer value   count
1      1     0 2274201
2      1     1      68
3      1     2      26
4      1     3      10
5      1     4       1
6      1     5       2
7      1     6       1
8      1     7       1
9      1     8       1
10     1     9       1
11     1    12       1
12     1    15       1
13     1    20       1
14     1    29       1
# total infections: 270

# I was using pops_mean_Year_202* in the all of the simulations as what I thought
# was the starting infection file because it had way more infections that the 
# original data for each year.

# All the data reductions (10%, etc.) start with the files on lines 28 and 29. 

# Then reproduction was all 0 and to fix that, I thought I was supposed to just use those same
# rasters for start, just with the reduced zeros, where the N of 0 matched N of
# locations with infection. 
# These are the reduced zero versions:
y21_redc <- rast(file.path(rootdir, "Late_blight/Manuscript_1_Data/simulation/LB/outputs21/locs10x2/upwdrev/rasts/mean2021redczero.tif"))
y22_redc <- rast(file.path(rootdir, "Late_blight/Manuscript_1_Data/simulation/LB/outputs22/locs10x2/upwdrev/rasts/mean2022redczero.tif"))

# So those have less zeros, but same N of infections as line 28 and 29



