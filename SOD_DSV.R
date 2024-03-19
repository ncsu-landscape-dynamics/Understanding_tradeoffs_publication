library(terra)
library(tigris)

# Calculations to set up the DSV scores.
# Formula from Meentemeyer et al. 2011
sodweather <- function(x){
  y = -0.066 + 0.056*x - 0.0036*(x - 15)^2 - 0.0003*(x - 15)^3
}

# Take some potential temps to calculate the weather coefficient
sodv2 <- lapply(seq(-2, 32, by=0.25), sodweather)

sodvf2 <- data.frame(t = seq(-2,32,by=0.25), dv = unlist(sodv2))

# The coefficient has a hump at small, positive temps.
# Remove the hump and use a line.
xs <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75)
yx <- c(0.1365, 0.12749, 0.11969, 0.11307, 0.1076, 0.10326,0.1,0.0978, 0.0967, 0.0966,0.0974,0.099255)
m1 <- lm(yx ~ xs)
summary(m1)

xs <- c(0,3.5)
yx <- c(0, 0.115)
m1 <- lm(yx ~ xs)
summary(m1)

# Simplify the SOD data frame, since all of the T < 0 have no impact on DSV.
sodvf2[sodvf2$t < 0, "dv"] <- 0 
sodvf2 <- sodvf2[sodvf2$dv > 0.0,]
sodvf2[1:15, "dv"] <- sodvf2[1:15, "t"] * 0.0338

# Calculate the thresholds
thresh1 <- 0.2*max(sodvf2$dv) # 0.1853247
thresh2 <- 0.4*max(sodvf2$dv) # 0.3706494
thresh3 <- 0.6*max(sodvf2$dv) # 0.5559741
thresh4 <- 0.8*max(sodvf2$dv) # 0.7412988

# DSV logic
dsv <- function(precip, tmean) {
  prec_val = ifel(precip >= 2, 1, 0)
  
  t_val = ifel(tmean > 13.75 & tmean < 24.75, 4, 
        ifel(tmean > 10.25 & tmean <= 13.75 | tmean >=24.75 & tmean < 26.5, 3,
             ifel(tmean <= 10.25 & tmean > 6.75 | tmean >= 26.5 & tmean < 28.0, 2,
                 ifel(tmean <= 6.75 | tmean >= 28.0 , 1, 0))))
  
  dsv_val = prec_val * t_val
  
  return(dsv_val)
}

# Geographical focus for SOD
co1 <- counties(state = c("CA","OR"))
co1 <- co1[co1$NAME %in% c("Coos","Curry","Douglas","Josephine"),]#,"Del Norte"),]

txl <- list.files("Q:/Shared drives/Data/Original/Daymet/tmax", full.names = T)
tnl <- list.files("Q:/Shared drives/Data/Original/Daymet/tmin", full.names = T)
prl <- list.files("Q:/Shared drives/Data/Original/Daymet/precip/", full.names = T)
txl <- txl[30:41]
tnl <- tnl[30:41]
prl <- prl[30:41]

# Set the focus to the right projection
r0 <- rast(txl[1])
co1 <- vect(co1)
co1 <- project(co1, y=crs(r0))

# Function for use in for loop
dsvf <- function(i){
  mx = rast(txl[i])
  mn = rast(tnl[i])
  pr = rast(prl[i])
  mx = crop(mx, co1, mask=T)
  mn = crop(mn, co1, mask=T)
  precip = crop(pr, co1, mask=T)
  
  tmean = (mx+mn)/2
  
  codsv = dsv(precip, tmean)
  return(codsv)
}

tdsv <- list()

tdsv <- for(i in 1:length(txl)){
  dsvf(i)
}

# Project the DSVs to EPSG 3857
for(i in 1:12){
indq = seq(1, 4380, by=365)
temr = r1[[indq[i]:(indq[i]+364)]]
r10[i] = project(temr, crs(r2))
}

# The resolution was wrong after projecting. 
# Used one of the weather rasters as a template for resampling in order
# to have right resolution.
res1 = list()
for(i in 1:12){
  temr = crop(rast(r10[i]), ext(r3), mask=T)
  res1[i] = resample(temr, r3, method="near")
}

# Fix some values that resampled as decimal
rm1 <- matrix(c(0,0.5,0,0.5,1,1,1,1.5,1,1.5,2,2,2,2.5,2,2.5,3,3,3,3.5,3,3.5,4,4),
              ncol=3, byrow = T)

temr2 = list()

for(i in 1:length(res1)){
  temr2[i] = classify(temr, rm1, right=T)
}

# Cumulative sum of the daily values in DSV.
# Cumulative counts of DSV/18
yseq <- 2009:2020

for(i in 1:12){
  l1 <- list.files(paste0("Q:/My Drive/popblightdata/newareas/DSV/",yseq[i],"/daily/"), full.names = T)
  l2 <- l1[c(1, 112, 223, 300, 311, 322, 333, 344, 355, 2, 13, 24, 35, 46, 57, 68, 79, 90, 101, 113, 124, 135, 146, 157, 168, 179, 190, 201, 212, 224, 235, 246, 257, 268, 279, 290, 297:299, 301:310, 312:321, 323:332, 334:343, 345:354, 356:365, 3:12, 14:23, 25:34, 36:45, 47:56, 58:67, 69:78, 80:89, 91:100, 102:111, 114:123, 125:134, 136:145, 147:156, 158:167, 169:178, 180:189, 191:200, 202:211, 213:222, 225:234, 236:245, 247:256, 258:267, 269:278, 280:289, 291:296)]
  cd1 <- app(rast(l1), cumsum)
  cd2 <- app(cd1, function(x) floor(x/18))
  for(j in 1:365){
    jpeg(paste0("Q:/My Drive/DSV/",yseq[i],"/cumulthreshold/cuthresh",yseq[i],"_",j,".jpg"))
    plot(cd2[[j]])
    dev.off()
  }
  writeRaster(cd2, paste0("Q:/My Drive/DSV/",yseq[i],"/cuthresh",yseq[i],"_",j,".tif"))
}
