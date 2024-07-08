library(data.table)
library(terra)
library(tidyverse)

# This works on HPC to take all cropped regions years and subset to the days
# ARBITRARILY PICKED for growing season. See the tomato_growing_seasons.R script
l1 <- list.files("/rs1/researchers/c/cmjone25/Late_blight/temp/seas_dsv", full.names=T)

seasdir <- read.csv("/rs1/researchers/c/cmjone25/Late_blight/relatedtables/seasondatesdirs.csv")

pixf <- function(x) {
  print(x)
  nm1 = str_sub(x, 55, 55)
  s_day = as.numeric(seasdir$season_st[seasdir$dirname == nm1])
  e_day = as.numeric(seasdir$season_end[seasdir$dirname == nm1])
  nyr = 12
  rast_p = list.files(x, full.names=T)
  print(paste0('rast_p ', rast_p))
  rdel1 = rast(rast_p[1])
  rdel1 = rdel1[[50]]
  nmcl = ncell(rdel1[!is.na(rdel1)])
  
  res_mat = matrix(nrow = nmcl * (e_day - s_day + 1), ncol = nyr + 1)
  
  row_ind = 1
  
  for (day in s_day:e_day) {
    for (year in 1:nyr) {
      rf1 = rast_p[year]
      r1 = rast(rf1)
      
      day_lay = r1[[day]]
      
      pixvalues = values(day_lay, na.rm=T)
      
      res_mat[row_ind:(row_ind + length(pixvalues) - 1), year] = pixvalues
      
      res_mat[row_ind:(row_ind + length(pixvalues) - 1), nyr + 1] = day
    }
    row_ind = row_ind + nmcl
  }
  colnames(res_mat) = c(paste0("year_", 1:nyr), "day")
  write.csv(res_mat, paste0("/rs1/researchers/c/cmjone25/Late_blight/temp/seas_dsv/",nm1 ,"_pix_sum.csv", row.names=F))
}

for (i in l1) {
  pixf(i)
}

===
# DOESN'T WORK. EVEN a month is too much data for ggplot to manage.
# Attempt to get each month of each cropped subset to create a strip plot

dtmo <- dt[moy == x,]
dtmlo <- melt(dtmo, measure.vars = patterns("^year_"), value.name = "pval")
dtml0[, c("day", "moy") := .(rep(day, times = .N), rep(moy, times = .N)), by = .(variable)]
dtmlo[, variable := NULL]

ggplot() + geom_jitter
# DOESN'T WORK. EVEN a month is too much data for ggplot to manage.
===

#l1 <- list.files("Z:/Late_blight/lbdsveconu/appscores/csv", pattern = "sub_", full.names=T)
l1 <- list.files("Z:/Late_blight/temp/seas_dsv/", full.names = T, pattern="csv$")
l1 <- gtools::mixedsort(l1)

chron_order <- month.name
# Manually create a data frame for summary() values of each month from the huge data.table
# Use in a loop
year_sum <- function(z, x) {
  j1 = z[moy == x,]
  j2 = j1 %>% pivot_longer(cols = starts_with("y"), names_to = "year", values_to = "values")
  j2sum = summary(j2$values)
  temdf = data.frame(v = rbind(j2sum[[1]], j2sum[[2]], j2sum[[3]], j2sum[[5]], j2sum[[6]]))
  return(temdf)
}

# Data frame with growing season description ("NOTE"), first day of season, last day of season, directory name,
# and the region designation
vdf <- read.csv("Z:/Late_blight/relatedtables/seasondatesdirs.csv")

# Get summary of values. Create a boxplot from the summary. Kind of short-circuits
# the ggplot boxplot. This could be useful in other places. Requires the year_sum()
sumgraphf <- function(x) {
  #dt1 <- fread("Z:/Late_blight/lbdsveconu/appscores/csv/sub_9.csv")
  tmpdt = fread(l1[x])
  print(x)
  print(x)
  
  nam1 = str_sub(l1[x], 30, 30)
  print(vdf[vdf$dirname == nam1,])
  
  # Classify all the days into months
  tmpdt[, moy := cut(day, breaks=c(0, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365), include.lowest = T, labels= month.name, right = F)]
  sumdf = data.frame(NA)
  
  ## The list of zones to supply a name
  #zlab = paste0("z", rep(3:12, each = 2), c("a", "b"))
  #zlab = zlab[1:19] 
  #nlab = gsub("z","",zlab)
  
  tempdf = list()
  
  for (j in month.name) {
    tempdf[j] = year_sum(tmpdt, j)
  }
  
  yr_sum_df = do.call(cbind, tempdf)
  yr_sum_df = as.data.frame(yr_sum_df)
  rownames(yr_sum_df) = c("ymin", "lower", "mid", "upper", "ymax")
  yr_sum_df = t(yr_sum_df)
  
  newlabs = c("January" = "Jan", "February" = "Feb", "March" = "Mar", "April" = "Apr", 
               "May" = "May", "June" = "Jun", "July" = "Jul", "August" = "Aug", 
               "September" = "Sep", "October" = "Oct", "November" = "Nov", 
               "December" = "Dec")
  
  reg_name = vdf$region[vdf$dirname == nam1]
  
  grow_seas1 = format(as_date(days(vdf$season_st[vdf$dirname == nam1])), "%d %B")
  
  grow_end = ifelse(vdf$season_end[vdf$dirname == nam1] == 365, 364, vdf$season_end[vdf$dirname == nam1])
  print(grow_end)
  grow_seas2 = format(as_date(days(grow_end)), "%d %B")
  
  p1 = ggplot(yr_sum_df, aes(x = factor(rownames(yr_sum_df), levels = chron_order))) + 
    geom_boxplot(aes(ymin = ymin, lower = lower, middle = mid, upper = upper, ymax = upper), 
    stat = "identity") + geom_point(aes(y = ymax)) + scale_y_continuous(limits = c(0, 4)) + 
    labs(x = "Month", y = "DSV", title = paste0("Boxplot of DSV Values in Region ",
         reg_name, ", 2009 - 2020\nGrowing Season Approximately: ", grow_seas1,
         " - ", grow_seas2)) +
    theme(panel.grid.major.x= element_blank(), panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(), axis.title.x = element_text(size = "14",
    face = "bold"), plot.title = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = "16", face = "bold"),
    axis.text.x = element_text(size = "11", face = "bold"),
    axis.text.y = element_text(size = "14", face = "bold") ) + scale_x_discrete(
      labels = newlabs)
  
  p1
  
  ggsave(paste0("Z:/Late_blight/lbdsveconu/appscores/csv/region",
                reg_name,"_plot.png"), width = 820, height = 620, units="px", dpi = 96, plot = p1)

}

#for (i in 1:length(l1)){
for (i in c(3,4,5,8,10:15)){
  sumgraphf(i)
}

tmpdel <- fread(l1[1])
tmpdel[, moy := cut(day, breaks=c(0, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365), include.lowest = T, labels= month.name, right = F)]
sumdf <- data.frame(NA)

p1 = ggplot(yrsumdf, aes(x = factor(rownames(yrsumdf), levels = chron_order))) + 
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = mid, upper = upper, ymax = upper), 
               stat = "identity") + geom_point(aes(y = ymax)) + scale_y_continuous(limits = c(0, 4)) + 
  labs(x = "Month", y = "DSV", title = paste0("Boxplot of DSV Values 2009 - 2020")) +
  theme(panel.grid.major.x= element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), axis.title.x = element_text(size = "14",
                                                                          face = "bold"), plot.title = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = "16", face = "bold"),
        axis.text.x = element_text(size = "11", face = "bold"),
        axis.text.y = element_text(size = "14", face = "bold") ) + scale_x_discrete(
          labels = newlabs)

