# Script name: moon_pred
# 
# Purpose of script: calculate the moon illumination
# 
# Author: Carlos Linares
# 
# Date Created: 07/29/2024
# 
# Email: carlosgarcialina@u.boisestate.edu
# 
# ---------------------------
#   
#   Notes: early in the analysis we calculated moon illumination using suncal and lunar. 11/21/2024 as recomended by Kayba we will use the moonlit. packae 


# inputs ------------------------------------------------------------------
# data_for_analysis/sites_coordinates.csv -> file with the coordinates for the field sites.
# data_for_analysis/prep_for_glmm/bat_combined.csv -> data base with the date and times


#script to calculate the average and median moon phase for use a predictor in a population model
#we will be using the suncalc package


# libraries ---------------------------------------------------------------



library(suncalc)
library(sf)
library(ggplot2)
library(lunar)
library(tidyverse)
library(moonlit)
library(beepr)


#load environment 

load('working_env/moon_pred.RData')


#we need the coordinates for the sites

sts<-read.csv("data_for_analysis/sites_coordinates.csv")

# load bat data 2021 for dates. 

kpro_2021_bat<-read.csv('data_for_analysis/kpro2021_v1.csv')
bm<-read.csv('data_for_analysis/data_glmm/bat_counts.csv', header = T) # bat counts by jday
# unique(bm$site)
#get dates from bat2021

datesite<-bm %>% select(site,jday)
datesite$date<-as_date(bm$jday, origin= "2021-01-01")
datesite$jday<-NULL

# generate dates from 2021 to 2023 just the summer months. 

start_date <- as.Date("2020-05-01")
end_date <- as.Date("2023-09-30")

# Generate the sequence of dates
all_dates <- seq.Date(from = start_date, to = end_date, by = "day")
# Filter the dates to include only May to September
dates<- all_dates[format(all_dates, "%m") %in% c("05", "06", "07", "08", "09")]
dates<-as.data.frame(dates)
colnames(dates)[1]<-"date"


# merge tdates# merge the dates with the site and the coordinates. 9180 obs because is 15 times the dates object. 

t1<-merge(dates, sts, by=NULL) 

# moon illumination, times, and pos ----------------------

t1<-left_join(datesite,sts, by = "site")
# t1<-rename(t1, date=noche) # rename noche as date but remember this when rejoining the other data.
t1<-t1 %>% select(date,lat,lon)

moon_pos<- unique(getMoonPosition(data = t1))
moon_pos$date<-as.Date(moon_pos$date)


t1$date<-as.Date(t1$date)
moon_times<-unique(getMoonTimes(data = t1))

t2 <- as.Date(unique(t1[, 1]), tz = "America/Denver")
moon_illumination <- getMoonIllumination(date = t2)


moon_pred <- left_join(left_join(moon_pos, moon_times), moon_illumination)


moon_pred<-moon_pred %>% select(date, lat, lon, altitude, rise,set,phase, fraction, parallacticAngle, angle)

moon_pred$above_horizon <- moon_pred$altitude > 0


#lunar illumination calculated with lunar package


attr(moon_pred$date, "tzone") <- "UTC" # make them UTM so the calculation is ok
# myshift <- as.numeric(format(moon_pred$Datetime, "%H")) - 12 # read https://stackoverflow.com/questions/71757462/calculate-lunar-illumination-using-the-lunar-package-in-r
moon_pred$l.illum<-lunar.illumination(moon_pred$noche)



# make moon pred 0 if not above the horizon.

moon.adj<-moon_pred %>% mutate(
  phase = ifelse(above_horizon==FALSE,0,phase),
  fraction= ifelse(above_horizon==FALSE,0,fraction),
  l.illum= ifelse(above_horizon==FALSE,0,l.illum)
)


colnames(moon_pred)[1]<-"date"




write.csv(moon_pred,file = 'data_for_analysis/moon_pred.csv', row.names = F) # ran once in case of needing to rewrite the data



# moonlit -----------------------------------------------------------------

# # install package
# install.packages("devtools")
# library(devtools)

# 
# install_github("msmielak/moonlit")

#load the moonlit library

# calculateMoonlightIntensity(lat, lon, date, e) Function requires dates and lat long data. 

# site coordinates
sts<-read.csv("data_for_analysis/sites_coordinates.csv")

# bat data set time and dates 2021-2023

bat_combined<-read_csv(file = 'data_for_analysis/prep_for_glmm/bat_combined.csv') # 70 rows are missing the time for the date time column. see below. 

#--- 1pm times -----
  # here we are exploring if there are any wrong dates in the original data set because there are wrong dates in the moon intensity. For some reason there are recordings made at 1 pm at one site but these are all noise and will be filter out at the end. 
  
rows_13pm <- bat_combined %>%
  filter(format(as.POSIXct(date_time, tz = "UTC"), "%H") == "13")

#------


site.date<- bat_combined %>% select(site, date_time)

sidalo<-left_join(site.date, sts, by="site")

# sidalo$date_time_prsed<- ymd_hms(sidalo$date_time, tz="UTC") #keep it UTC.

# failed_rows <- which(is.na(sidalo$date_time_prsed))
# failed_dates <- sidalo[failed_rows, ]

sidalo$date_time<- force_tz(sidalo$date_time, tzone = "America/Denver")

summary(sidalo)




# calculate moon intensity
# we use the combine data set to extract the bat calls time and dates for the 2021-2023 time frame. 

moon.int <- calculateMoonlightIntensity(
  sidalo$lat,
  sidalo$lon,
  sidalo$date_time ,
  e = 0.16
)

summary(moon.int) # yes there is no NAs now we need to merge it with sidalo. 



s.sidalo <- sidalo[sample(nrow(t), size = 5000, replace = TRUE), ]
s.sidalo<- head(sidalo)
saveRDS(sidalo, file = "data_for_analysis/moon_pred/sidalo.csv")

moon.stat <- calculateMoonlightStatistics(
  s.sidalo$lat,
  s.sidalo$lon,
  s.sidalo$date_time,
  e = 0.16,
  t = "1 hour",
  timezone = "America/Denver"
)
beep()


Rprof("profile.out")
moon.stat <- calculateMoonlightStatistics(
  s.sidalo$lat,
  s.sidalo$lon,
  as.POSIXct(s.sidalo$date_time, tz = "America/Denver"),
  e = 0.16,
  t = "1 hour",
  timezone = "America/Denver"
)
Rprof(NULL)
summaryRprof("profile.out")

stats <- apply(moon.int, 1, function(row) {
  calculateMoonlightStatistics(
    lat = as.numeric(row["lat"]),
    lon = as.numeric(row["lon"]),
    date = as.POSIXct(row["date"], tz = "America/Denver"),
    e = 0.26,
    t = "15 mins",
    timezone = "America/Denver"
  )
})
beep()

stats <- calculateMoonlightStatistics(sidalo$lat, sidalo$lon, min(sidalo$date_time), e = 0.16, t = "1 hour", timezone = "UTC")


# parallel ----------------------------------------------------------------
library(parallel)

# # Preprocess date_time column
# sidalo$date_time <- as.POSIXct(sidalo$date_time, tz = "UTC")

# Define function to process chunks
process_chunk <- function(lat, lon, date_time) {
  calculateMoonlightStatistics(lat, lon, date_time, e = 0.16, t = "1 hour", timezone = "UTC")
}

# Split into chunks
num_cores <- detectCores() - 2
rows_per_chunk <- ceiling(nrow(sidalo) / num_cores)
chunks <- split(1:nrow(sidalo), ceiling(seq_along(1:nrow(sidalo)) / rows_per_chunk))

# Create a cluster
cl <- makeCluster(num_cores)

# Export necessary variables and functions to the cluster
clusterExport(cl, c("sidalo", "process_chunk", "calculateMoonlightStatistics"))

# Apply function in parallel
results <- parLapply(cl, chunks, function(idx) {
  process_chunk(
    sidalo$lat[idx],
    sidalo$lon[idx],
    sidalo$date_time[idx]
  )
})

# Stop the cluster
stopCluster(cl)

# Combine results
moon.stat <- do.call(rbind, results)

#----


# outputs -----------------------------------------------------------------

write.csv(moon.int, file = 'data_for_analysis/moon_pred_out/moon.ill.csv', row.names = F) # data 


# save --------------------------------------------------------------------

save.image(file = "working_env/moon_pred.RData")



#-----------------plots-----------------

#coordinates 
library(sf)
library(ggplot2)

# Create a data frame with the coordinates
coords <- data.frame(lat = 43.5418, lon = -113.7331)
points_sf <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)

# Plot the point on a map
ggplot() +
  geom_sf(data = points_sf, color = "red", size = 3) +
  borders("state") +
  coord_sf(xlim = c(-114, -112), ylim = c(42, 44), expand = FALSE) +
  theme_minimal()


ggplot(moon.adj, aes(x=altitude, y=phase))+
  geom_point()+
  geom_vline(xintercept = 0)

ggplot(moon.adj, aes(x=altitude, y=l.illum))+
  geom_point()+
  geom_vline(xintercept = 0)

hist(moon_pred$phase)
hist(moon_pred$fraction)
hist(moon_pred$l.illum)

plot(moon_pred$phase ~ moon_pred$l.illum)# no linear relationship between moon illumination and moon phase. 


# Create a new column indicating moon above horizon
moon_pos$above_horizon <- moon_pos$altitude > 0

# Create the plot
ggplot(moon_pos, aes(x = yday(moon_pos$date), y = altitude, color = above_horizon)) +
  geom_line() +
  # scale_x_datetime(breaks = date_breaks("12 hours")) +  # Adjust breaks as needed
  labs(x = "Time", y = "Moon Altitude", color = "Moon Above Horizon?") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray")  # Optional: Show threshold lin



# examples ----------------------------------------------------------------
# 
# Here is code we don't need or use anymore but that helped solve a problem in the past
# 
# Sample of sites 

# I coded this section to learn how ti fix the issue withe the time zone changing the time when passing it to the moonlit package functions. Once I learned I moved it to the end of the code. It should not be needed to run again. It should be used as an example.
t_sampled <- t[sample(nrow(t), size = 10, replace = TRUE), ]

date_time_utc <- ymd_hms(t_sampled$date_time, tz = "UTC") # convert to UTC 

# Change the time zone to "America/Denver", but keep the exact same clock time
# To keep the same time, we use `force_tz`.
date_time_denver <- force_tz(date_time_utc, tzone = "America/Denver") # this forces the time zone to MDT

# Show the result
date_time_denver


# trash

# 
# # Generate the sequence of dates
# all_dates <- seq.Date(from = start_date, to = end_date, by = "day")
# 
# # Filter the dates to include only May to September
# dates <- all_dates[format(all_dates, "%m") %in% c("05", "06", "07", "08", "09")]
# dates <- as.data.frame(dates)
# colnames(dates)[1] <- "date"
# 
# # Generate times for each date from 6 PM to 6 AM the next day
# time_sequences <- lapply(dates$date, function(date) {
#   start_time <- as.POSIXct(paste(date, "18:00:00"), tz = "America/Denver")
#   end_time <- as.POSIXct(paste(date + 1, "06:00:00"), tz = "America/Denver")
#   seq(from = start_time, to = end_time, by = "min")
# })
# 
# # Combine the list into a single data frame
# times <- do.call(rbind, lapply(seq_along(time_sequences), function(i) {
#   data.frame(date = dates$date[i], time = time_sequences[[i]])
# }))
# 
# head(times)
# 
# # merge the dates with the site and the coordinates. 
# 
# t1<-merge(times, sts, by=NULL)

