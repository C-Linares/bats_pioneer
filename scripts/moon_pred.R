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
#   Notes: early in the analysis we calculated moon illumination using suncal and lunar. 11/21/2024 but as recomended by Cris Kayba we will use the moonlit. package 


# inputs ------------------------------------------------------------------
# data_for_analysis/sites_coordinates.csv -> file with the coordinates for the field sites.
# data_for_analysis/prep_for_glmm/bat_combined.csv -> data base with the date and times

# outputs -----------------------------------------------------------------
# fwrite(combined_data, file = "data_for_analysis/moon_pred/moon.stat.csv", row.names = FALSE)
# 
# fwrite(combined, file = "data_for_analysis/moon_pred/moon.int.csv", row.names = FALSE)

#Notes:
#script to calculate the average and median moon phase for use a predictor in models.
#we will be using the suncalc, lunar, and moonlit package.


# libraries ---------------------------------------------------------------


if (!require("pacman")) install.packages("pacman")

pacman::p_load(suncalc, sf, ggplot2, lunar, tidyverse, moonlit, beepr, data.table, parallel, future)
library(data.table)

# moonlit  recommended by Kayba to Jesse Barber.


#load environment 

load('working_env/moon_pred.RData') # this is outdated. I need to load the data from the data base.


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
moon_illumination <- getMoonIllumination(date = t2) # moon illumination as fraction of the moon 0(new moon) to 1(full moon) 


moon_pred <- left_join(left_join(moon_pos, moon_times), moon_illumination)


moon_pred<-moon_pred %>% select(date, lat, lon, altitude, rise,set,phase, fraction, parallacticAngle, angle)

moon_pred$above_horizon <- moon_pred$altitude > 0


#lunar illumination calculated with lunar package
# this package needs some kind of shift in the time zone. 
# The number of hours by which to shift the calculation of lunar phase. By default lunar phase is calculated at 12 noon UT.


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


#setting up moonlit. 
# # install package
# install.packages("devtools")
# library(devtools)

# 
# install_github("msmielak/moonlit")

#load the moonlit library
#
library(moonlit)

#site coordinates
sts<-read.csv("data_for_analysis/sites_coordinates.csv")

# bat data set time and dates 2021-2023

bat_combined<-read_csv(file = 'data_for_analysis/prep_for_glmm/bat_combined.csv') # 70 rows are missing the time for the date time column. see below. 

# filter noise out

bat_combined<- bat_combined %>% filter(sp !="Noise")

#--- 1pm times -----
  # here we are exploring if there are any wrong dates in the original data set because there are wrong dates in the moon intensity. For some reason there are recordings made at 1 pm at one site but these most likely is an error while setting the clock of the sm3 and these calls will be left out of the analysis. 
  
rows_13pm <- bat_combined %>%
  filter(format(as.POSIXct(date_time, tz = "UTC"), "%H") == "13")

#------ site date location 

# below we extract sites and dates from the combined data sets from 2021-2023. 

site.date<- bat_combined %>% select(site, date_time)

sidalo<-left_join(site.date, sts, by="site")

sidalo$date_time<- force_tz(sidalo$date_time, tzone = "America/Denver")

summary(sidalo)
str(sidalo)

# check 
unique(lubridate::year(sidalo$date_time)) # all years present
sum(is.na(sidalo$date_time)) # no NAs in date time


# write site and dates. 

dir.create("data_for_analysis/moon_pred", recursive = TRUE, showWarnings = FALSE)
fwrite(sidalo, file = "data_for_analysis/moon_pred/sidalo.csv", row.names = FALSE) # fwrite does better than write.csv for large files and is faster. 


# below is the non-reproducible way to get the results I need. 
sidalo<-fread("data_for_analysis/moon_pred/sidalo.csv") # read the file to check it is ok.

# Not suere i need the code here until line 201
# # extract date time for bm2
# nightsummary<-bm %>% select(site, noche) # get dates for summarized bat data from glmm_v2
# 
# # merge with the coordinates
# 
# sidalo<-left_join(nightsummary, sts, by="site") # merge with the coordinates
# str(sidalo) # check the data set.)
# sidalo$date_time<- as.POSIXct(sidalo$noche) # convert to date time 
# summary(sidalo)


# code suggestions by Kyle Shanon ---------------------------------------------------
library(future)
#
print(nrow(sidalo))

print("calculating intensity...")
system.time(moon.int <- calculateMoonlightIntensity(
  sidalo$lat,
  sidalo$lon,
  sidalo$date_time ,
  e = 0.16,
))


# Optionally, retrieve the result once it is completed
result <- value(moon_int_future)
print(result)

# this is the metadata for the moon.int data set.
metadata <- data.frame(
  night = "Logical, TRUE when sun is below the horizon",
  sunAltDegrees = "Solar altitude in degrees",
  moonlightModel = "Predicted moonlight illumination, relative to an 'average' full moon",
  twilightModel = "Predicted twilight illumination in lux",
  illumination = "Combined moon and twilight intensity, in lux",
  moonPhase = "Lunar phase as a numerical value (% of moon face illuminated)"
)
# Convert the metadata data frame to a data table while preserving variable names
metadata_dt <- as.data.table(t(metadata), keep.rownames = "Variable")
setnames(metadata_dt, "V1", "Description")

# Combine the metadata and moon.int data frames

t <- left_join(moon.int, sidalo, by = c("date" = "date_time", "lat" = "lat", "lon" = "lon"))

combined <- rbindlist(list(metadata_dt, moon.int), use.names = FALSE, fill = TRUE)

# There are some observations where the moon was below the horizon and the sun was barely out or out

moon.int.night<-moon.int %>% filter(night==TRUE)

# create year col

moon.int.night$year<-lubridate::year(moon.int.night$date)
unique(moon.int.night$year) # check years are there






# moon.stat ---------------------------------------------------------------




cl <- makeCluster(detectCores())
print("calculating parApply statistics")
moon.stat <- NULL
system.time(moon.stat <- parRapply(cl, sidalo, function(row) {
  library(moonlit)
  stats <- calculateMoonlightStatistics(
    as.double(row[3]),
    as.double(row[4]), row[2],
    e = 0.16,
    t = "1 hour",
    timezone = "America/Denver")
  return(stats)
}
))
stopCluster(cl)

moon.df <- do.call(rbind, moon.stat)

#make it a dataframe with data.table 
moon.df<-as.data.frame(moon.df)

#write it as csv
# fwrite(moon.df, file = "data_for_analysis/moon_pred/moon.stat.csv", row.names = FALSE)
# we read back the data set. 
# moon.df<-fread("data_for_analysis/moon_pred/moon.stat.csv")
# now we include metadata. 
metadata <- data.frame(
  sunset = "Time of sunset",
  sunrise = "Time of sunrise",
  meanMoonlightIntensity = "Mean value of modeled illumination for the night",
  minMoonlightIntensity = "Min value of modeled illumination for the night",
  maxMoonlightIntensity = "Max value of modeled illumination for the night",
  meanMoonPhase = "Mean value of moon phase (% of moon illuminated)",
  minMoonPhase = "Min value of moon phase (% of moon illuminated)",
  maxMoonPhase = "Max value of moon phase (% of moon illuminated)"
)

metadata_dt <- as.data.table(t(metadata), keep.rownames = "Variable")

# Combine the metadata and moon.stat data frames
combined_data <- rbindlist(list(metadata_dt, moon.df), use.names = FALSE, fill = TRUE)

# check combined_data

print(combined_data)


summary(moon.df)









# outputs -----------------------------------------------------------------
fwrite(combined_data, file = "data_for_analysis/moon_pred/moon.stat.csv", row.names = FALSE)

fwrite(moon.int, file = "data_for_analysis/moon_pred/moon.int.csv", row.names = FALSE)

fwrite(moon.int.night, file = "data_for_analysis/moon_pred/moon.int.night.csv", row.names = FALSE)

# 

# Create a README  to explain the moon.stat and moon.int data sets.
readme_content <- "Carlos Linares 2/11/2025
This directory contains the files produce with the moonlit package in collaboration with Kyle Shannon.
-moon.int.csv file contains the moon illumination for each date and time of the bat data set.
-moon.stat.csv file contains the sunrise and sunset of and moon phase and intensity valuesfor each night of the bat data set.
-moon.int.night.csv file contains the moon illumination for each date and time of the bat data set but only at night.
-sidalo file contains the coordinates and date time for each bat call. in bat_combined but I filtered noise out

#  moon.stat
# sunset = Time of sunset,
#   sunrise = Time of sunrise,
#   meanMoonlightIntensity = Mean value of modeled illumination for the night,
#   minMoonlightIntensity = Min value of modeled illumination for the night,
#   maxMoonlightIntensity = Max value of modeled illumination for the night,
#   meanMoonPhase = Mean value of moon phase (% of moon illuminated),
#   minMoonPhase = Min value of moon phase (% of moon illuminated),
#   maxMoonPhase = Max value of moon phase (% of moon illuminated)

moon. int
night = Logical, TRUE when sun is below the horizon,
  sunAltDegrees = Solar altitude in degrees,
  moonlightModel = Predicted moonlight illumination, relative to an 'average' full moon,
  twilightModel = Predicted twilight illumination in lux,
  illumination = Combined moon and twilight intensity, in lux,
  moonPhase = Lunar phase as a numerical value (% of moon face illuminated)

"




# Write the README content to a file
writeLines(readme_content, "data_for_analysis/moon_pred/README.txt")





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




# moonlit graphs ----------------------------------------------------------


head(moon.int)

# plot illumination vs jday

a<-ggplot(moon.int, aes(x = yday(date), y = illumination)) +
  geom_point() +
  labs(title = "Moonlight Intensity vs. Julian Day",
       x = "Julian Day",
       y = "Moonlight Intensity (lux)") +
  theme_minimal()
a

# seems like we have very large values at some point... makes me unsure we can trust the data.
# show me the 10 larges values for illumination. The object below showed me there were observations in the mideld of the day and I filter them out. created object moon.int.night
maxillum<-moon.int %>% arrange(desc(illumination)) %>% head(10)

# a graph with lots of data but not much information. 
a1 <- ggplot(moon.int.night, aes(yday(date), illumination)) +
  geom_point(aes(color = as.factor(year))) +  # Use the new year column for coloring
  labs(x = "Day of Year", y = "Illumination", color = "Year") +
  theme_minimal()

a1

# The graph below shows we have several observations below the horizon. Thus we need to make those zeros because the moon illumination is not relevant when the moon is below the horizon. However we don't need to do this anymore because the moonlit package does this automatically. Illumination is the sum of moon and twilight illumination.Thus no need to modify. 

b<-ggplot(moon.int.night, aes(x=moonAltDegrees, y=illumination))+
  geom_point()+
  geom_vline(xintercept = 0)
b

# now we plot just moon illumination.

c<-ggplot(moon.int, aes(x = yday(date), y = moonlightModel)) +
  geom_point() +
  labs(title = "Moonlight Intensity vs. Julian Day",
       x = "Julian Day",
       y = "Moonlight Intensity (lux)") +
  theme_minimal()
c

# saving plots. 

figures_dir <- "figures/moon_pred"


plots <- list(a1, b, c)
filenames <- c("a1.tiff", "b.tiff", "c.tiff")

# Save all plots
lapply(seq_along(plots), function(i) {
  ggsave(filename = file.path(figures_dir, filenames[i]), plot = plots[[i]], device = "tiff", width = 15, height = 10, units = "cm")
})


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

# Create directory if it doesn't exist
dir.create("data_for_analysis/moon_pred", recursive = TRUE, showWarnings = FALSE)

# saveRDS(sidalo, file = "data_for_analysis/moon_pred/sidalo.csv") this saves it in a non human readable format.

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