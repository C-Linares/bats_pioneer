
# R Script: **sunset and delta sunset** (color: blue)

#---------------------------------------------------------#

# Author: **Carlos**
# Date: **4/8/2024**
# Purpose: **calculate the sun and moon coovariates**
# 


library(suncalc)

#load bat data 

kpro_2021_bat<-read.csv(file = 'data_for_analysis/kpro2021_v1.csv') # raw kpro data
kpro_2021_bat$noche<-as.Date(kpro_2021_bat$noche)

datetime<-paste(kpro_2021_bat$DATE, kpro_2021_bat$TIME)#merge date and time
datetime.parse<-lubridate::ymd_hms(datetime) # parse as date time
kpro_2021_bat$date_time<-datetime.parse # add to data. 
###-------------------
# now we need times since sunset. I calculated these by the camp site coordinates and not each site coordinates
# we calculate sunset = sun disappears below the horizon.


 
sun<-getSunlightTimes(date = seq.Date(from = as.Date("2021-06-22"),to = as.Date("2021-08-24"), by=1),
                      keep = c("sunrise", "sunriseEnd", "sunset", "sunsetStart"),
                      lat = 43.5479918,
                      lon = -113.7362989,
                      tz= "America/Denver")

sun2019 <-
  getSunlightTimes(
    date = seq.Date(
      from = as.Date("2019-05-08"),
      to = as.Date("2019-07-24"),
      by = 1
    ),
    keep = c("sunrise", "sunriseEnd", "sunset", "sunsetStart"),
    lat = 43.5479918,
    lon = -113.7362989,
    tz = "America/Denver"
  )

moon<- getMoonIllumination(date = seq.Date(from = as.Date("2021-06-22"),to = as.Date("2021-08-24"),
                                           by=1),
                           keep = c("phase", "fraction"))

moon2019 <-
  getMoonIllumination(
    date = seq.Date(
      from = as.Date("2019-05-08"),
      to = as.Date("2019-07-24"),
      by = 1
    ),
    keep = c("phase", "fraction")
  )

# now we join these with the bat data. 

moonsun <- left_join(x = moon, y = sun, "date")
colnames(moonsun)[1]<-"noche" # we need to make it noche so it can left join by="noche"
moonsun$noche<-as.Date(moonsun$noche)
moonsun2019 <- left_join(x = moon2019, y = sun2019, "date")
colnames(moonsun2019)[1]<-"noche"

#join with the rest of data 
kpro_2021_bat <- left_join(x = kpro_2021_bat, y = moonsun, "noche")

#now we calculate the time since sunset with "sunset" from the package suncal


kpro_2021_bat$sunset<-ymd_hms(kpro_2021_bat$sunset,tz="UTC")
kpro_2021_bat$dlt.sunset.hrs<-time_length(kpro_2021_bat$date_time- kpro_2021_bat$sunset, "hours") 
# check for negative times in the dlt.sunset col
check1<-subset(kpro_2021_bat,subset = kpro_2021_bat$dlt.sunset<0) #491 bat detection occurred before sunset. now 234 with new data

