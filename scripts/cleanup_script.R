
# This script is for cleaning the data bases for the bat acoustic data Years processed: 2019(not pioneer data), 2021, 2022(pending),2023(pending) 

# aims:
 
# Compile all the output files and contruct a single data base. 

# create a table with the meaning of each of the columns outputed from SonoBat

# create a date column but the date should be from the same 6pm to 5 am so it correspond to one night. 

# create a site column


#libraries 

library(stringr)
library(lubridate)
library(tidyverse)
library(data.table)
library(suncalc) # calculate sunset and sunrise.


# load the data fieles

bat2021_raw<- read.csv(file = 'data2021_db_only/bats2021.csv', header = T)
bat2019_raw<- read.csv(file = 'data2019_db_only/bats2019.csv', header = T, check.names = T)

# we check there's no NAs

colSums(is.na(bat2021_raw)) # there are tons but that's before cleanup 
colSums(is.na(bat2019_raw))

#clean up data

keep<- c("Path","Filename","HiF","LoF", "SppAccp", "Prob", "calls.sec", "X1st")

bat2019_v1 <- bat2019_raw %>% select(all_of(keep)) # removes uncecessary columns 

bat2019_v1[bat2019_v1==""]<- NA # makes the empty spaces NAs

bat2019_v2 <- bat2019_v1 %>%  filter(!is.na(SppAccp)) # removes NAs 

bat2021_v2 <- bat2021_raw %>% select(all_of(keep)) # removes uncecessary columns 

bat2021_v2[bat2021_v2==""]<- NA # makes the empty spaces NAs

# bat2021_v2 <- bat2021_v2 %>%  filter(!is.na(SppAccp)) #leaves us with just bat observations

#We want to extract the data and time from the Filename and Path columns. There must be a more efficient way fo doing this but I don't know how to. 

bat2021_v2$Filename<- as.character(bat2021_v2$Filename) # we make them character. 
bat2021_v2$Path<- as.character(bat2021_v2$Path)         # we make them character.

bat2019_v2$Filename<-as.character(bat2019_v2$Filename)
bat2019_v2$Path<- as.character(bat2019_v2$Path)

bat2021_v2$date <-
  unlist(lapply(strsplit(bat2021_v2$Filename, "_"), function(x)
    # this splits the string by _ and give us the date[5]
    x[5]))

bat2019_v2$date<-
  unlist(lapply(strsplit(bat2019_v2$Filename, "_"), function(x)
    x[5]))

bat2021_v2$hms_sp <-
  unlist(lapply(strsplit(bat2021_v2$Filename, "_"), function(x)
    # this splits the string by _ and give us the time-species[6]
    x[6]))

bat2019_v2$hms_sp<-
  unlist(lapply(strsplit(bat2019_v2$Filename, "_"), function(x)
    # this splits the string by _ and give us the time-species[6]
    x[6]))


bat2021_v2$hms <-
  unlist(lapply(strsplit(bat2021_v2$hms_sp, "-"), function(x)
    # this splits the string by _ and give us the time[6]
    x[1]))

bat2019_v2$hms<-
  unlist(lapply(strsplit(bat2019_v2$hms_sp,"\\."), function(x)
    x[1]))

# bat2021_v2$hms2 <-
#   unlist(lapply(strsplit(bat2021_v2$hms, "1"), function(x)
#     # this splits the string by _ and give us the time[6]
#     x[1]))


bat2021_v2$site <-
  unlist(lapply(strsplit(bat2021_v2$Filename, "_"), function(x) # we are missing VIZ02 and LON02
    x[1]))

bat2019_v2$site<-
  unlist(lapply(strsplit(bat2019_v2$Path, "\\\\"), function(x) # how to do this
    x[5]))
  

# here we fix some problems with the sites misspellings

bat2021_v2$site = ifelse(bat2021_v2$site %in% c("LON01-"),"LON01", bat2021_v2$site)

bat2021_v2$site = ifelse(bat2021_v2$site %in% c("LON05-"),"LON05", bat2021_v2$site)

bat2021_v2$site = ifelse(bat2021_v2$site %in% c("VIS02"),"VIZ02", bat2021_v2$site)

bat2021_v2$site = ifelse(bat2021_v2$site %in% c("VIS01"),"VIZ01", bat2021_v2$site)

bat2021_v2$site = ifelse(bat2021_v2$site %in% c("VIS04"),"VIZ04", bat2021_v2$site)

unique(bat2021_v2$site) # we are missing viz03 and long02



#checking Nas

nas<- bat2021_v2$hms[!complete.cases(bat2021_v2)] # this should be empty.

nas2019<- bat2019_v2$hms[!complete.cases(bat2019_v2)] # other way of doing the same for the above file. 

# remover extra columns not needed. 

bat2021_v2<-bat2021_v2 %>% select(-Path, -Filename, -hms_sp)
bat2019_v2<- bat2019_v2 %>% select(-Path, -Filename, -hms_sp)

# we need a column that assigns the the data of of the night. For example if we have 6pm-4am recording then the night should be the date of the start of the recording. 

# First date as date with lubridate 

bat2021_v2$date<-ymd(bat2021_v2$date)

bat2019_v2$date<-ymd(bat2019_v2$date)

# hms as numeric or actual hour? I am struggling with what is the best way. It seems like I need to create a new variable that is date and time together in a single column. 

bat2021_v2$datetime<- paste(bat2021_v2$date, bat2021_v2$hms) %>% 
  ymd_hms() 

bat2019_v2$datetime<- paste(bat2019_v2$date, bat2019_v2$hms) %>% 
  ymd_hms()

# bat2021_v2$hms2<- as.numeric(bat2021_v2$hms) # this does not work 
# bat2021_v2$hms2<- hms::as_hms(bat2021_v2$hms2) # this does not work

which(is.na(bat2021_v2)) # shows the location of the NAs 

# here we make an hour column.

bat2021_v2$hr<-lubridate::hour(bat2021_v2$datetime)
bat2021_v2$min<-lubridate::minute(bat2021_v2$datetime)

bat2019_v2$hr<-lubridate::hour(bat2019_v2$datetime)
bat2019_v2$min<-lubridate::minute(bat2019_v2$datetime)
# if hms<090000 then date-1

bat2021_v2$noche <-
  if_else(bat2021_v2$hr < 9, # if it is less than 9 put the date of the previous day
          true =  (bat2021_v2$date - ddays(1)),
          false = bat2021_v2$date) # this created the "noche" column that is the date for the night of the start of the recording.

bat2019_v2$noche <-
  if_else(bat2019_v2$hr < 9,
          true =  (bat2019_v2$date - ddays(1)),
          false = bat2019_v2$date)

### EFFORT ---------------------------------

#we need to calculate the number of hrs each site was working. 
# we will do this by getting difference between the first recording by site and the last one for that day. that will give us an estimate for the number of hrs each recording was working

str(bat2021_v2)

effort<- bat2021_v2 %>% 
  group_by(site, noche) %>% 
  summarise(stard= min(datetime), endd= max(datetime)) %>% 
  mutate(n.hrs= time_length(endd-stard, unit="hours"))

# drop unnecessary variables to join with bat2021_v2

effort <- effort %>% select("site","noche", "n.hrs") 

t<- left_join(x = bat2021_v2, y = effort, by= c("noche","site"))

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
moonsun2019 <- left_join(x = moon2019, y = sun2019, "date")
colnames(moonsun2019)[1]<-"noche"

#join with the rest of data 
bat2021_v2 <- left_join(x = bat2021_v2, y = moonsun, "noche")
bat2019_v2<- left_join(x=bat2019_v2, y=moonsun2019, "noche")

#now we calculate the time since sunset with "sunset" from the package suncal

#first we transform the sunset time into UTC 
toutc<- c("sunrise","sunriseEnd","sunset","sunsetStart")


bat2021_v2$sunset<-ymd_hms(bat2021_v2$sunset,tz="UTC")
bat2021_v2$dlt.sunset.hrs<-time_length(bat2021_v2$datetime- bat2021_v2$sunset, "hours") 
# check for negative times in the dlt.sunset col
check1<-subset(bat2021_v2,subset = bat2021_v2$dlt.sunset<0) #491 bat detection occurred before sunset. 


# lit and dark sites -----------------

stnames<-unique(bat2021_v2$site)
litcond<-c("lit", "dark")
litsites<-c("IRON01","IRON03","IRON05","LON01","LON03")


bat2021_v2$treatmt<-ifelse(bat2021_v2$site %in% litsites , "lit", "dark") # this makes a treatment variable.

bat2021_v2$trmt_bin<- ifelse(bat2021_v2$treatmt== "lit", 1, 0) #this makes that variable binary

# bat2021_v2 <- bat2021_v2 %>%  filter(!is.na(SppAccp)) #leaves us with just bat observations

# we write the data until we get the data for temperature

write.csv(bat2021_v2,file = 'data_analysis/bat2021_v2.csv') 
write.csv(bat2019_v2,file = "data_analysis/bat2019_v2.csv")


# here we add the lat long for each site to the bat2021_v2
sts<-read.csv('data_analysis/sites_coordinates.csv')
bat2021_v2<-read.csv('data_analysis/bat2021_v2.csv')

bat2021_v2<-subset(bat2021_v2, select = -c(lat,lon))#drop the bat2021_v2 lat and lon to update
bat2021_v2<-left_join(bat2021_v2,sts) # adds the correct lat and long coordinates for the sites. 

#------elevation -------
#added elevation to points worked 10/15/2023 but apparently the 

library(raster)
library(rgdal)
sts<-read.csv('data_analysis/sites_coordinates.csv')
dem<-raster('data_analysis/elev/USGS_13_n44w114_20130911.tif')
coordinates<-data.frame(lon=sts$lon, lat=sts$lat)
coordinates_sp <- SpatialPoints(coordinates, proj4string = CRS(proj4string(dem)))
sts$elevation <- extract(dem, coordinates_sp)

# we add elevation to the bat data.

# bat2021_v2<-read.csv('data_analysis/bat2021_v2.csv')
bat2021_v2<-left_join(bat2021_v2, sts, by="site")

write.csv(bat2021_v2,file = 'data_analysis/bat2021_v2.csv') # we update the data base

# now we need a temperature... where do I get the temperature? Apparently there is the PRISM data that could be useful. 



# percentage of riparian area?
# We will use NDVI data.

ndvi_2021<-raster('data_analysis/NDVI/MYD13Q1.A2021185.h09v04.061.2021202231848.hdf',crs=)
coordinates_sp <- SpatialPoints(coordinates, proj4string = CRS(proj4string(ndvi_2021)))

buffer_distance <-100

# Create a buffer around the point
buffer <- raster::buffer(coordinates_sp, width = buffer_distance)

# Extract the NDVI values for the buffer
ndvi_values <- extract(ndvi_2021, buffer)

dvi_stats <- lapply(ndvi_values, function(x) {
  if (length(x) > 0) {
    # Calculate mean NDVI within the buffer
    mean_ndvi <- mean(x, na.rm = TRUE)
    return(mean_ndvi)
  } else {
    # If no NDVI values were found within the buffer, return NA
    return(NA)
  }
})
