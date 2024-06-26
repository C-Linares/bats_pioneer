
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



# load 2021 2019 ----------------------------------------------------------


bat2021_raw<- read.csv(file = 'data_for_analysis/data2021_db_only/bats2021_v4.csv')
bat2019_raw<- read.csv(file = 'data2019_db_only/bats2019.csv', header = T, check.names = T)

# we check there's no NAs

colSums(is.na(bat2021_raw)) # there are tons but that's before cleanup 
colSums(is.na(bat2019_raw))

#clean up data

keep<- c(".id","Path","Filename","HiF","LoF", "SppAccp", "Prob","X.Spp", "X.Prob","X1st")

bat2019_v1 <- bat2019_raw %>% select(all_of(keep)) # removes unnecessary columns 

bat2019_v1[bat2019_v1==""]<- NA # makes the empty spaces NAs

bat2019_v2 <- bat2019_v1 %>%  filter(!is.na(SppAccp)) # removes NAs 

bat2021_v2 <- bat2021_raw %>% select(all_of(keep)) # this stopped working if you load the raster package
# bat2021_v2 <- bat2021_raw[, keep, drop = FALSE]

bat2021_v2[bat2021_v2==""]<- NA # makes the empty spaces NAs

# bat2021_v2 <- bat2021_v2 %>%  filter(!is.na(SppAccp)) #leaves us with just bat observations


# time extraction 2021 2019 ---------------------------------------------------------

bat2021_v2 <- bat2021_v2 %>%
  mutate(date_time = ymd_hms(str_extract(Filename, "\\d{8}_\\d{6}"), tz = "America/Denver")) 

bat2021_v2 <- bat2021_v2 %>%
  mutate(
    date_time = tryCatch( # this puts an NA when the date failed to parce. 
      ymd_hms(str_extract(Filename, "\\d{8}_\\d{6}"), tz = "America/Denver"),
      error = function(e) NA  # Set to NA if an error occurs
    )
  )

rows_with_na <- subset(bat2021_v2, is.na(date_time))

bat2021_v2 <- bat2021_v2[complete.cases(bat2021_v2$date_time), ] # remove NAs from the date_time col

bat2021_v2$hrs<- hour(bat2021_v2$date_time) #extract the hour. 

#We want to extract the data and time from the Filename and Path columns. There must be a more efficient way fo doing this but I don't know how to. 

bat2019_v2$Filename<-as.character(bat2019_v2$Filename)
bat2019_v2$Path<- as.character(bat2019_v2$Path)


bat2019_v2$date<-
  unlist(lapply(strsplit(bat2019_v2$Filename, "_"), function(x)
    x[5]))

bat2019_v2$hms_sp<-
  unlist(lapply(strsplit(bat2019_v2$Filename, "_"), function(x)
    # this splits the string by _ and give us the time-species[6]
    x[6]))


bat2019_v2$hms<-
  unlist(lapply(strsplit(bat2019_v2$hms_sp,"\\."), function(x)
    x[1]))

# site 2021 ---------------------------------------------------------------

bat2021_v2$site<-str_extract(bat2021_v2$Filename, "^[A-Za-z]{3,4}\\d{2}")

# site_info <- regmatches(bat2021_v2$Filename, regexpr("[A-Za-z]\\d{2}", bat2021_v2$Filename)) another way to do it


bat2019_v2$site<-
  unlist(lapply(strsplit(bat2019_v2$Path, "\\\\"), function(x) # how to do this
    x[5]))
  
unique(bat2021_v2$site)

# here we fix some problems with the sites misspellings


bat2021_v2$site = ifelse(bat2021_v2$site %in% c("LON01-"),"LON01", bat2021_v2$site)

bat2021_v2$site = ifelse(bat2021_v2$site %in% c("LON05-"),"LON05", bat2021_v2$site)

bat2021_v2$site = ifelse(bat2021_v2$site %in% c("VIS02"),"VIZ02", bat2021_v2$site)

bat2021_v2$site = ifelse(bat2021_v2$site %in% c("VIS01"),"VIZ01", bat2021_v2$site)

bat2021_v2$site = ifelse(bat2021_v2$site %in% c("VIS04"),"VIZ04", bat2021_v2$site)

bat2021_v2$site = ifelse(bat2021_v2$site %in% c("VIS03"),"VIZ03", bat2021_v2$site)


unique(bat2021_v2$site) # we are missing long02 

# checking NAs 2021 2029 --------------------------------------------------

na_rows <- bat2021_v2[is.na(bat2021_v2$hrs) | is.na(bat2021_v2$X1st), ] # shows the NA rows

nas2019<- bat2019_v2$hms[!complete.cases(bat2019_v2)] # other way of doing the same for the above file. 

# remover extra columns not needed. 

bat2021_v2<-bat2021_v2 %>% select(-Path, -Filename)
bat2019_v2<- bat2019_v2 %>% select(-Path, -Filename, -hms_sp)


# Noche  ------------------------------------------------------------------

bat2021_v2$noche <-
  if_else(bat2021_v2$hrs < 9, # if it is less than 9 put the date of the previous day
          true =  (date(bat2021_v2$date_time) - ddays(1)),
          false = date(bat2021_v2$date_time))


bat2019_v2$noche <-
  if_else(bat2019_v2$hr < 9,
          true =  (bat2019_v2$date - ddays(1)),
          false = bat2019_v2$date)

### EFFORT ---------------------------------

#we need to calculate the number of hrs each site was working. 
# we will do this by getting difference between the first recording by site and the last one for that day. that will give us an estimate for the number of hrs each recording was working


effort <- bat2021_v2 %>%
  group_by(site, noche) %>%
  summarise(stard = min(date_time), endd = max(date_time)) %>%
  mutate(eff.hrs = time_length(endd - stard, unit = "hours")) %>% 
  mutate(wk=week(noche))




effort_days <- bat2021_v2 %>%
  group_by(site) %>%
  summarise(
    stard = min(date_time),
    endd = max(date_time),
    eff.days = as.numeric(difftime(max(date_time), min(date_time), units = "days"))
  )

write.csv(effort,file = 'data_for_analysis/effort2021.csv',)


table(bat2021_v2$site,bat2021_v2$SppAccp) #how many rows pe

# drop unnecessary variables to join with bat2021_v2

effort <- effort %>% select("site","noche", "n.hrs") 

bat2021_v2<- left_join(x = bat2021_v2, y = effort, by= c("noche","site"))
bat2021_v2<- left_join(x = bat2021_v2, y = effort_days, by= c("site"))

bat2021_v2<-bat2021_v2 %>% select(-stard, -endd)

# suncal ------------------------------------------------------------------

###-------------------
# now we need times since sunset. I calculated these by the camp site coordinates and not each site coordinates
# we calculate sunset = sun disappears below the horizon. 
# 
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
bat2021_v2$dlt.sunset.hrs<-time_length(bat2021_v2$date_time- bat2021_v2$sunset, "hours") 
# check for negative times in the dlt.sunset col
check1<-subset(bat2021_v2,subset = bat2021_v2$dlt.sunset<0) #491 bat detection occurred before sunset. now 234 with new data

# joining sppacc and ~spp -------------------------------------------------

# to increase the data available, we will mixed the spp accp with the 

bat2021_v2 <- mutate(bat2021_v2, Sp.plus = coalesce(SppAccp, X.Spp))



# lit and dark sites -----------------

stnames<-unique(bat2021_v2$site)
litcond<-c("lit", "dark")
litsites<-c("IRON01","IRON03","IRON05","LON01","LON03")


bat2021_v2$treatmt<-ifelse(bat2021_v2$site %in% litsites , "lit", "dark") # this makes a treatment variable.

bat2021_v2$trmt_bin<- ifelse(bat2021_v2$treatmt== "lit", 1, 0) #this makes that variable binary

# bat2021_v2 <- bat2021_v2 %>%  filter(!is.na(SppAccp)) #leaves us with just bat observations

# we write the data until we get the data for temperature

write.csv(bat2021_v2,file = 'data_for_analysis/bat2021_v3.csv') # updated the bat2021 data base. 
write.csv(bat2019_v2,file = "data_for_analysis/bat2019_v2.csv")



# sites lat long ----------------------------------------------------------


# here we add the lat long for each site to the bat2021_v2
 sts<-read.csv('data_for_analysis/sites_coordinates.csv')
# bat2021_v2<-read.csv('data_for_analysis/bat2021_v2.csv')

bat2021_v2<-subset(bat2021_v2, select = -c(lat,lon))#drop the bat2021_v2 lat and lon to update
bat2021_v2<-left_join(bat2021_v2,sts) # adds the correct lat and long coordinates for the sites. 


write.csv(bat2021_v2,file = 'data_for_analysis/bat2021_v3.csv',) # updated the bat2021 data base. 


#---- comparing the 2021 data produce with the v3 and v4 of the 
#
#Issue: it seems like the original data base has more data than the new one. I am not getting all the sonobat files 


a<-fread('data_for_analysis/bats2021_update.csv',drop = "V1")
b<-fread('data_for_analysis/data2021_db_only/bats2021.csv')

setkey(a,"Path")
setkey(b,"Path")

# Find rows present in only one of the data frames
diff_rows <- a[!b]
diff_rows <- rbind(diff_rows, b[!a], fill=T)  # Combine differences from both sides

library(DT)
datatable(diff_rows)



#  2022 Sonobat data  -----------------------------------------------------

# load data.

bat2022_raw<-read.csv('data_for_analysis/data2022_db_only/bats2022.csv')
 
keep<- c("Path","Filename","HiF","LoF", "SppAccp", "Prob", "calls.sec", "X1st", "X2nd")# cols to keep

bat2022_v1 <- bat2022_raw %>% select(all_of(keep)) # removes unnecessary columns 

colSums(is.na(bat2022_raw))# the raw file apparently has no NA's

bat2022_v1[bat2022_v1==""]<- NA # makes the empty spaces NA's

bat2022_v2 <- bat2022_v1 %>%  filter(!is.na(X1st )) # removes NAs


# sites 2022 -------------------------------------------------------------------

sites<-str_extract(bat2022_raw$Filename, "^[A-Za-z]{3,4}\\d{2}")


# time extraction ---------------------------------------------------------


bat2022_v2 <- bat2022_v2 %>%
  mutate(date_time = ymd_hms(str_extract(Filename, "\\d{8}_\\d{6}"), tz = "America/Denver"))

# bat2022_v2$t <- str_extract(string = bat2022_v2$Filename,pattern = "\\d{8}_\\d{6}") # gets the date and time string
# bat2022_v2$date_time<- ymd_hms(bat2022_v2$t,tz ="America/Denver") # creates a data time

bat2022_v2<-bat2022_v2[,-10] # remove the t col 

# make hour 

bat2022_v2$hr<-hour(bat2022_v2$date_time)

# create a data for the night. 

bat2022_v2$noche <-
  if_else(bat2022_v2$hr < 9, # if it is less than 9 put the date of the previous day
          true =  (date(bat2022_v2$date_time) - ddays(1)),
          false = date(bat2022_v2$date_time))























# trash -------------------------------------------------------------------

# here we remove the NAs for the date time. 
bat2021_v2 <- bat2021_v2[complete.cases(bat2021_v2$date_time), ]

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




# kpro clean up -----------------------------------------------------------


kpro2021_raw<- read.csv(file = 'data_for_analysis/2021_kpro_raw/bats2021_kpro_v1.csv',header = T)

keep<- c(".id","INDIR","OUTDIR","FOLDER","IN.FILE","DURATION","DATE","TIME","HOUR","AUTO.ID.")

kpro2021_raw <- kpro2021_raw %>% select(all_of(keep))



# site 

kpro2021_raw$site<-str_extract(kpro2021_raw$OUTDIR, "[A-Za-z]{3,4}\\d{2}")

unique(kpro2021_raw$site) # we have all sites

# date week

kpro2021_raw$DATE<-lubridate::ymd(kpro2021_raw$DATE)

# noche/night

kpro2021_raw$noche <-
  if_else(kpro2021_raw$HOUR < 9, # if it is less than 9 put the date of the previous day
          true =  (date(kpro2021_raw$DATE) - ddays(1)),
          false = date(kpro2021_raw$DATE))
#week
kpro2021_raw$wk<- week(kpro2021_raw$noche)


# date time col. 
datetime<-paste(kpro2021_raw$DATE, kpro2021_raw$TIME)#merge date and time
datetime.parse<-lubridate::ymd_hms(datetime) # parse as date time
kpro2021_raw$date_time<-datetime.parse # add to data. 


# 
litsites<-c("iron01","iron03","iron05","long01","long03")


kpro2021_raw$treatmt<-ifelse(kpro2021_raw$site %in% litsites , "lit", "dark") # this makes a treatment variable.

kpro2021_raw$trmt_bin<- ifelse(kpro2021_raw$treatmt== "lit", 1, 0)


kpro2021_raw$jday<-lubridate::yday(kpro2021_raw$noche) # julian day


# write the data. 

write.csv(kpro2021_raw, file='data_for_analysis/kpro2021_v1.csv',  row.names = F)



