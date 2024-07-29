# 
# Script name: prep_for_glmm
# 
# Purpose of script: Combine the 2021, 2022, 2023, data sets and get them ready to go throuhg glmm_v1
# 
# Author: Carlos Linares
# 
# Date Created: 07/29/2024
# 
# Email: carlosgarcialina@u.boisestate.edu
# 
# ---------------------------
#   
#   Notes: we also should check we have the same temporal window for all variables.  
# ---------------------------


# inputs ------------------------------------------------------------------
-data_for_analysis/2021_kpro_raw/bats2021_kpro_v1.csv
-data_for_analysis/2022_kpro_raw/bat2022_kpr.csv
-data_for_analysis/2023_kpro_raw/bat2023_kpr.csv

# outputs ----------------------

# this should be a database ready to analyze with the glmm_v1 script. 

# libraries

library(tidyverse)
library(beepr)


# kpro_data --------------------------------------------------------------

# load

kpro_2021_bat<- read.csv(file = 'data_for_analysis/2021_kpro_raw/bats2021_kpro_v1.csv',header = T)
kpro_2022_bat<-read.csv(file = 'data_for_analysis/2022_kpro_raw/bat2022_kpr.csv',header = T)
kpro_2023_bat<-read.csv(file = 'data_for_analysis/2023_kpro_raw/bat2023_kpr.csv',header = T)

bat_combined <- bind_rows(kpro_2021_bat, kpro_2022_bat, kpro_2023_bat)
str(bat_combined)

keep<- c(".id","INDIR","OUTDIR","FOLDER","IN.FILE","DURATION","DATE","TIME","HOUR","AUTO.ID.") # cols to keep

bat_combined <- bat_combined %>% select(all_of(keep))


# site 

bat_combined$site<-str_extract(bat_combined$OUTDIR, "[A-Za-z]{3,4}\\d{2}")

unique(kpro2021_raw$site) # site labells

# date week

bat_combined$DATE<-lubridate::ymd(bat_combined$DATE)
sum(is.na(bat_combined$DATE)) # check for NAs. 

# noche/night

bat_combined$noche <-
  if_else(bat_combined$HOUR < 9, # if it is less than 9 put the date of the previous day
          true =  (date(bat_combined$DATE) - ddays(1)),
          false = date(bat_combined$DATE))



# date time col. 
datetime<-paste(bat_combined$DATE, bat_combined$TIME)#merge date and time
datetime.parse<-lubridate::ymd_hms(datetime) # parse as date time
bat_combined$date_time<-datetime.parse # add to data. 
sum(is.na(bat_combined$date_time)) # check for NAs. 

#year

bat_combined$yr<-year(bat_combined$DATE)
unique(bat_combined$yr)

litsites<-c("iron01","iron03","iron05","long01","long03")


bat_combined$treatmt<-ifelse(bat_combined$site %in% litsites , "lit", "dark") # this makes a treatment variable.

bat_combined$trmt_bin<- ifelse(bat_combined$treatmt== "lit", 1, 0)


bat_combined$jday<-lubridate::yday(bat_combined$noche) # julian day


summary(bat_combined)




# add effort

effort_days <- bat_combined %>%
  group_by(site) %>%
  summarise(
    stard = min(noche),
    endd = max(noche),
    eff.days = as.numeric(difftime(max(noche), min(noche), units = "days"))
  )

effort_hrs <- bat_combined %>%
  group_by(site, noche, jday) %>%
  summarise(stard = min(date_time), endd = max(date_time)) %>%
  mutate(eff.hrs = time_length(endd - stard, unit = "hours"))

# merge effort with bat combined 

bat_combined<- left_join(bat_combined, effort_hrs, by=c("site", "jday"))









# write the data. 

write.csv(kpro2021_raw, file='data_for_analysis/kpro2021_v1.csv',  row.names = F)