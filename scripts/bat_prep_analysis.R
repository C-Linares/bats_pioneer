
# Script: bat_pot_analysis.R
# The purpose of this script is to prepare the data before modeling. 
# 
# Carlos Linares 2023
# 



# 
# #
# libraries 
library(tidyverse)
library(lubridate)
library(magrittr)
library(elevatr)

# Data load-------------####

bat2021<-read.csv('data_for_analysis/bat2021_v3.csv',header = T,stringsAsFactors = F) # we load the data we created in the clean up script. This data has several variables added to a clean up raw data is the product of clean_up script
table(is.na(bat2021$date_time))

# make dates-------------
bat2021$date_time <- lubridate::ymd_hms(bat2021$date_time,tz = "America/Denver") # makes dates as dates
failed_rows <- is.na(bat2021$date_time)
failed_dates <- bat2021[failed_rows, ]

bat2021$noche <-
  ymd(bat2021$noche)# makes noche as date. 

# noche is the date for all the calls that come from the same sampling night.

bat2021$jday <- yday(bat2021$noche)
bat2021$wk <-
  week(bat2021$noche)# we need to calculate a week column.
bat2021$yr <-
  year(bat2021$noche) #year added for when we will have multiple years.

summary(bat2021) # there are no NAs in this database


# Sp.plus cleanup  --------------------------------------------------------
# lets talk to jesse about this one.

table(bat2021$Sp.plus)
table(unique(bat2021$Sp.plus))
t <- mutate(bat2021, Sp.plus = coalesce(SppAccp, X1st))



#sites
unique(bat2021$site)# tell us what sites we have

# Matrix building  ------------------

bat1 <- bat2021 %>%
  group_by(site, treatmt) %>%
  count(SppAccp, jday, hrs,drop=F) %>%
  ungroup()
summary(bat1)

#we filter by just one sp.


bat_js <- bat2021 %>%
  group_by(site, SppAccp) %>% # I don't include year because it is a single year
  count(wk, .drop = FALSE) %>%  # we might have to include the argument .drop=false to count the NAs and the zeros
  pivot_wider(names_from = wk, values_from = n) %>%
  ungroup()
summary(bat_js)


bat_js<-bat_js[,sort(colnames(bat_js))] # sort the cols
bat_js<-replace(bat_js, is.na(bat_js),0)

lano_js <- bat_js[bat_js$SppAccp == "Lano",] # filter data to just Lano.
lano_js <-lano_js %>%  select(!c("34", "site", "SppAccp")) # remove site, sp, and the last week because it has  data for just one site.

lano_js <- replace(lano_js, is.na(lano_js), 0) # NAs to zeros 
lano_js <- lano_js[,sort(colnames(lano_js))] # sort the cols

mylu_js <- bat_js[bat_js$SppAccp == "Mylu", ]
mylu_js <- mylu_js[!is.na(mylu_js$SppAccp), ]
mylu_js <- replace(mylu_js, is.na(mylu_js), 0) # NAs to zeros 
mylu_js <-mylu_js %>%  select(!c("34", "site", "SppAccp"))

# we filter by just one sp.
# what species has the more calls. Lano Mylu
table(bat2021$SppAccp)

ggplot(bat2021,aes(x=SppAccp, fill=treatmt))+
  geom_bar()
  

# funtion to filter bat data ----------------------------------------------

filter_and_clean_data <- function(data, species) {
  filtered_data <- data %>%
    filter(SppAccp == species) %>%
    drop_na(SppAccp) %>%
    mutate(across(where(is.numeric), ~if_else(is.na(.), 0, .))) %>%
    select(-c("34", "site", "SppAccp"))
  
  return(filtered_data)
}

myvo_js<-filter_and_clean_data(bat_js, "Myvo")
#rename week columns as numbers just numbers 

write.csv(bat_js,file = 'data_analysis/bat_prep_analysis/bat_js.csv',
          row.names = F) 

#write single species df

write.csv(lano_js,file = 'data_for_analysis/bat_pop_analysis/lano_js.csv', row.names = F) 
write.csv(mylu_js,file = 'data_for_analysis/bat_pop_analysis/mylu_js.csv', row.names = F) 


# site level cov. -------------
# pending veg, temp and rain.

s.l.c <- bat2021 %>% dplyr::distinct(site,elevation,trmt_bin)
s.l.c<-s.l.c[,-1] # remove sites column

write.csv(s.l.c,file = 'data_analysis/bat_prep_analysis/slc.csv',
          row.names = F) 

# obs cov ----------------


moon_pred <-
  read.csv('data_for_analysis/moon_pred.csv') #loads the data from the moon_pred script

moon_pred <- moon_pred %>%  select(!c("X.1", "X"))# removes random row col names added when writing the file.


obs.cov<- moon_pred %>%  # now there's NA's that we might need to double check.
  select(site,phase,wk, l.illum, fraction) %>% 
  group_by(site, wk) %>% 
  summarize(av_phase= mean(phase)) %>% 
  pivot_wider(names_from = wk, values_from = av_phase)

obs.cov <- replace(obs.cov, is.na(obs.cov), 0) # replace NAs with 0

obs.cov<- obs.cov %>%  select(!c("site", "34")) # remove site and last week like lano_js

obs.cov2 <- moon_pred %>% select(site, phase, wk, l.illum) %>% # making another obs.cov with illumination
  group_by(site, wk) %>%
  summarise(av_m.ill = mean(l.illum)) %>%
  pivot_wider(names_from = wk, values_from = av_m.ill) %>% 
ungroup()

obs.cov2<- replace(obs.cov2, is.na(obs.cov2),0)
obs.cov<- obs.cov2 %>% select(!c("site","34")) # removes site and week 34 columns
         
write.csv(obs.cov,file = 'data_for_analysis/bat_pop_analysis/obs.cov.csv',
          row.names = F) 
write.csv(obs.cov2,file = 'data_for_analysis/bat_pop_analysis/obs.cov2.csv',
          row.names = F)







# kpro_ data --------------------------------------------------------------

# load

kpro_2021_bat<-read.csv(file = 'data_for_analysis/kpro2021_v1.csv')

kpro_2021_bat$jday<-lubridate::yday(kpro_2021_bat$DATE) # julian day

# date time col. 
datetime<-paste(kpro_2021_bat$DATE, kpro_2021_bat$TIME)#merge date and time
datetime.parse<-lubridate::ymd_hms(datetime) # parse as date time
kpro_2021_bat$date_time<-datetime.parse # add to data. 




# kpro_2021_bat$hora<- as.POSIXct(kpro_2021_bat$TIME, format = "%H:%M:%S") # possibly remove. 

# kpro_2021_bat$roun.min <- round(kpro_2021_bat$hora, units = "mins")
kpro_2021_bat$rmins<-round(kpro_2021_bat$date_time, units="mins") #rounds to the nearest min



# building matrix of days 
sites<-unique(kpro_2021_bat$site)

for (sites in i:length(sites)){
  print(sites)
}

effort_days <- kpro_2021_bat %>%
  group_by(site) %>%
  summarise(
    stard = min(noche),
    endd = max(noche),
    eff.days = as.numeric(difftime(max(noche), min(noche), units = "days"))
  )



bmat <- kpro_2021_bat %>%
  group_by(site, AUTO.ID.) %>% # I don't include year because it is a single year
  count(wk, .drop = FALSE) %>%  # we might have to include the argument .drop=false to count the NAs and the zeros
  pivot_wider(names_from = wk, values_from = n) %>%
  ungroup()
summary(bat_js)

bm <- kpro_2021_bat %>%
  group_by(site, AUTO.ID.) %>% 
  count(jday, .drop = FALSE) %>%  # we might have to include the argument .drop=false to count the NAs and the zeros
  ungroup()

bmat2 <- kpro_2021_bat %>%
  group_by(site, AUTO.ID.) %>% # 
  count(jday,wk, .drop = FALSE) %>%  # this is per day and week 
  pivot_wider(names_from = c(jday,wk), values_from = n) %>%
  ungroup()

# by species 

no_id <- bmat2[bmat2$AUTO.ID. == "NoID",] # filter data to just Unkown calls

byspecies <- function(data, species) { # function that creates one matrix for spp 
  filtered_data <- data %>%
    filter(AUTO.ID. == species) %>%
    drop_na(AUTO.ID.) %>%
    mutate(across(where(is.numeric), ~if_else(is.na(.), 0, .))) 
    # select(-c("34", "site", "SppAccp"))
  
  return(filtered_data)
}

mylu_w<-byspecies(bmat, "MYOLUC") # myly by week
mylu_d.w<-byspecies(bmat2,"MYOLUC") # mylu by day and week



# finding true zeros 






# Miller ------------------------------------------------------------------

#How to calculate number of minutes of activity per night?

# Group by site, bat species, date, and minute block, then count the number of rows
# presence_summary <- kpro_2021_bat %>%
#   group_by(site, AUTO.ID., noche, roun.min) %>%
#   summarize(presence_count = n()) %>%
#   ungroup()

presence_min<-kpro_2021_bat %>% #min of activity 
  group_by(site, AUTO.ID., noche, rmins) %>% 
  summarize(activity_min= n()) %>% 
  ungroup()
  
  
species_summary <- presence_summary %>% # number of minutes by day. 
  group_by(site, noche, AUTO.ID.) %>%
  summarize(activity_min = sum(presence_count))


#get week and jday

species_summary$wk<-week(species_summary$noche)
species_summary$jday<-yday(species_summary$noche)

sp_wk <- species_summary %>%
  group_by(site, AUTO.ID.,wk) %>% # 
  count(activity_min, .drop = FALSE) %>%  # this is per day and week 
  pivot_wider(names_from = c(wk), values_from = n()) %>%
  ungroup()


t <- species_summary %>%
  group_by(site, AUTO.ID., wk) %>%
  summarize(activity_minutes = sum(activity_min)) %>%
  pivot_wider(names_from = wk,
              values_from = activity_minutes) %>% 
  ungroup()

# Pivot the data to wide format
wide_data <- summary_data %>%
  pivot_wider(
    names_from = wk,
    values_from = total_activity_minutes,
    names_prefix = "week_"
  )













#  PLOTS ----------

ggplot(bat1, aes(jday, n, fill = treatmt)) +
  geom_col() +
  facet_wrap( ~ SppAccp, scales = "free")

ggplot(data = bat1, aes(factor(hr, levels = c(
  20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6
)), n, fill = treatmt)) +
  geom_col() +
  facet_grid(. ~ treatmt) +
  theme_classic() +
  xlab("hours") +
  ylab("Calls counts 2021")
# scale_fill_manual(values=Blues)

ggplot(data = bat2021, aes(dlt.sunset, fill=treatmt))+
  geom_histogram(position = "dodge")

ggplot(bat2021, aes(phase, hr)) +
  geom_histogram()

# graph of predictors -------------

hist(moon_pred$phase) # we want to see if there is a good spread
hist(moon_pred$fraction) # values are sparcer. 
hist(bat2021$dlt.sunset) # time since sunset 
hist(bat2021$SppAccp)






# plot Activity index -------------------------------------------------------------




ggplot(species_summary, aes(yday(noche), activity_min))+
  geom_col()+
  facet_wrap(~ AUTO.ID.,scales = "free")

ggplot(species_summary, aes(yday(noche), activity_min))+
  geom_col()+
  facet_wrap(~ site,scales = "free")

ggplot(species_summary, aes(activity_min))+
  geom_histogram(bins = 30)

ggplot(bm, aes(jday, n))+
  geom_col()+
  facet_wrap(~ site,scales = "free")



















#junk







lano <- bat2021[bat2021$SppAccp == "Lano",]
filtered_data <-
  lano[(lano$site == "IRON03") &
         (lano$wk == 25),] # I checked to see if the NA in Lano were actually missing values and to miss calculations

lano_js<-lano %>% 
  group_by(site) %>% 
  count(wk) %>% 
  pivot_wider(names_from = wk, values_from = n) %>% 
  ungroup()

lano_js<-lano_js[,-1] # remove col one that has just consecutive numbers
names(lano_js)
