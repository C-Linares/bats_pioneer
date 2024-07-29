
# Script: bat_prep_analysis.R
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

slc<-species_summary %>% distinct(site,trmt_bin)
slc<-slc[,-1]
write.csv(slc,file = 'data_for_analysis/bat_pop_analysis/slc.csv',row.names = F)



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







# kpro_data --------------------------------------------------------------

# load

kpro_2021_bat<-read.csv(file = 'data_for_analysis/kpro2021_v1.csv') # data loading is the product of the script cleanup_script_v2




kpro_2021_bat$date_time<-ymd_hms(kpro_2021_bat$date_time) # for some reason 22 fail to parse. 
# kpro_2021_bat$jday<-lubridate::yday(kpro_2021_bat$DATE) # julian day

#date time col. 
datetime<-paste(kpro_2021_bat$DATE, kpro_2021_bat$TIME)#merge date and time
datetime.parse<-lubridate::ymd_hms(datetime) # parse as date time
kpro_2021_bat$date_time<-datetime.parse # add to data. 


kpro_2021_bat$rmins<-round(kpro_2021_bat$date_time, units="mins") #rounds to the nearest min
kpro_2021_bat$rmins<-round_date(kpro_2021_bat$date_time, unit = "minute")

# building matrix of days 


effort_days <- kpro_2021_bat %>%
  group_by(site) %>%
  summarise(
    stard = min(noche),
    endd = max(noche),
    eff.days = as.numeric(difftime(max(noche), min(noche), units = "days"))
  )

effort_hrs <- kpro_2021_bat %>%
  group_by(site, noche, jday) %>%
  summarise(stard = min(date_time), endd = max(date_time)) %>%
  mutate(eff.hrs = time_length(endd - stard, unit = "hours")) %>% 
  mutate(wk=week(noche)) # calculates the week too. 

write.csv(effort_hrs,file = 'data_for_analysis/data_glmm/effort_hrs.csv',
          row.names = F)

# calls by week

bmat <- kpro_2021_bat %>% # count of calls 
  group_by(site, AUTO.ID.) %>% # I don't include year because it is a single year
  count(wk, .drop = FALSE) %>%  # we might have to include the argument .drop=false to count the NAs and the zeros
  pivot_wider(names_from = wk, values_from = n) %>%
  ungroup()

# calls by day 

bm <- kpro_2021_bat %>% # 
  group_by(site, AUTO.ID.) %>% 
  count(jday, .drop = FALSE) %>%  
  ungroup()

# add treatment 
litsites<-c("iron01","iron03","iron05","long01","long03")


bm$treatmt<-ifelse(bm$site %in% litsites , "lit", "dark") # this makes a treatment variable.

bm$trmt_bin<- ifelse(bm$treatmt== "lit", 1, 0)
summary(bm) #check for NAs
write.csv(bm,file = 'data_for_analysis/data_glmm/bat_counts.csv',
          row.names = F)

# calls by day and week 
bmat2 <- kpro_2021_bat %>% 
  group_by(site, AUTO.ID.) %>% # 
  count(jday,wk, .drop = FALSE) %>%  # this is per day and week 
  # pivot_wider(names_from = c(jday,wk), values_from = n) %>%
  ungroup()

# finding true NAs


# by species 

no_id <- bmat2[bmat2$AUTO.ID. == "NoID",] # filter data to just Unknown calls

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










# Miller ------------------------------------------------------------------

#How to calculate number of minutes of activity per night?

# # Group by site, bat species, date, and minute block, then count the number of rows
# presence_summary <- kpro_2021_bat %>%
#   group_by(site, AUTO.ID., noche, rmins) %>%
#   summarize(presence_count = n()) %>%
#   ungroup()

presence_min<-kpro_2021_bat %>% #min of activity 
  group_by(site, AUTO.ID., noche, rmins) %>% 
  summarize(activity_min= n()) %>%  #calculate the num of min activte
  ungroup()
  
  
species_summary <- presence_min %>% # number of minutes by day. 
  group_by(site, noche, AUTO.ID.) %>%
  summarize(activity_min = sum(activity_min))

#treatment 

litsites<-c("iron01","iron03","iron05","long01","long03")


species_summary$treatmt<-ifelse(species_summary$site %in% litsites , "lit", "dark") # this makes a treatment variable.

species_summary$trmt_bin<- ifelse(species_summary$treatmt== "lit", 1, 0)


#get week and jday

species_summary$wk<-week(species_summary$noche)
species_summary$jday<-yday(species_summary$noche)

# standardize by effort hrs.

species_summary<-left_join(effort_hrs,species_summary, by = c("site","noche","wk"))

# AI standardize by hrs sampled. 
zeros <- filter(species_summary, eff.hrs <= 0) # some obs have 0 hrs of effort

species_summary$AI.st<-species_summary$activity_min/species_summary$eff.hrs # need to fix the efforst that have less than an hour of effort see(zeros)


bat.x.wk <- species_summary %>%
  group_by(site, AUTO.ID., wk) %>%
  summarize(activity_minutes = sum(activity_min)) %>%
  pivot_wider(names_from = wk,
              values_from = activity_minutes) %>% 
  ungroup()

# by sp

lit_brw<-byspecies(bat.x.wk,"MYOLUC") # made with AI
lit_brw <- lit.brw[, !(names(lit.brw) %in% c("site", "AUTO.ID.","25","34"))] #filter out cols

write.csv(x = lit_brw,file = "data_for_analysis/bat_pop_analysis/lit_brw.csv")


write.csv(species_summary, file = "data_for_analysis/data_glmm/AI_sp_summary.csv")








#  PLOTS ----------

# sonobat output plots
ggplot(bat1, aes(jday, n, fill=treatmt))+
  geom_col()+
  facet_wrap(~ treatmt)+
  labs(title="sonobat, #calls")

ggplot(bat1, aes(jday, n, fill=treatmt))+
  geom_col()+
  facet_wrap(~ site)+
  labs(title="sonobat, #calls")

ggplot(bat1, aes(jday, n, fill = treatmt)) +
  geom_col(position="dodge") +
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



# plot kpro call counts

# Filter the bm data set to exclude Noise rows
filtered_bm <- bm[bm$AUTO.ID. != "Noise", ] # if noise is not filtered there is more records for dark sites. 

ggplot(filtered_bm, aes(jday, n, col=treatmt))+
  geom_point()+
  facet_wrap(~ treatmt)+
  labs(title = "n calls all bat sp (kpro output)")+
  geom_vline(xintercept = 180, linetype = "dashed", color = "red")

ggplot(filtered_bm, aes(jday, n, col=treatmt))+
  geom_point()+
  facet_wrap(~ site)+
  labs(title = "kpro # calls")+
  geom_vline(xintercept = 180, linetype = "dashed", color = "red")


ggplot(filtered_bm, aes(jday, n, col=treatmt))+
  geom_point(alpha=.5)+
  facet_wrap(~ AUTO.ID.,scales = "free")+
  labs(title = "kpro # calls")+
  geom_vline(xintercept = 180, linetype = "dashed", color = "red")



# plot Activity index -------------------------------------------------------------
filter.species_summary <- species_summary[species_summary$AUTO.ID. != "Noise", ]

ggplot(filter.species_summary, aes(yday(noche), activity_min, fill=treatmt))+
  geom_col(position = "dodge")+
  facet_wrap(~ treatmt,)+
  labs(title="#min active")

ggplot(filter.species_summary, aes(yday(noche), AI.st, fill=treatmt))+
  geom_col(position = "dodge")+
  facet_wrap(~ treatmt,scales = "free")+
  labs(title="#min standardize /eff.hrs")

ggplot(filter.species_summary, aes(yday(noche), activity_min, fill=treatmt))+
  geom_col(position = "dodge")+
  facet_wrap(~ site,scales = "free")+
  labs(title="#min active")

ggplot(filter.species_summary, aes(yday(noche), activity_min, fill=treatmt))+
  geom_col(position = "dodge")+
  facet_wrap(~ AUTO.ID.,scales = "free")+
  labs(title="activity minutes")




















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



