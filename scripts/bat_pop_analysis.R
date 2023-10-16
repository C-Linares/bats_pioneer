
# Carlos Linares
# 
# Script para hacer modelos de población con los datos de 2021.

# 
# #
# libraries 
library(tidyverse)
library(lubridate)
library(magrittr)
library(elevatr)

# Data load-------------####

bat2021<-read.csv('data_analysis/bat2021_v2.csv',check.names = T) # we load the data we created in the clean up script. This data has several variables added to a clean up raw data is the product of clean_up script

# make dates-------------
bat2021$datetime<-ymd_hms(bat2021$datetime) # makes dates as dates 
bat2021$noche<-ymd(bat2021$noche)# makes noche as date. noche is the date for all the calls that come from the same night instead of having half one day and the others the other day.
bat2021$jday<-yday(bat2021$noche)
bat2021$wk<-week(bat2021$noche)# we need to calculate a week column.
bat2021$yr<-year(bat2021$noche) #year added for when we will have multiple years.

#sites
unique(bat2021$site)# tell us what sites we have

# Matrix building  ------------------

bat1<- bat2021 %>% # bat matrix
  group_by(site,treatmt) %>% 
  count(SppAccp, jday, hr) %>% 
  ungroup()

# Predictors matrix
moon_pred<-read.csv('data_analysis/moon_pred.csv')


bat_pred<- bat2021 %>% 
  group_by(site,SppAccp,jday) %>% 
  select(phase,dlt.sunset,trmt_bin,wk,yr) %>% 
  ungroup()

# Wrangling for modeling

bat_js<- bat2021 %>% 
  group_by(site,SppAccp) %>% # I don't include year because it is a single year 
  count(wk) %>%
  pivot_wider(names_from = wk, values_from = n) %>%
  ungroup()

t<- bat_pred %>% 
  group_by(site, SppAccp, wk) %>% 
  
  pivot_wider(names_from= wk, values_from = c(phase, dlt.sunset, trmt_bin))













#some PLOTS ##-----------#####

ggplot(bat1, aes(jday, n, fill=treatmt))+
  geom_col()+
  facet_wrap(~ SppAccp, scales = "free")

ggplot(data = bat1, aes(factor(hr, levels = c(20,21,22,23,0,1,2,3,4,5,6)), n, fill=treatmt))+
  geom_col()+
  facet_grid(.~treatmt)+
  theme_classic()+
  xlab("hours")+
  ylab("Calls counts 2021")
  # scale_fill_manual(values=Blues)

ggplot(data = bat2021, aes(dlt.sunset, fill=treatmt))+
  geom_histogram(position = "dodge")

ggplot(bat2021, aes(phase, hr)) +
  geom_histogram()

# graph of predictors -------------

hist(bat2021$phase) # we want to see if there is a good spread
hist(bat2021$fraction) # values are sparcer. 
hist(bat2021$dlt.sunset) # time since sunset 


