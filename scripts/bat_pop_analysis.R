
# Carlos Linares
# 
# Script para empezar a trabajar con Jen Cruz on los datos de 2021
# 
# #
# libraries 
library(tidyverse)
library(lubridate)
library(magrittr)
# Data load-------------

bat2021<-read.csv('data_analysis/bat2021_v2.csv',check.names = T) # we load the data we created in the clean up script. this data has several variables added to a clean up raw data


bat2021$datetime<-ymd_hms(bat2021$datetime) # makes dates as dates 
bat2021$noche<-ymd(bat2021$noche)# makes dates as dates
bat2021$jday<-yday(bat2021$noche)

unique(bat2021$site)# tell us what sites we have

# we want a matrix of species counts by jday 

bat1<- bat2021 %>% 
  group_by(site,treatmt) %>% 
  count(SppAccp, jday, hr) %>% 
  ungroup()

# no se como incluir el treatment en bat1

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


# graph of bats counts vs jday 

hist(bat2021$phase) # we want to see if there is a good spread
hist(bat2021$fraction) # values are sparcer. 
hist(is.na(bat2021$sunrise)) # check why there's NA's
hist(bat2021$sunset) # check why there' NA's or not numeric.

hist(bat2021$dlt.sunset)

