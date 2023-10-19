

#script to calculate the average and median moon phase for use a predictor in a population model
#we will be using the suncalc package

library(suncalc)
library(sf)
library(ggplot2)
library(lunar)
library(tidyverse)

#we need the coordinates for the sites

sts<-read.csv("data_analysis/sites_coordinates.csv")
bat2021<-read.csv('data_analysis/bat2021_v2.csv',check.names = T)
#get dates from bat2021

datesite<-bat2021 %>% select(site,noche)

# moon illumination, times, and pos ----------------------

t1<-left_join(datesite,sts, by = "site")
t1<-rename(t1, date=noche) # rename noche as date but remember this when rejoining the other data.
t1<-t1 %>% select(date,lat,lon)

moon_pos<- unique(getMoonPosition(data = t1))
moon_pos$date<-as.Date(moon_pos$date)

t1$date<-as.Date(t1$date)
moon_times<-unique(getMoonTimes(data = t1))

t2<-t1[,1]
t2<-unique(t2)
t2<-as.Date(t2)
attr(t2, "tzone")<-"America/Denver"


moon_illumination<-getMoonIllumination(date = t2)


t<- left_join(moon_pos,moon_times)
moon_pred<- left_join(t,moon_illumination)

moon_pred<-moon_pred %>% select(date, lat, lon, altitude, rise,phase, fraction)

#lunar illumination calculated with lunar package

moon_pred$l.illum<-lunar.illumination(moon_pred$date, shift = -6)
moon_pred<-rename(moon_pred, noche=date)#re rename date as noche so it can be merge with the other tables

# add sites labels back
moon_pred<-left_join(moon_pred,sts)

#add week 

moon_pred$wk<-week(moon_pred$noche)

# adding if above horizon marker. 




# write.csv(moon_pred,file = 'data_analysis/moon_pred.csv') # ran once in case of needing to rewrite the data



#-----------------plots-----------------
ggplot(moon_pred, aes(x=altitude, y=phase))+
  geom_point()+
  geom_vline(xintercept = 0)

ggplot(moon_pred, aes(x=altitude, y=l.illum))+
  geom_point()+
  geom_vline(xintercept = 0)

hist(moon_pred$phase)
hist(moon_pred$fraction)
hist(moon_pred$l.illum)

plot(moon_pred$phase ~ moon_pred$l.illum)# no linear relationship between moon illumination and moon phase. 


  