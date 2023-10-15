

#script to calculate the average and median moon phase for use a predictor in a population model
#we will be using the suncalc package

library(suncalc)
library(sf)
library(ggplot2)
library(lunar)

#we need the coordinates for the sites

sts<-read.csv("data_analysis/sites_coordinates.csv")

#get dates from bat2021

datesite<-bat2021 %>% select(site,noche)

t1<-left_join(datesite,sts, by = "site")
t1<-rename(t1, date=noche) # rename noche as date but remember this when rejoining the other data.
t2<-select(t1, date)
t2<-unique(t2)
t2<-as.Date(t2$date)
attr(t2, "tzone")<-"America/Denver"

moon_pos<- unique(getMoonPosition(data = t1))
moon_illumonation<-getMoonIllumination(date = t2)
moon_times<-unique(getMoonTimes(data = t1))

t<- left_join(moon_pos,moon_times)
moon_pred<- left_join(t,moon_illumonation)

moon_pred<-moon_pred %>% select(date, lat, lon, altitude, rise,phase)

#lunar illumination calculated with lunar package

moon_pred$l.illum<-lunar.illumination(moon_pred$date, shift = -6)
moon_pred<-rename(moon_pred, noche=date)#re rename date as noche so it can be merge with the other tables

write.csv(moon_pred,file = 'data_analysis/moon_pred.csv')


#-----------------plots-----------------
ggplot(moon_pred, aes(x=altitude, y=phase))+
  geom_point()+
  geom_vline(xintercept = 0)

ggplot(moon_pred, aes(x=altitude, y=l.illum))+
  geom_point()+
  geom_vline(xintercept = 0)

