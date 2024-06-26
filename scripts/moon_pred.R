

#script to calculate the average and median moon phase for use a predictor in a population model
#we will be using the suncalc package

library(suncalc)
library(sf)
library(ggplot2)
library(lunar)
library(tidyverse)

#we need the coordinates for the sites

sts<-read.csv("data_for_analysis/sites_coordinates.csv")
sts<-sts[1,]


# load bat data 2021 for dates. 

kpro_2021_bat<-read.csv('data_for_analysis/kpro2021_v1.csv')
bm<-read.csv('data_for_analysis/data_glmm/bat_counts.csv', header = T) # bat counts by jday
# unique(bm$site)
#get dates from bat2021

datesite<-bm %>% select(site,jday)
datesite$date<-as_date(bm$jday, origin= "2021-01-01")
datesite$jday<-NULL

# generate dates from 2021 to 2023 just the summer months. 

start_date <- as.Date("2020-05-01")
end_date <- as.Date("2023-09-30")

# Generate the sequence of dates
all_dates <- seq.Date(from = start_date, to = end_date, by = "day")
# Filter the dates to include only May to September
dates<- all_dates[format(all_dates, "%m") %in% c("05", "06", "07", "08", "09")]
dates<-as.data.frame(dates)
colnames(dates)[1]<-"date"


# merge tdates# merge the dates with the site and the coordinates. 

t1<-merge(dates, sts, by=NULL)

# moon illumination, times, and pos ----------------------

t1<-left_join(datesite,sts, by = "site")
t1<-rename(t1, date=noche) # rename noche as date but remember this when rejoining the other data.
t1<-t1 %>% select(date,lat,lon)

moon_pos<- unique(getMoonPosition(data = t1))
moon_pos$date<-as.Date(moon_pos$date)


t1$date<-as.Date(t1$date)
moon_times<-unique(getMoonTimes(data = t1))

t2 <- as.Date(unique(t1[, 1]), tz = "America/Denver")
moon_illumination <- getMoonIllumination(date = t2)


moon_pred <- left_join(left_join(moon_pos, moon_times), moon_illumination)


moon_pred<-moon_pred %>% select(date, lat, lon, altitude, rise,set,phase, fraction, parallacticAngle, angle)

moon_pred$above_horizon <- moon_pred$altitude > 0


#lunar illumination calculated with lunar package


attr(moon_pred$date, "tzone") <- "UTC" # make them UTM so the calculation is ok
# myshift <- as.numeric(format(moon_pred$Datetime, "%H")) - 12 # read https://stackoverflow.com/questions/71757462/calculate-lunar-illumination-using-the-lunar-package-in-r
moon_pred$l.illum<-lunar.illumination(moon_pred$noche)



# make moon pred 0 if not above the horizon.

moon.adj<-moon_pred %>% mutate(
  phase = ifelse(above_horizon==FALSE,0,phase),
  fraction= ifelse(above_horizon==FALSE,0,fraction),
  l.illum= ifelse(above_horizon==FALSE,0,l.illum)
)


colnames(moon_pred)[1]<-"date"




write.csv(moon_pred,file = 'data_for_analysis/moon_pred.csv', row.names = F) # ran once in case of needing to rewrite the data



#-----------------plots-----------------
ggplot(moon.adj, aes(x=altitude, y=phase))+
  geom_point()+
  geom_vline(xintercept = 0)

ggplot(moon.adj, aes(x=altitude, y=l.illum))+
  geom_point()+
  geom_vline(xintercept = 0)

hist(moon_pred$phase)
hist(moon_pred$fraction)
hist(moon_pred$l.illum)

plot(moon_pred$phase ~ moon_pred$l.illum)# no linear relationship between moon illumination and moon phase. 


# Create a new column indicating moon above horizon
moon_pos$above_horizon <- moon_pos$altitude > 0

# Create the plot
ggplot(moon_pos, aes(x = yday(moon_pos$date), y = altitude, color = above_horizon)) +
  geom_line() +
  # scale_x_datetime(breaks = date_breaks("12 hours")) +  # Adjust breaks as needed
  labs(x = "Time", y = "Moon Altitude", color = "Moon Above Horizon?") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray")  # Optional: Show threshold lin
