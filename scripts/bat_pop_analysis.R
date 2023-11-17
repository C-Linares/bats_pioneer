
# Carlos Linares
# 
# Population modelling for 2021 data.
# Script para hacer modelos de población con los datos de 2021.

# 
# #
# libraries 
library(tidyverse)
library(lubridate)
library(magrittr)
library(elevatr)
library(lubridate)

# Data load-------------####

bat2021<-read.csv('data_analysis/bat2021_v2.csv',check.names = T) # we load the data we created in the clean up script. This data has several variables added to a clean up raw data is the product of clean_up script

# make dates-------------
bat2021$datetime <- ymd_hms(bat2021$datetime) # makes dates as dates
bat2021$noche <-
  ymd(bat2021$noche)# makes noche as date. noche is the date for all the calls that come from the same night instead of having half one day and the others the other day.
bat2021$jday <- yday(bat2021$noche)
bat2021$wk <-
  week(bat2021$noche)# we need to calculate a week column.
bat2021$yr <-
  year(bat2021$noche) #year added for when we will have multiple years.

#sites
unique(bat2021$site)# tell us what sites we have

# Matrix building  ------------------

bat1 <- bat2021 %>%
  group_by(site, treatmt) %>%
  count(SppAccp, jday, hr) %>%
  ungroup()

bat_js <- bat2021 %>%
  group_by(site, SppAccp) %>% # I don't include year because it is a single year
  count(wk) %>%
  pivot_wider(names_from = wk, values_from = n) %>%
  ungroup()

lano_js <- bat_js[bat_js$SppAccp == "Lano",] # filter data to just Lano.
lano_js <-
  lano_js[, c(-1, -2, -12)] # remove site, sp, last week because it has few data points.

t<- lano_js[is.na(lano_js)]<-0

# we filter by just one sp.
# what species has the more calls. Lano 

ggplot(bat2021,aes(x=SppAccp, fill=treatmt))+
  geom_bar()
  



#rename week columns as numbers just numbers 

write.csv(bat_js,file = 'data_analysis/bat_pop_analysis/bat_js.csv') 

#write single species df

write.csv(lano_js,file = 'data_analysis/bat_pop_analysis/lano_js.csv') 


# site level cov -------------

s.l.c <- bat2021 %>% dplyr::distinct(site,treatmt,elevation)

write.csv(s.l.c,file = 'data_analysis/bat_pop_analysis/slc.csv') 

# obs cov ----------------

# moon phase not changes by site so we just need one by site. 



moon_pred<-read.csv('data_analysis/moon_pred.csv') #loads the data from the moon_pred script

obs.cov<- moon_pred %>%  # now there's NA's that I am not sure where they come from
  select(site,phase,wk, l.illum, fraction) %>% 
  group_by(site, wk) %>% 
  summarize(av_phase= mean(phase)) %>% 
  pivot_wider(names_from = wk, values_from = av_phase)


obs.cov2<-moon_pred %>% select(site, phase, wk, l.illum) %>% 
  group_by(site,wk) %>% 
  summarise(av_m.ill=mean(l.illum)) %>% 
  pivot_wider(names_from = wk, values_from = av_m.ill)


write.csv(obs.cov,file = 'data_analysis/bat_pop_analysis/obs.cov.csv') 










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










# junk ----------
# 



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
