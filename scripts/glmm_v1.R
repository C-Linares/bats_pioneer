# ---------------------------
##
## Script name:  glmm_v1
##
## Purpose of script: Start running models for the 2021 data with glmm
##
## Author: Carlos Linares
##
## Date Created: 05/21/2024
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes: need to anotate what script produced each data set inputed. 
##   
##
## ---------------------------
## 


# libraries  --------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(lme4)
library(sjPlot)
library(ggeffects)
library(car)
library(glmmTMB)


# data --------------------------------------------------------------------


bm<-read.csv('data_for_analysis/data_glmm/bat_counts.csv', header = T) # bat counts by jday
bm$datefromyday<-as_date(bm$jday, origin= "2021-01-01")
filtered_bm <- bm[bm$AUTO.ID. != "Noise", ] # if noise is not filtered there is more records for dark sites. 
colnames(filtered_bm)[3]<-"sp" # change Auto.ID to sp

btrait<-read.csv('data_for_analysis/Bat_trait.csv', header = T)
btrait$Species<-toupper(btrait$Species) # makes all caps
batnames<-read.csv('data_for_analysis/Species_bats.csv')
colnames(batnames)[3]<-"sp"
colnames(btrait)[2]<-"sp"

elev<-read.csv('data_for_analysis/elev/elevation.csv', header = T)
elev<-elev %>% rename("site" = "name")
elev$site <-tolower(elev$site)
#Use gsub to replace 'viz' with 'vizc' in the 'site' column of df1
elev$site <- gsub("viz(\\d{2})", "vizc\\1", elev$site)


ndvi<-read.csv('data_for_analysis/NDVI/NDVI_of_rip2021.csv', header = T)
ndvi<- ndvi %>% dplyr::select(-c("ele", "time", "magvar", "geoidheigh", "dgpsid") ) # remove unncessary cols
ndvi<-ndvi %>% rename("site" = "name")
ndvi<-ndvi %>% rename("ndvi_mean" = "X_mean")
ndvi$site<-tolower(ndvi$site)
#Use gsub to replace 'viz' with 'vizc' in the 'site' column of df1
ndvi$site <- gsub("viz(\\d{2})", "vizc\\1", ndvi$site)

weather<-read.csv('data_for_analysis/weather/dailyavg.csv', header = T)
weather$date<-as_date(weather$date)

# moon clean-up
moon<-read.csv('data_for_analysis/moon_pred.csv')
# change noche to date
moon<-moon %>% rename("date" = "noche")
moon$date<- as_date(moon$date)
moon$site <-tolower(moon$site)

moon<- moon %>% dplyr::select(-c("X.1","X") ) # remove unncessary cols




#merge 
namestraits<-left_join(btrait, batnames, by= "sp")
colnames(namestraits)[2]<-"sp4"
traits<- namestraits %>% dplyr::select(c("sp4","Six.letter.species.code","Mass","Aspect","PeakFreq","Loading" ) )
colnames(traits)[2]<-"sp"
rm(namestraits) 

filtered_bm<- left_join(filtered_bm, traits, by="sp")

bm2<-left_join(filtered_bm, elev, by="site" )
bm2 <- left_join(bm2, ndvi, by = "site") %>%
  select(-c("X", "time", "buff_area.x", "buff_area.y"))

colnames(bm2)[7]<-"date"

bm2$date<-as_date(bm2$date)
bm2<- left_join(bm2, weather, by="date")


t<- left_join(bm2, moon, by=c("site","date"))


# Identify rows with NA values
rows_with_na <- bm2[rowSums(is.na(bm2)) > 0, ] # some wind and temp have NAs because we are missing agusot 2021 weather data.


# scale data 

bm2$elev_mean_s<- scale(bm2$elev_mean, center = T, scale = T)
bm2$ndvi_mean_s<- scale(bm2$ndvi_mean)
bm2$percent_s<-scale(bm2$percent)
bm2$jday_s<-scale(bm2$jday)

# explore data ------------------------------------------------------------

# Plot the distribution of the count data
ggplot(filtered_bm, aes(x = n)) +
  geom_histogram(binwidth = 40, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Bat Calls", x = "Number of Calls", y = "Frequency")

ggplot(filtered_bm, aes(x=jday, y=n, col=treatmt))+
  geom_point()+
  facet_wrap("site")


ggplot(bm2, aes(x = elev_mean)) +
  geom_histogram(fill = "blue", color = "black") +
  theme_minimal() 

ggplot(bm2, aes(x = percent_s)) +
  geom_histogram( fill = "blue", color = "black") +
  theme_minimal() 

# models -------------------------------------------------------------------

# Fit a Poisson model intercept model. 
m1 <- glmer(n ~ (1|site),
                     data = bm2, 
                     family = poisson)
exp(m1@beta)# base line for bats at all sites. 

summary(m1)



m1.1<- glmer(n ~ trmt_bin +(1|site),
           data = bm2, 
           family = poisson)
exp(m1.1@beta)
summary(m1.1)



m1.2<- glmer(n ~ trmt_bin * jday_scaled + percent_s +elev_mean_s+  (1|site),
                   data = bm2, 
                   family = poisson)

plot_model(m1.2)

exp(m1.2@beta)
summary(m1.2)

plot(fitted(m1.2), residuals(m1.2), main = "Residuals vs Fitted", 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(residuals(m1.2), main = "Normal Q-Q Plot")
qqline(residuals(m1.2), col = "red")

# Histogram of residuals
hist(residuals(m1.2), breaks = 30, main = "Histogram of Residuals", 
     xlab = "Residuals")

# Scale-Location Plot
plot(fitted(m1.2), sqrt(abs(residuals(m1.2))), main = "Scale-Location Plot",
     xlab = "Fitted values", ylab = "Square Root of |Residuals|")
abline(h = 0, col = "red")

# Calculate VIF
vif(m1.2)



# quasipoisson -------------------------------------------------------



