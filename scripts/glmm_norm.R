# ---------------------------
##
## Script name:  glmm_norm
##
## Purpose of script: run models on the normalized data
##
## Author: Carlos Linares, \
## Date Created: 12/13/2024
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes: 
##   
##
## ---------------------------
## # inputs ------------------------------------------------------------------
#   data_for_analysis/prep_for_glmm/normalized_bm.csv  #product of the prep_for_glmm.R script. summary daily of bat call counts


# outputs ----------------------

# figures and model data


# libraries  --------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(lme4)
library(sjPlot)
library(ggeffects)
library(car)
library(glmmTMB)
library(corrplot)
library(effects)
library(reshape2)


# data --------------------------------------------------------------------


normalized_bm<-read_csv('data_for_analysis/prep_for_glmm/normalized_bm.csv')
head(normalized_bm)
normalized_bm$date<- ymd(normalized_bm$noche)


moon<-read.csv('data_for_analysis/moon_pred.csv')
moon$date<- as_date(moon$date)
moon.adj<-moon %>% mutate(
  phase = ifelse(above_horizon==FALSE,0,phase),
  fraction= ifelse(above_horizon==FALSE,0,fraction),
  l.illum= ifelse(above_horizon==FALSE,0,l.illum)
)
head(moon.adj)
moon.adj$date<-ymd(moon.adj$date)

# dayly average of l.illum

l.illum <- moon.adj %>% group_by(date) %>% summarise(l.illum = mean(l.illum))
l.illum$date<-ymd(l.illum$date)
head(l.illum)

weather<-read.csv('data_for_analysis/weather/nigh_averages.csv', header = T) #load nightly averages
weather$date<-as_date(weather$date)
head(weather)
weather$date<-ymd(weather$date)


#merge data by date
bm<-left_join(normalized_bm, l.illum, by='date')
bm<-left_join(bm, weather, by='date')
head(bm)



#correlation matrix just numeric variables  (not including date)
#how do I first remove the NA falues from bm? 




cor_matrix <- cor(bm[, sapply(bm, is.numeric)])
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", tl.cex = 0.7, tl.col = "black", tl.srt = 45)
  

m1.1 <- glmmTMB(normalized_activity ~ + jday + I(jday^ 2) + l.illum +  avg_wind_speed + avg_temperature,
  data = bm,
  gaussian(link = "identity"))
summary(m1.1)
