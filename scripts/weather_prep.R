
#----------------- weather data preparation
# ------- purpose : generate weekly estimates for temp and rain. 
# 
# author: Carlos Linares
# date: dec 2023
# 

#------------  load library 
#
library(tidyverse)
library(magrittr)
library(data.table)
library(purrr)




#----------- load data ---------------
wetfls<-list.files('data_for_analysis/weather/',pattern = "*.csv",full.names = T) # list csvs

# this creates the data frame list but not the data base. 
wetdb <- lapply(wetfls, function(wetfls) {
  read.csv(wetfls, skip =  7) # Skips 7 rows (header + 4 data rows)
})


wetdb <- purrr::map_df(wetfls, ~ read.csv(.x, skip = 7)) # new function from purr to read multiple files 
# data cleaning --------------

wetdb<- setnames(wetdb, old = "X.1", new = "date.time") # change col name to date

wetdb$date.time<- mdy_hms(wetdb$date.time)

wetdb$wk<- week(wetdb$date.time)

wetdb <- wetdb %>% select(-"Millimeters")# remove empty colum

# -- calculate the mean temp and rain -----

wtemp <- wetdb %>%
  group_by(wk) %>%
  summarize(mean_tem = mean(Celsius), mean_wind = mean(m.s))# there are NA's I don't know from where.

# Filter and save rows with NAs
na_rows <- wetdb[!complete.cases(wetdb), ]
