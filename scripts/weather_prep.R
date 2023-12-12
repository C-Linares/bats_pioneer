
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



#----------- load data
wetfls<-list.files('data_for_analysis/weather/',pattern = "*.csv",full.names = T) # list csvs

wetdb <- lapply(wetfls, function(wetfls) {
  read.csv(wetfls, skip =  7) # Skips 7 rows (header + 4 data rows)
})

library(purrr)

wetdb <- purrr::map_df(wetfls, ~ read.csv(.x, skip = 7))
