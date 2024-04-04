## ---------------------------
##
## Script name:  2019 bat data explore
##
## Purpose of script: we will cleanup and analyze the bat data from 2019 to see if it has similar patterns to 2021-2022
##
## Author: Carlos LInares
##
## Date Created: 4/4/2024
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------
## 


# libraries  --------------------------------------------------------------


library(tidyverse)


# data --------------------------------------------------------------------


bat2019_raw<- read.csv(file = 'data_for_analysis/data2019_db_only/bats2019.csv', header = T, check.names = T) # this is sonobat data 


# clean_up ----------------------------------------------------------------

keep<- c("Path","Filename","HiF","LoF", "SppAccp", "Prob","X.Spp", "X.Prob","X1st")

bat2019_v1 <- bat2019_raw %>% select(all_of(keep)) # removes unnecessary columns 

bat2019_v1[bat2019_v1==""]<- NA # makes the empty spaces NAs

#bat2019_v2 <- bat2019_v1 %>%  filter(!is.na(SppAccp)) # we don't want to remove NAs yet

# time extract

bat2019_v1 <- bat2019_v1 %>%
  mutate(date_time = ymd_hms(str_extract(Filename, "\\d{8}_\\d{6}"), tz = "America/Denver")) 

# site


split_text <- str_split(bat2019_v1$Path, "\\\\", simplify = TRUE)

# Extract the 4th element (index 4)
site_names <- split_text[,5]

bat2019_v1$site<-str_extract(bat2019_v1$Path, "^[A-Za-z]{2}\\d+")
