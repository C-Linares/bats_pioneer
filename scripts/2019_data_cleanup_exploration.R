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


split_text <- str_split(bat2019_v1$Path, "\\\\", simplify = TRUE) # split the path string

# Extract the 4th element (index 4)
site_names <- split_text[,5] # keep the string for site

bat2019_v1$site<-site_names # add column to data. 

#noche

bat2019_v1$noche <-
  if_else(hour(bat2019_v1$date_time) < 9, # if it is less than 9 put the date of the previous day
          true =  (date(bat2019_v1$date_time) - ddays(1)),
          false = date(bat2019_v1$date_time))

#week

bat2019_v1$wk<-week(bat2019_v1$noche)

bmat <- bat2019_v1 %>%
  group_by(site, SppAccp) %>% # I don't include year because it is a single year
  count(wk, .drop = FALSE) %>%  # we might have to include the argument .drop=false to count the NAs and the zeros
  pivot_wider(names_from = wk, values_from = n) %>%
  ungroup()
