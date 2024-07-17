
#----------------- weather data preparation
# ------- purpose : generate weekly estimates for temp and rain. 
# 
# author: Carlos Linares
# date: dec 2023
# 

#------------  load library 

library(tidyverse)
library(magrittr)
library(data.table)
library(purrr)




#----------- load data ---------------
#Data was obtainde from https://mesowest.utah.edu/ statition KSUN in Idaho. 
#list data
wetfls<-list.files('data_for_analysis/weather/',pattern = "^KSUN.*\\.csv$",full.names = T)

# this creates the data frame list but not the data base. 
wetdb <- lapply(wetfls, function(wetfls) {
  read.csv(wetfls, skip =  7) # Skips 7 rows (header + 4 data rows)
})


wetdb <- purrr::map_df(wetfls, ~ read.csv(.x, skip = 7)) # new function from purr to read multiple files 
str(wetdb)

# data cleaning --------------

wetdb<- setnames(wetdb, old = c("X.1","X"), new = c("date.time", "station")) # change col name to date

wetdb$date.time<- mdy_hm(wetdb$date.time, tz = "America/Denver")

wetdb$wk<- week(wetdb$date.time)

wetdb$hr<-hour(wetdb$date.time)

# wetdb <- wetdb %>% select(-"Millimeters")# remove empty colum

# just date column

wetdb$date<- date(wetdb$date.time)

# -- calculate the mean temp and rain -----

daily_averages <- wetdb %>%
  
  group_by(date) %>%
  summarize(
    avg_temperature = mean(Celsius, na.rm = TRUE),
    avg_wind_speed = mean(m.s, na.rm = TRUE)
  )
  
# calculate night averages

nigh_averages <- wetdb %>%
  # Filter for nighttime records (example: 6 PM to 6 AM)
  filter(hr >= 18 | hr < 6) %>%
  # Group by date
  group_by(date) %>%
  # Calculate nightly averages
  summarize(
    avg_temperature = mean(Celsius, na.rm = TRUE),
    avg_wind_speed = mean(m.s, na.rm = TRUE)
  )

#filter clouds just at night. 
cls<- wetdb %>% 
  filter(hr >= 18 | hr < 6) %>% 
  mutate(cld= str_extract(cls$code, "\\d$"))

# extract the last digit of the cloud code. This correspond to the cloud magnitude. 
cls.code<-str_extract(cls$code, "\\d$")

write.csv(daily_averages, file = "data_for_analysis/weather/dailyavg.csv",row.names = F)
write.csv(daily_averages, file = "data_for_analysis/weather/nigh_averages.csv",row.names = F)




#write data

write.csv(wetdb,file = "data_for_analysis/weather/rawweather.csv",row.names = F)



# figure ------------------------------------------------------------------

ggplot(wetdb, aes(x = date.time, y = Celsius)) +
  geom_point(color = "blue") +
  labs(title = "Time Series of Temperature",
       x = "Date and Time",
       y = "Temperature (Celsius)") +
  theme_minimal()

ggplot(daily_averages, aes(x= date, y=avg_temperature ))+
  geom_point()+
  labs(title = "time series temp",
       x="date",
       y= "temp C")


ggplot(cls, aes(x = cld)) +
  geom_bar(stat = "count", fill = "blue", color = "black") + # there's almost no clouds at night. 
  theme_minimal() +
  labs(title = "", x = "", y = "Frequency")



# junk --------------------------------------------------------------------

# week weather.
#   possibly not necessary 
mtempwind<- wetdb %>%
  group_by(wk) %>%
  summarize(mean_tem = mean(Celsius, na.rm = T), mean_wind = mean(m.s, na.rm = T))# there are NA's I don't know from where. I know now that we need to include the na.rm=T argument so there are no 

# Filter and save rows with NAs
na_rows <- wetdb[!complete.cases(wetdb), ]

