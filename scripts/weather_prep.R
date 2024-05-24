
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
wetfls<-list.files('data_for_analysis/weather/',pattern = "*.csv",full.names = T) # list csvs

# this creates the data frame list but not the data base. 
wetdb <- lapply(wetfls, function(wetfls) {
  read.csv(wetfls, skip =  7) # Skips 7 rows (header + 4 data rows)
})


wetdb <- purrr::map_df(wetfls, ~ read.csv(.x, skip = 7)) # new function from purr to read multiple files 
str(wetdb)

# data cleaning --------------

wetdb<- setnames(wetdb, old = "X.1", new = "date.time") # change col name to date

wetdb$date.time<- mdy_hm(wetdb$date.time, tz = "America/Denver")

wetdb$wk<- week(wetdb$date.time)

wetdb <- wetdb %>% select(-"Millimeters")# remove empty colum

# just date column

wetdb$date<- date(wetdb$date.time)

# -- calculate the mean temp and rain -----

daily_averages <- wetdb %>%
  
  group_by(date) %>%
  summarize(
    avg_temperature = mean(Celsius, na.rm = TRUE),
    avg_wind_speed = mean(m.s, na.rm = TRUE)
  )

write.csv(daily_averages, file = "data_for_analysis/weather/dailyavg.csv",row.names = F)


mtempwind<- wetdb %>%
  group_by(wk) %>%
  summarize(mean_tem = mean(Celsius, na.rm = T), mean_wind = mean(m.s, na.rm = T))# there are NA's I don't know from where. I know now that we need to include the na.rm=T argument so there are no 

# Filter and save rows with NAs
na_rows <- wetdb[!complete.cases(wetdb), ]


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
