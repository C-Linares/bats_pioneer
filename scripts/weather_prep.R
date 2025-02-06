
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



# inputs ------------------------------------------------------------------

# data source: Data was obtainde from https://mesowest.utah.edu/ statition KSUN in Idaho. 
# Craters of the moon data: obtain directly from Craters. 


# outputs -----------------------------------------------------------------

# data with weather data summarized. 


#----------- load data ---------------

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

# extnullfile()# extract the last digit of the cloud code. This correspond to the cloud magnitude. 
cls.code<-str_extract(cls$code, "\\d$")

write.csv(daily_averages, file = "data_for_analysis/weather/dailyavg.csv",row.names = F)
write.csv(daily_averages, file = "data_for_analysis/weather/nigh_averages.csv",row.names = F)




#write data

write.csv(wetdb,file = "data_for_analysis/weather/rawweather.csv",row.names = F)




# crater data -------------------------------------------------------------

# here we load the creates of the moon data. There are two types of data sets one called export.csv and another called Creaters_of_Moon_2023. they both have complementary data.
# 

cmon<-read_csv(file = "data_for_analysis/weather/craters_weater/CRMO_met_2021-2023.csv", skip = 9, col_names = T)

# check for NAs
sum(is.na(cmon)) # none

# check for duplicate
sum(duplicated(cmon)) # none

# change col names
cmon <- setnames(
  cmon,
  old = c(
    "DATE_TIME",
    "CRMO-VC_SWS_M_S",
    "CRMO-VC_VWS_M_S",
    "CRMO-VC_VWD_DEG",
    "CRMO-VC_TMP_DEGC",
    "CRMO-VC_SOL_W_M2"
  ),
  new = c(
    "date.time",
    "wnspd_ms",
    "vwndspd_ms",
    "wnddir_deg",
    "temp_degc",
    "sol_wm2"
  )
) # change col name to date

# date and time 

cmon$date_time<- mdy_hm(cmon$date.time, tz = "America/Denver") # seems 3 failed to parse but it is ok because they are in march.

# Extract the date part
cmon$date <- as.Date(cmon$date_time)

# Extract the time part

cmon$hour <- hour(cmon$date_time)


# filter months 5-9 may to sept

cromo_wtr <- cmon %>%
  filter(month(date) %in% 5:9)

summary(cromo_wtr)

# summarize data by night. 

cromo_night <- cromo_wtr %>%
  filter(hour >= 18 | hour < 6) %>%
  group_by(date) %>%
  summarize(
    nit_avg_tempC = mean(temp_degc, na.rm = TRUE),
    nit_avg_wspm.s = mean(wnspd_ms, na.rm = TRUE)
  )

summary(cromo_night)

# summarize data by day

cromo_day <- cromo_wtr %>%
  filter(hour >= 6 & hour < 18) %>%
  group_by(date) %>%
  summarize(
    day_avg_tempC = mean(temp_degc, na.rm = TRUE),
    day_avg_wspm.s = mean(wnspd_ms, na.rm = TRUE)
  )

# Write data --------------------------------------------------------------

# night and day summaries
write.csv(cromo_night, file = "data_for_analysis/weather/craters_weater/craters_night.csv", row.names = F)
write.csv(cromo_day, file = "data_for_analysis/weather/craters_weater/craters_day.csv", row.names = F)
# hourly data raw
write.csv(cromo_wtr, file = "data_for_analysis/weather/craters_weater/craters_wtr.csv", row.names = F)

# Create a README file with information about the script
readme_content <- "Carlos Linares, 2/5/2025
the firle craters_wtr.csv contains the weather data from the Craters of the Moon National Monument and Preserve. The data was collected from 2021 to 2023. The data was filtered to include only the months of May to September. The data includes the following columns: date.time, wnspd_ms, vwndspd_ms, wnddir_deg, temp_degc, sol_wm2, date, and hour. It was produced with the weather_prep.R script.
 "
# Write the README content to a file
writeLines(readme_content, "data_for_analysis/weather/craters_weater/README.txt")



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

