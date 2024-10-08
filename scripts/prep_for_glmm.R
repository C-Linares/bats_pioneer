# 
# Script name: Prep_for_glmm
# 
# Purpose of script: Combine the 2021, 2022, 2023, data sets and get them ready to go throuhg glmm_v1
# 
# Author: Carlos Linares
# 
# Date Created: 07/29/2024
# 
# Email: carlosgarcialina@u.boisestate.edu
# 
# ---------------------------
#   
#   Notes: we also should check we have the same temporal window for all variables.  
# sessionInfo() at end of script


# inputs ------------------------------------------------------------------
# -data_for_analysis/2021_kpro_raw/bats2021_kpro_v1.csv
# -data_for_analysis/2022_kpro_raw/bat2022_kpr.csv
# -data_for_analysis/2023_kpro_raw/bat2023_kpr.csv

# outputs ----------------------

# bat_combined.csv - process data no counts
# bm.csv - counts of bat calls by day from 2021 to 2023 all sites
# bm.miller.day - number of minutes of activity by day  for 2021-2023 data all sites

# this should be a database ready to analyze with the glmm_v1 script. 

# libraries

library(tidyverse)
library(beepr)


# kpro_data --------------------------------------------------------------

# load

kpro_2021_bat<- read.csv(file = 'data_for_analysis/2021_kpro_raw/bats2021_kpro_v1.csv',header = T)
kpro_2022_bat<-read.csv(file = 'data_for_analysis/2022_kpro_raw/bat2022_kpr.csv',header = T)
kpro_2023_bat<-read.csv(file = 'data_for_analysis/2023_kpro_raw/bat2023_kpr.csv',header = T)



bat_combined <- bind_rows(kpro_2021_bat, kpro_2022_bat, kpro_2023_bat)
summary(bat_combined) # check the output

keep<- c(".id","INDIR","OUTDIR","FOLDER","IN.FILE","DURATION","DATE","TIME","HOUR","AUTO.ID.", "PULSES") # cols to keep

bat_combined <- bat_combined %>% select(all_of(keep)) # keeps variables of interest


# site 

bat_combined$site<-str_extract(bat_combined$OUTDIR, "[A-Za-z]{3,4}\\d{2}")

unique(bat_combined$site) # site labels

# here we fix some problems with the sites misspellings

bat_combined$site = ifelse(bat_combined$site %in% "Iron02","iron02", bat_combined$site)
bat_combined$site = ifelse(bat_combined$site %in% "viz01","vizc01", bat_combined$site)
bat_combined$site = ifelse(bat_combined$site %in% "viz02","vizc02", bat_combined$site)
bat_combined$site = ifelse(bat_combined$site %in% "viz03","vizc03", bat_combined$site)
bat_combined$site = ifelse(bat_combined$site %in% "viz04","vizc04", bat_combined$site)

unique(bat_combined$site) # site labels


# date 

bat_combined$DATE<-lubridate::ymd(bat_combined$DATE)
sum(is.na(bat_combined$DATE)) # check for NAs. 

# noche/night
# this variable assigns the same date to calls that belong to a single night. If a call happened in the July 6 at 3 am it assigns it to July 5 

bat_combined$noche <-
  if_else(bat_combined$HOUR < 9, # if it is less than 9 put the date of the previous day
          true =  (date(bat_combined$DATE) - ddays(1)),
          false = date(bat_combined$DATE))



# date time col. 
datetime<-paste(bat_combined$DATE, bat_combined$TIME)#merge date and time # this time is wrong and needs to be extracted from the file itself. ############################!!!!!!!!!!!!!!!!
datetime.parse<-lubridate::ymd_hms(datetime) # parse as date time
bat_combined$date_time<-datetime.parse # add to data. 
sum(is.na(bat_combined$date_time)) # check for NAs. 

#year

bat_combined$yr<-year(bat_combined$DATE)
unique(bat_combined$yr) #we check the three years are present.


# treatment column 

litsites<-c("iron01","iron03","iron05","long01","long03")


bat_combined$treatmt<-ifelse(bat_combined$site %in% litsites , "lit", "dark") # this makes a treatment variable.

bat_combined$trmt_bin<- ifelse(bat_combined$treatmt== "lit", 1, 0)


bat_combined$jday<-lubridate::yday(bat_combined$noche) # julian day


summary(bat_combined)




# effort ------------------------------------------------------------------


effort_days <- bat_combined %>%
  group_by(site, yr) %>%
  summarise(
    stard = min(noche),
    endd = max(noche),
    eff.days = as.numeric(difftime(max(noche), min(noche), units = "days"))
  )

effort_hrs <- bat_combined %>%
  group_by(site, noche, jday, yr) %>%
  summarise(stard = min(date_time), endd = max(date_time)) %>%
  mutate(eff.hrs = time_length(endd - stard, unit = "hours"))

# merge effort with bat combined 

bat_combined<- left_join(bat_combined, effort_hrs, by=c("site", "jday", "yr", "noche"))


keep<- c("AUTO.ID.", "PULSES", "site","noche","date_time", "yr","treatmt","trmt_bin","jday","eff.hrs") # cols to keep

bat_combined <- bat_combined %>% select(all_of(keep))




summary(bat_combined)



# count matrix ------------------------------------------------------------
# this is a matrix where we create a n column that tells us how many calls for each bat are there.
# daily counts

bm <- bat_combined %>% # 
  group_by(noche, AUTO.ID., site,yr, trmt_bin, treatmt) %>% 
  summarise(n = n(), .groups = 'drop') 

summary(bm)

# miller matrix -----------------------------------------------------------



# minutes activity 
# in here we calculate the minutes of activity insipired by Miller 2001 paper. 

bat_combined$rmins<-round(bat_combined$date_time, units="mins") #rounds to the nearest min

bm.miller<-bat_combined %>% #min of activity 
  group_by(site, AUTO.ID., noche, rmins) %>% 
  summarize(activity_min= n()) %>%  #calculate the num of min activity
  ungroup()


bm.miller.day <- bm.miller %>% # number of minutes active  by night. 
  group_by(site, noche, AUTO.ID.) %>%
  summarize(activity_min = sum(activity_min))

summary(bm.miller.day)


# outputs -----------------------------------------------------------------


# dir.create("data_for_analysis/prep_for_glm", showWarnings = FALSE) # just run if the dir is abscent

write.csv(bat_combined, file = 'data_for_analysis/prep_for_glmm/bat_combined.csv', row.names = F) # raw combine data 
write.csv(bm, file = 'data_for_analysis/prep_for_glmm/bm.csv', row.names = F) #daily counts
write.csv(bm.miller.day, file = "data_for_analysis/prep_for_glmm/bm.miller.day.csv") # miller Ai index data



# Create a README file with information about the script
readme_content <- "Carlos Linares 8/01/2024 
This directory contains the bat_combined.csv file which was created using the script prep_for_glmm.R combines bat species call abundance data. This script merges 2021-23 data that was previously scanned with Kaleidoscope pro

bat_combined.csv - process data no counts
bm.csv - counts of bat calls by day from 2021 to 2023 all sites
bm.miller.day - number of minutes of activity by day  for 2021-2023 data all sites"

# Write the README content to a file
writeLines(readme_content, "data_for_analysis/prep_for_glmm/README.txt")







# plots -------------------------------------------------------------------

#-------------------------- calls by site and year


bat_summary <- bat_combined %>%
  filter(!AUTO.ID. %in% c("Noise", "NoID")) %>%  # Filter out Noise and NoID tags
  group_by(site, yr, AUTO.ID.) %>%
  summarise(count = n()) %>%
  ungroup()


ggplot(bat_summary, aes(x = yr, y = count, fill = AUTO.ID.)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ site, scales = "free_y") +
  labs(title = "species by sites and years",
       x = "Year",
       y = "Count",
       fill = "Tag (AUTO.ID.)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#-------------------------- calls by jday and treatment

filtered_bat_combined <- bat_combined %>%
  filter(!AUTO.ID. %in% c("Noise", "EUDMAC"))

# Summarize the number of pulses per Julian day and treatment
summary_data <- filtered_bat_combined %>%
  group_by(jday, treatmt, yr) %>%
  summarise(count = n()) %>%
  ungroup()

# Create the plot
ggplot(summary_data, aes(x = jday, y = count, col = treatmt)) +
  geom_point() +
  facet_wrap(~ yr+ treatmt , scales = "free_y") +  labs(title = "Call activity by Julian Day and Treatment",
       x = "Julian Day",
       y = "Number of calls") +
  geom_vline(xintercept = 180, linetype = "dashed", color = "red") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#-------------------------- species calls by Julian day and treatment

summary_data <- filtered_bat_combined %>%
  group_by(jday, treatmt, AUTO.ID.) %>%
  summarise(count = n()) %>%
  ungroup()


# Create the plot
ggplot(summary_data, aes(x = jday, y = count, col = treatmt)) +
  geom_point() +
  facet_wrap( ~  AUTO.ID. + treatmt, scales = "free_y") +  labs(title = "Call activity by Julian Day and Treatment", x = "Julian Day", y = "call counts") +
  geom_vline(xintercept = 180,
             linetype = "dashed",
             color = "red") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))





summary_data <- filtered_bat_combined %>%
  group_by(jday, treatmt, AUTO.ID., yr) %>%
  summarise(count = n()) %>%
  ungroup()

# Create the plot
ggplot(summary_data, aes(x = jday, y = count, col = treatmt)) +
  geom_point() +
  facet_wrap(~ AUTO.ID.+ yr , scales = "free_y") +  labs(title = "Call activity by Sp and year",
                                                               x = "Julian Day",
                                                               y = "Number of Pulses") +
  geom_vline(xintercept = 180, linetype = "dashed", color = "red") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


# pulses by treatment 

ggplot(bat_combined, aes(x = treatmt, y = scale(PULSES), fill = treatmt)) +
  geom_boxplot() +
  labs(title = "Distribution of Bat Pulses by Treatment",
       x = "Treatment",
       y = "Number of Pulses") +
  theme_minimal()



# Session info ------------------------------------------------------------



# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 22631)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: America/Denver
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] beepr_2.0       lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5    
# [8] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.5       compiler_4.4.1     tidyselect_1.2.1   scales_1.3.0       R6_2.5.1           generics_0.1.3    
# [7] sjlabelled_1.2.0   knitr_1.46         sjPlot_2.8.16      insight_0.20.2     munsell_0.5.1      pillar_1.9.0      
# [13] tzdb_0.4.0         sjstats_0.19.0     rlang_1.1.3        utf8_1.2.4         stringi_1.8.4      performance_0.11.0
# [19] xfun_0.44          audio_0.1-11       ggeffects_1.6.0    datawizard_0.12.1  timechange_0.3.0   cli_3.6.2         
# [25] withr_3.0.0        magrittr_2.0.3     grid_4.4.1         rstudioapi_0.16.0  hms_1.1.3          lifecycle_1.0.4   
# [31] sjmisc_2.8.10      vctrs_0.6.5        glue_1.7.0         fansi_1.0.6        colorspace_2.1-0   tools_4.4.1       
# [37] pkgconfig_2.0.3   
# --------------------------- trash ----------------
