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
library(lubridate)


# kpro_data --------------------------------------------------------------

# load
# 
load("working_env/prep_for_glm.RData")




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

sum(is.na(bat_combined$noche)) # check for NAs. 


# date time col. 
bat_combined$DATE<-as.character(bat_combined$DATE)
bat_combined$TIME<-as.character(bat_combined$TIME)

datetime<-paste(bat_combined$DATE, bat_combined$TIME)#merge date and time 

bat_combined$datetime<- datetime

bat_combined <- bat_combined %>%
  mutate(
    datetime = ymd_hms(datetime, tz = "America/Denver"), # Parse as POSIXct
    # Add 10 seconds to midnight times
    datetime = if_else(
      format(datetime, "%H:%M:%S") == "00:00:00", 
      datetime + seconds(10), 
      datetime
    )
  )

# View result

midnight_times <- bat_combined$datetime %>%
  as_tibble() %>%
  filter(hour(value) == 0 & minute(value) == 0 & second(value) == 10)

# View the rows that match
midnight_rows


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
  summarise(stard = min(datetime), endd = max(datetime)) %>%
  mutate(eff.hrs = time_length(endd - stard, unit = "hours"))

# merge effort with bat combined 

bat_combined<- left_join(bat_combined, effort_hrs, by=c("site", "jday", "yr", "noche"))


keep<- c("AUTO.ID.", "PULSES", "site","noche","datetime", "yr","treatmt","trmt_bin","jday","eff.hrs") # cols to keep

bat_combined <- bat_combined %>% select(all_of(keep))

bat_combined <- bat_combined %>% rename(sp = AUTO.ID.)# change the auto.id to sp 



summary(bat_combined)



# count matrix ------------------------------------------------------------
# this is a matrix where we create a n column that tells us how many calls for each bat are there.
# daily counts

bm <- bat_combined %>% # 
  group_by(noche, sp, site,yr, trmt_bin, treatmt) %>% 
  summarise(n = n(), .groups = 'drop') 

summary(bm)
head(bm)

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


bm.miller<- bat_combined %>%
  # Extract the relevant columns and round to minute level
  mutate(rmins = round(date_time, units = "mins")) %>%
  # Remove duplicate entries for the same site, date, and minute
  distinct(site, sp, noche, rmins, .keep_all = TRUE) %>%
  # Group by site, date, and minute
  group_by(site,sp, noche, rmins) %>%
  # Summarize to count the number of unique minutes
  summarize(activity_min = n_distinct(rmins), .groups = 'drop')

bm.miller.day <- bm.miller %>% # number of minutes active  by night. 
  group_by(site, noche, sp) %>%
  summarize(activity_min = sum(activity_min))

head(bm.miller.day)
summary(bm.miller.day)



# correcting for abundance ------------------------------------------------

# This code maps sites into defined pairs, groups data by night, site pair, and species, calculates normalized bat activity by dividing experimental activity (lit treatment) by the mean control activity (dark treatment), and filters out noise and unidentified species. It also adds a Julian day column to the resulting dataset.

# Define the pairs
site_pairs <- list(
  c("long01", "long02"),
  c("long03", "long04"),
  c("iron01", "iron02"),
  c("iron03", "iron04"),
  c("iron05", "iron06")
)

# Add a column to map pairs
bm <- bm %>%
  mutate(
    pair_group = case_when(
      site %in% c("long01", "long02") ~ "long01:long02",
      site %in% c("long03", "long04") ~ "long03:long04",
      site %in% c("long05") ~ "long05",
      site %in% c("iron01", "iron02") ~ "iron01:iron02",
      site %in% c("iron03", "iron04") ~ "iron03:iron04",
      site %in% c("iron05", "iron06") ~ "iron05:iron06",
      TRUE ~ NA_character_
    )
  )

# Normalize activity
normalized_bm <- bm %>%
  group_by(noche, pair_group, sp) %>%  # Add species (sp) to the grouping
  summarize(
    control_mean = mean(n[treatmt == "dark"], na.rm = TRUE),  # Average control activity per species
    experimental_activity = sum(n[treatmt == "lit"], na.rm = TRUE),  # Total experimental activity per species
    .groups = "drop"
  ) %>%
  mutate(
    normalized_activity = experimental_activity / control_mean  # Calculate normalized activity
  )

normalized_bm <- normalized_bm[!(normalized_bm$sp %in% c("Noise","NoID")), ]# filter out Noise and NoID rows. 
normalized_bm$jday<-yday(normalized_bm$noche)

# View the result
normalized_bm
summary(normalized_bm)


# outputs -----------------------------------------------------------------


# dir.create("data_for_analysis/prep_for_glm", showWarnings = FALSE) # just run if the dir is abscent

write.csv(bat_combined, file = 'data_for_analysis/prep_for_glmm/bat_combined.csv', row.names = F) # raw combine data 
write.csv(bm, file = 'data_for_analysis/prep_for_glmm/bm.csv', row.names = F) #daily counts
write.csv(bm.miller.day, file = "data_for_analysis/prep_for_glmm/bm.miller.day.csv") # miller Ai index data
write.csv(normalized_bm, file = "data_for_analysis/prep_for_glmm/normalized_bm.csv")


# Create a README file with information about the script
readme_content <- "Carlos Linares 8/01/2024, 12/12/2024 
This directory contains the bat_combined.csv file which was created using the script prep_for_glmm.R combines bat species call abundance data. This script merges 2021-23 data that was previously scanned with Kaleidoscope pro

bat_combined.csv - process data no counts (Update: 12/2/2024 some modifications to date time column.) 
bm.csv - counts of bat calls by day from 2021 to 2023 all sites
bm.miller.day - number of minutes of activity by day  for 2021-2023 data all sites (last update 9/23/2024)
normalized_bm.csv - normalized bat activity by dividing experimental activity (lit treatment) by the mean control activity (dark treatment) by species and noche(date) "
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


# plots_normalized_activity  ----------------------------------------------

p7<-ggplot(normalized_bm, aes(x=jday, y=normalized_activity, color=pair_group))+
  geom_point()
p7

p7 <- ggplot(normalized_bm, aes(x = jday, y = normalized_activity, color = pair_group)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(color = pair_group)) +  # Adds a linear trend line
  facet_wrap(~ sp, scales = "free") +
  labs(
    title = "Normalized Bat Activity by Species with Trend Lines",
    x = "Julian Day (jday)",
    y = "Normalized Activity",
    color = "Site Pair Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_color_viridis_d()

# Print the plot


# save --------------------------------------------------------------------

save.image("working_env/prep_for_glm.RData")

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



midnight_strings <- datetime[grepl("^\\d{4}-\\d{2}-\\d{2} 00:00:00$", datetime)] # make sure the midnight strings are there. 

t<-ymd_hms(head(midnight_strings), tz = "America/Denver")

formatted_times <- format(t, "%Y-%m-%d %H:%M:%S %Z")
print(formatted_times)\\\



bat_combined$date_time<-ymd_hms(datetime) # add to data. 
sum(is.na(bat_combined$date_time)) # check for NAs. 