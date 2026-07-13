

#=============================================================================
# Script Name:    effort.R
# Purpose:        calculate the effort for site and night for each sm3 deployeded. 
#
# Version:        v1
#
# Author:         Carlos Linares
# Created:        2025/7/11
# Contact:        carlosgarcialina@u.boisestate.edu
#
# Notes:
# 
#
# Inputs:
#   - G:\PioneerLights_2021\iron01\sm3\
#   - G:\PioneerLights_2022\sm3
#   - G:\PioneerLights_2023\sm3\sm3
#
# Outputs:
#   - Visualizations of model predictions
#   - Diagnostic plots for model assumptions
# =============================================================================

# libraries 

library(tidyverse)
library(lubridate)
library(magrittr)


# paths -------------------------------------------------------------------

# Folder containing all site folders for 2021
project_root <- "G:/"



# Identify the experimental-year folders
year_folders <- list.dirs(
  path = project_root,
  recursive = FALSE,
  full.names = TRUE
) %>%
  keep(~ str_detect(
    basename(.x),
    regex("^PioneerLights_20(21|22|23)$", ignore_case = TRUE)
  ))

year_folders


# Explicitly define the three yearly folders
year_folders <- tibble(
  year = 2021:2023,
  year_path = paste0("G:/PioneerLights_", year)
)

year_folders


# find all the ultrasonic calls 

ultrasonic_files <- year_folders %>%
  mutate(
    files = map(
      year_path,
      ~ list.files(
        path = .x,
        pattern = "__1__",
        recursive = TRUE,
        full.names = TRUE,
        include.dirs = FALSE,
        ignore.case = TRUE
      )
    )
  ) %>%
  select(year, files) %>%
  unnest(files) %>%
  rename(file_path = files) %>%
  mutate(
    filename = basename(file_path),
    file_path = str_replace_all(file_path, "\\\\", "/")
  )


ultrasonic_files %>%
  count(year)



# extract the site. 


valid_sites <- c(
  "iron01", "iron02", "iron03", "iron04", "iron05", "iron06",
  "long01", "long02", "long03", "long04", "long05",
  "vizc01", "vizc02", "vizc03", "vizc04"
)

site_pattern <- paste(valid_sites, collapse = "|")

# extract the site from the file path
ultrasonic_files <- ultrasonic_files %>%
  mutate(
    file_path_lower = str_to_lower(file_path),
    
    site = str_extract(
      file_path_lower,
      paste0("(?i)(?<![a-z0-9])(", site_pattern, ")(?![a-z0-9])")
    ),
    
    site = str_to_lower(site)
  )

ultrasonic_files %>%
  count(year, site) %>%
  arrange(year, site) %>%
  print(n = Inf)


# missing site 

missing_sites <- ultrasonic_files %>%
  filter(is.na(site)) %>%
  distinct(year, file_path)

missing_sites %>%
  print(n = 50)




site_variants <- c(
  valid_sites,
  "lon03",
  "viz01", "viz02", "viz03", "viz04"
)

site_pattern <- paste(site_variants, collapse = "|")

ultrasonic_files <- ultrasonic_files %>%
  mutate(
    file_path_lower = str_to_lower(file_path),
    
    site = str_extract(
      file_path_lower,
      paste0("(?i)(?<![a-z0-9])(", site_pattern, ")(?![a-z0-9])")
    ),
    
    site = case_when(
      site == "lon03" ~ "long03",
      site == "viz01" ~ "vizc01",
      site == "viz02" ~ "vizc02",
      site == "viz03" ~ "vizc03",
      site == "viz04" ~ "vizc04",
      TRUE ~ site
    )
  )


table(ultrasonic_files$site, useNA = "ifany")



#extract the date and time


ultrasonic_files <- ultrasonic_files %>%
  mutate(
    timestamp_text = str_extract(
      filename,
      "\\d{8}[_-]\\d{6}"
    ),
    
    datetime = ymd_hms(
      str_replace(timestamp_text, "_", " "),
      tz = "America/Denver",
      quiet = TRUE
    )
  )

ultrasonic_files %>%
  select(year, site, filename, timestamp_text, datetime) %>%
  slice_head(n = 30)


ultrasonic_files %>%
  filter(is.na(datetime)) %>%
  distinct(year, filename) %>%
  print(n = 50)


# remove duplicates 

duplicate_recordings <- ultrasonic_files %>%
  filter(!is.na(site), !is.na(datetime)) %>%
  count(year, site, datetime, name = "n_copies") %>%
  filter(n_copies > 1) %>%
  arrange(desc(n_copies))

ultrasonic_files %>%
  semi_join(
    duplicate_recordings,
    by = c("year", "site", "datetime")
  ) %>%
  arrange(year, site, datetime) %>%
  select(year, site, datetime, filename, file_path)





#retain only unique files 

ultrasonic_unique <- ultrasonic_files %>%
  filter(
    !is.na(site),
    !is.na(datetime)
  ) %>%
  distinct(
    year, site, datetime,
    .keep_all = TRUE
  )

#how many dupplicate copies removed about 664k
tibble(
  files_before = nrow(ultrasonic_files),
  unique_recordings = nrow(ultrasonic_unique),
  copies_removed = nrow(ultrasonic_files) -
    nrow(ultrasonic_unique)
)


# Remove recordings made during daytime -------------------------------

# Keep an unfiltered copy for checking
ultrasonic_unique_all <- ultrasonic_unique

ultrasonic_unique <- ultrasonic_unique %>%
  mutate(
    recording_hour = hour(datetime),
    
    # Plausible nocturnal recording:
    # 18:00–23:59 or 00:00–07:59
    nocturnal_recording =
      recording_hour >= 18 | recording_hour < 8
  ) %>%
  filter(nocturnal_recording)




# check what we are removing. 

daytime_recordings <- ultrasonic_unique_all %>%
  mutate(
    recording_hour = hour(datetime)
  ) %>%
  filter(
    recording_hour >= 8,
    recording_hour < 18
  ) %>%
  arrange(year, site, datetime)

daytime_recordings %>%
  select(
    year, site, datetime,
    filename, file_path
  ) %>%
  print(n = 50)


# calculate effort per night and per year.
effort_hrs <- ultrasonic_unique %>%
  mutate(
    monitoring_night = as_date(datetime - hours(12)), # this produces the correct monitoring night. 
    jday = yday(monitoring_night)
  ) %>%
  group_by(
    site,
    year,
    monitoring_night,
    jday
  ) %>%
  summarise(
    stard = min(datetime),
    endd = max(datetime),
    n_recordings = n_distinct(datetime),
    
    eff.hrs = as.numeric(
      difftime(endd, stard, units = "hours")
    ),
    
    .groups = "drop"
  ) %>%
  arrange(year, site, monitoring_night)

effort_hrs

effort_days <- effort_hrs %>%
  group_by(site, year) %>%
  summarise(
    first_night = min(monitoring_night),
    last_night = max(monitoring_night),
    eff.days = n_distinct(monitoring_night),
    total_eff.hrs = sum(eff.hrs, na.rm = TRUE),
    total_recordings = sum(n_recordings),
    .groups = "drop"
  ) %>%
  arrange(year, site)




effort_hrs %>%
  summarise(
    sites = n_distinct(site),
    nights = n_distinct(monitoring_night),
    recordings = sum(n_recordings),
    total_effort = sum(eff.hrs),
    .by = year
  )


# I think its ready to export. all times are less than 14 hours for the time of the year. 



readme_content<- "this folder contains the effort calculations taken as reference the recordings from the sm3 recorders stored in file. 
time and dates where stracted from these to be able to estimate nightly effort and yearly effort. 

A tibble: 3 × 5
   year sites nights recordings total_effort (hrs)
  <int> <int>  <int>      <int>        <dbl>
1  2021    15     64     696061        6629.
2  2022    15     65     785079        5668.
3  2023    15     77     702420        7891.

"

# Write the README content to a file
writeLines(readme_content, "data_for_analysis/effort/README.txt")


# write effort files 

write.csv(effort_days, "data_for_analysis/effort/effort_days.csv", row.names = FALSE)
write.csv(effort_hrs, "data_for_analysis/effort/effort_hrs.csv", row.names = FALSE)










# junk --------------------------------------------------------------------


# now we list all the files 

ultrasonic_files <- map_dfr(
  ultrasonic_dirs,
  ~ tibble(
    file_path = list.files(
      path = .x,
      recursive = TRUE,
      full.names = TRUE,
      include.dirs = FALSE
    )
  )
) %>%
  mutate(
    filename = basename(file_path)
  ) %>%
  # Keep only ultrasonic-channel recordings
  filter(str_detect(filename, fixed("__1__")))

glimpse(ultrasonic_files)
nrow(ultrasonic_files)


# extract the site

ultrasonic_files <- ultrasonic_files %>%
  mutate(
    site = str_match(
      file_path,
      regex("[/\\\\]([^/\\\\]+)[/\\\\]sm3[/\\\\]", ignore_case = TRUE)
    )[, 2],
    site = str_to_lower(site)
  )

ultrasonic_files %>%
  count(site, sort = TRUE)

unique(ultrasonic_files$site)


# any files failed?
ultrasonic_files %>%
  filter(is.na(site)) %>%
  select(file_path)


ultrasonic_files %>%
  distinct(filename) %>%
  slice_head(n = 30) %>%
  pull(filename)


# extract time stamps 
ultrasonic_files <- ultrasonic_files %>%
  mutate(
    date_text = str_extract(filename, "\\d{8}(?=[_-]\\d{6})"),
    time_text = str_extract(filename, "(?<=\\d{8}[_-])\\d{6}"),
    
    datetime = ymd_hms(
      paste(date_text, time_text),
      tz = "America/Denver",
      quiet = TRUE
    )
  )

# we check for NAs
ultrasonic_files %>%
  select(site, filename, date_text, time_text, datetime) %>%
  slice_head(n = 30)

sum(is.na(ultrasonic_files$datetime)) 

ultrasonic_files %>%
  filter(is.na(datetime)) %>%
  distinct(filename) %>%
  slice_head(n = 50)

duplicate_recordings <- ultrasonic_files %>%
  count(site, datetime, name = "n_copies") %>%
  filter(n_copies > 1) %>%
  arrange(desc(n_copies))

duplicate_recordings # we have duplicates but that was expected because we have the all_ultrasonic and the regular data files. 

# now I want to know where the copies are

ultrasonic_files %>%
  semi_join(
    duplicate_recordings,
    by = c("site", "datetime")
  ) %>%
  arrange(site, datetime) %>%
  select(site, datetime, filename, file_path)


# keep just one recording per site and datetime. 
ultrasonic_unique <- ultrasonic_files %>%
  filter(!is.na(site), !is.na(datetime)) %>%
  distinct(site, datetime, .keep_all = TRUE)


# noon to noon to get it ride of the issue with the night spliting into two dates 
ultrasonic_unique <- ultrasonic_unique %>%
  mutate(
    monitoring_night = as_date(datetime - hours(12)),
    year = year(monitoring_night)
  )

# se we calculated effort and we found some issues. 
#ther are nights where it says that the recording was 18  hours long. That can be as it is from sunset to sunrise. 
#I check those site and it seems that the timing was accidentally set at 12 noon for those ones. 

effort_night <- ultrasonic_unique %>%
  mutate(
    # Assign after-midnight recordings to the preceding night
    monitoring_night = as_date(datetime - hours(12)),
    year = year(monitoring_night)
  ) %>%
  group_by(site, year, monitoring_night) %>%
  summarise(
    first_recording = min(datetime),
    last_recording  = max(datetime),
    n_recordings    = n_distinct(datetime),
    eff.hrs = as.numeric(
      difftime(last_recording, first_recording, units = "hours")
    ),
    .groups = "drop"
  )

effort_year <- effort_night %>%
  group_by(site, year) %>%
  summarise(
    first_night = min(monitoring_night),
    last_night  = max(monitoring_night),
    eff.days    = n_distinct(monitoring_night),
    eff.hrs     = sum(eff.hrs, na.rm = TRUE),
    total_recordings = sum(n_recordings),
    mean_hrs_per_night = mean(eff.hrs, na.rm = TRUE),
    median_hrs_per_night = median(eff.hrs, na.rm = TRUE),
    .groups = "drop"
  )

effort_year
