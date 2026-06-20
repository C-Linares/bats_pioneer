
# =======================================================================
# Script Title:    rob_spkr_prep_v2.R
# Project:         robomoth and speaker data analysis
# =======================================================================

# Description:
# this script helps prep the data for analysis, standardize, add predictors, etc. 
# We updated the robomoth and speaker data to add an amplitude measurement. 
# 

# Author: Carlos Linares
# Date:   2026-17-06
# Contact: carlosgarcialina@u.boisestate.edu

# =======================================================================
# Inputs:
# - bat_robomoth.csv (home/r/bats_pioneer/data_for_analysis/robomoth_build/bat_robomoth_amp.csv)
# - bat_speakr.csv (home/r/bats_pioneer/data_for_analysis/robomoth_build/bat_speakr_amp.csv)

#  Outputs:
#  - ready to analyze data base ("data_for_analysis/rob_spkr_prep/rob_db_v2.csv") this is tentative. 
#  
#  
# =======================================================================

#  Load Libraries ------------------------------------------------------


# libraries 

if(!require("pacman")) install.packages("pacman")

p_load(
  "tidyverse",
  "emmeans",
  "ggplot2",
  "ggeffects",
  "sjPlot",
  "patchwork",
  "marginaleffects",
  "glmmTMB",
  "performance",
  "viridis",
  "corrplot",
  "DHARMa",
  "janitor"
)


# functions ---------------------------------------------------------------


  
run_correlation_analysis <- function(data, plot = TRUE) {
  numeric_cols <- sapply(data, is.numeric)
  cor_data <- data[, numeric_cols]
  zero_variance_cols <- sapply(cor_data, function(x) sd(x, na.rm = TRUE) == 0)
  cor_data_filtered <- cor_data[, !zero_variance_cols]
  
  cor_matrix <- cor(cor_data_filtered, use = "pairwise.complete.obs")
  
  if (plot) {
    corrplot(cor_matrix, order = 'AOE')
    corrplot(cor_matrix, order = 'AOE', addCoef.col = "black", 
             tl.pos = "lt", method = "color", cl.cex =.5, 
             tl.cex = 1, number.cex = .5)
  }
  
  return(cor_matrix)
}


# funtion to scale. 
scale_by_2sd_tidy <- function(data, variables_to_scale) {
  # Keep only variables that exist and are numeric
  valid_vars <- variables_to_scale[variables_to_scale %in% names(data) & sapply(data[variables_to_scale], is.numeric)]
  
  if (length(valid_vars) == 0) {
    warning("No valid numeric variables found to scale.")
    return(data)
  }
  
  data <- data %>%
    mutate(across(all_of(valid_vars),
                  ~ (. - mean(., na.rm = TRUE)) / (2 * sd(., na.rm = TRUE)),
                  .names = "{.col}_s"))
  
  return(data)
}



# load data ---------------------------------------------------------------


robomoth <- read_csv("data_for_analysis/robomoth_build/bat_robomoth_amp.csv") %>% 
  clean_names()

str(robomoth)
colSums(is.na(robomoth)) # there are many NAs in the max-dbfs and avg_dbfs becaue we have not filter out 2021 data 

# remove 2021 data 

# filter out 2021 from rob_db because we don't have speaker data for that year

robomoth <- robomoth %>% # we will do this in the prep script 
  filter(year != 2021)

colSums(is.na(robomoth)) # there's still 1661 NAs in the amplitude. some files are missing. I need to figure out why. we will proceed as this for now. 


# predictors --------------------------------------------------------------

# Craters weather (night)
crmo.wet.night <- read_csv("data_for_analysis/weather/craters_weater/craters_night.csv") %>% 
  clean_names() # load craters night weather


# Moon
moon.int <- read_csv('data_for_analysis/moon_pred/moon.int.csv') %>%
  clean_names()
summary(moon.int) # check the structure of the moon data

# convert date times from UTC to Denver/America
moon.int$denver_time <- with_tz(moon.int$date, tzone = "America/Denver") # Convert to Denver time zone
attr(moon.int$denver_time, "tzone") # check the timezone is correct

# Step 1: Filter data where the moon is above the horizon
moon_filtered <- moon.int %>%
  filter(moon_alt_degrees > 0)

# Step 2: Create a new 'noche' variable (just the date part of the time stamp)
# but also makes nights any time stamps that are less than 9 am.
moon_filtered <- moon_filtered %>%
  mutate(
    hour = hour(denver_time),  # Extract hour from datetime
    noche = if_else(hour < 9,
                    true = as_date(denver_time) - days(1),
                    false = as_date(denver_time))
  )

# Step 3: Group by 'noche', then summarize the average values
moon_daily_avg <- moon_filtered %>%
  group_by(noche) %>%
  summarise(
    avg_moonlight = mean(moonlight_model, na.rm = TRUE),
    avg_twilight = mean(twilight_model, na.rm = TRUE),
    avg_illumination = mean(illumination, na.rm = TRUE),
    n_obs = n()  # optional: number of observations per night
  )

# insects 
# I am not using c_bugs anymore, I switched to ins_bm that is just lepidoptera
# c_bugs <- read_csv("data_for_analysis/insect_wranglin/c_bugs.csv") %>%  # load insect data
#   clean_names() %>%
#   rename(yr = yrs) # safe rename
# # add treatment. 
# 
# litsites<-c("iron01","iron03","iron05","long01","long03")
# 
# 
# c_bugs$treatmt<-ifelse(c_bugs$site %in% litsites , "lit", "dark") # this makes a treatment variable.

ins_bm<-read.csv("data_for_analysis/insect_wranglin/ins_bm.csv") %>% 
  clean_names() %>% 
  rename(yr = yrs) # load insect biomass data and rename year column

# calculate mean by yr and site. I will use this to substitute the NA valeus from tha appeared when merging the bat data. 
# 
# c_bugs_mean <- c_bugs %>%
#   group_by(yr, site) %>%
#   summarise(
#     t_insect = mean(t_insect, na.rm = TRUE),
#     t_lepidoptera = mean(t_lepidoptera, na.rm = TRUE)
#   ) %>%
#   ungroup()

ins_bm_mean <- ins_bm %>%
  group_by(yr, site) %>%
  summarise(
    t_leps_mean = mean(t_leps, na.rm = TRUE)
  ) %>%
  ungroup()

# elevation 

elev <- read_csv("data_for_analysis/elev/elevation.csv") %>% 
  clean_names()  
# rename site
elev <- elev %>%
  rename(site = name) %>% 
  select(site, elev_max) # select only site and elevation columns

# merge -------------------------------------------------------------------


# merge weather
# Merge crmo.wet.night into filtered_bm by matching dates
rob_db <- robomoth %>%
  left_join(crmo.wet.night, by = c("noche" = "date"))
summary(rob_db)

# merge moon

rob_db <- rob_db %>%
  left_join(moon_daily_avg, by = "noche") 

summary(rob_db) #

# merge with insects
# I calculate week for insects to merge with robomot data. 
rob_db<- rob_db %>%
  mutate(wk = week(noche)) # add week and year to rob_db for merging

rob_db <- rob_db %>%
  left_join(ins_bm, by = c("site", "wk", "year"= "yr")) # merge

# check for NAs
summary(rob_db)  #lots of NAs 13 columns for the t_leps column. there are many days where there's no data for insects.  

# we remove the some of the nas from insects by adding the average year insect count by site 

rob_db <- rob_db %>%
  left_join(ins_bm_mean, 
            by = c("year"="yr", "site"),
            suffix = c("", "_mean")) %>%  # rename mean columns directly
  mutate(
    t_leps = coalesce(t_leps, t_leps_mean)
  ) %>%
  select( -t_leps_mean)

summary(rob_db) # 

# merge with elevation

rob_db<- rob_db %>%
  left_join(elev, by = "site")

summary(rob_db)


# correlation -------------------------------------------------------------

# run correlation analysis on the rob_db data frame.
# we first add 0 to hif or lof if NA
rob_db <- rob_db %>%
  mutate(
    hi_f = ifelse(is.na(hi_f), 0, hi_f),
    lo_f = ifelse(is.na(lo_f), 0, lo_f)
  )

cor.rob<-run_correlation_analysis(rob_db, plot = T)


# last check and cleanup

# remove unnecessary columns

names(rob_db)

# # there's 3 diff treatment col keep just one.
# 
# rob_db <- rob_db %>%
#   select( -treatmt, -treatmt.y) %>% 
#   rename( trmt = treatmt.x) # rename the treatment column to be consistent with the speakr data.

summary(rob_db)


# standardize -------------------------------------------------------------

#calculate Julian day from 'noche' and scale variables
rob_db<-rob_db %>% 
  mutate(
    jday = lubridate::yday(noche),  # Calculate Julian day from 'noche'
  )

# check there's no NAs for jday
sum(is.na(rob_db$jday)) # no NAs for jday)

variables_to_scale <- c(
  "avg_moonlight",
  "avg_twilight",
  "avg_illumination",
  "nit_avg_temp_c",
  "nit_avg_wspm_s",
  "t_leps",
  "jday",
  "elev_max",
  "max_dbfs",
  "avg_dbfs"
)

rob_db <- rob_db %>%
  scale_by_2sd_tidy(variables_to_scale)

summary(rob_db)

# year standardize. 

# make year between -1:1
rob_db <- rob_db %>%
  mutate(yr_s = case_when(
    year == 2022 ~ -1, # this will give NAs for the 2021 data but I am leaving those out of the analysis.
    year == 2023 ~ 1
  ))

summary(rob_db$yr_s)

# species names for graphs

species <- data.frame(
  sp = c("ANTPAL", "CORTOW", "EPTFUS", "EUDMAC", "LASCIN", "LASNOC",
         "MYOCAL", "MYOCIL", "MYOEVO", "MYOLUC", "MYOTHY", "MYOVOL",
         "MYOYUM", "PARHES"),

  species_name = c("Antrozous pallidus", "Corynorhinus townsendii", "Eptesicus fuscus", "Euderma maculatum",
                   "Lasiurus cinereus", "Lasiurus noctivagans", "Myotis californicus", "Myotis ciliolabrum",
                   "Myotis evotis", "Myotis lucifugus", "Myotis thysanodes", "Myotis volans",
                   "Myotis yumanensis", "Parastrellus hesperus")
)

# species to small caps 

species <- species %>%
  mutate(sp_s = tolower(sp)) # convert to sentence case

species <- species %>%
  mutate(
    genus = word(species_name, 1),
    species = word(species_name, 2),
    sp_label = paste0(substr(genus, 1, 1), ".", species)
  )

# 4 letter code

species <- species %>%
  mutate(
    sp_code = (tolower(substr(genus, 1, 2)) %>% paste0(substr(species, 1,2 )))
  )

# merge the rob_db and by sp_code

rob_db <- rob_db %>%
  left_join(species %>% select(sp_code, sp_label), by = c("sp" = "sp_code") )


glimpse(rob_db)



#  speakr data  -----------------------------------------------------------

# there's an error in the tags for sites on vizcaine. this needs to be fixed on robomoth_build.R script. (fixed)


speakr <- read_csv("data_for_analysis/robomoth_build/spkr_all_amp.csv") %>% 
  clean_names()

str(speakr)
colSums(is.na(speakr)) # no missing data for the amplitude columns in this data. 

unique(speakr$site) # check the unique sites in speakr

# merge -------------------------------------------------------------------

# 1. merge weather
spkr_db <- speakr %>%
  left_join(crmo.wet.night, by = c("noche" = "date"))

summary(spkr_db)

# 2. merge moon
spkr_db <- spkr_db %>%
  left_join(moon_daily_avg, by = "noche")

summary(spkr_db)

# 3. add week for insect merge
spkr_db <- spkr_db %>%
  mutate(wk = week(noche))

# 4. merge insect data
spkr_db <- spkr_db %>%
  left_join(ins_bm, by = c("site", "wk", "year" = "yr"))

summary(spkr_db)

# 5. fill missing insect values with annual site mean
spkr_db <- spkr_db %>%
  left_join(ins_bm_mean, 
            by = c("year"="yr", "site"),
            suffix = c("", "_mean")) %>%  # rename mean columns directly
  mutate(
    t_leps = coalesce(t_leps, t_leps_mean)
  ) %>%
  select( -t_leps_mean)

summary(spkr_db) # 


# 6. merge elevation
spkr_db <- spkr_db %>%
  left_join(elev, by = "site")

summary(spkr_db)

# correlation -------------------------------------------------------------

# if speaker data also has hi_f and lo_f, replace NAs with 0
if (all(c("hi_f", "lo_f") %in% names(spkr_db))) {
  spkr_db <- spkr_db %>%
    mutate(
      hi_f = ifelse(is.na(hi_f), 0, hi_f),
      lo_f = ifelse(is.na(lo_f), 0, lo_f)
    )
}

cor.spkr <- run_correlation_analysis(spkr_db, plot = TRUE)

# last check and cleanup



summary(spkr_db)
str(spkr_db)


# standardize -------------------------------------------------------------

#calculate Julian day from 'noche' and scale variables
spkr_db<-spkr_db %>% 
  mutate(
    jday = lubridate::yday(noche),  # Calculate Julian day from 'noche'
  )

# check there's no NAs in the jday column 

sum(is.na(spkr_db$jday)) # no NAs in jday column, good to go

variables_to_scale <- c(
  "avg_moonlight",
  "avg_twilight",
  "avg_illumination",
  "nit_avg_temp_c",
  "nit_avg_wspm_s",
  "t_leps",
  "jday",
  "elev_max",
  "max_dbfs",
  "avg_dbfs"
)

spkr_db <- spkr_db %>%
  scale_by_2sd_tidy(variables_to_scale)

summary(spkr_db)

# year standardize. 

# make year between -1:1
spkr_db <- spkr_db %>%
  mutate(yr_s = case_when(
    year == 2022 ~ -1,
    year == 2023 ~ 1
  ))

summary(spkr_db)
# species names for graphs


# merge the ps_label with the spkr and by sp_code

spkr_db <- spkr_db %>%
  left_join(species %>% select(sp_code, sp_label), by = c("sp" = "sp_code") )


glimpse(spkr_db)
summary(spkr_db)

#---------------------------------------------------------------------------------------------
# add zeros 
#---------------------------------------------------------------------------------------------
## now we are going to clean up the data to make it ready for analysis. 
## start with the robomoth data. 


# correct treatment binomial form 0,1 to -1 = dark, 1 = lit

rob_db <- rob_db %>%
  mutate(
    trmt_bin = if_else(treatmt == "dark", -1, 1)
  )

glimpse(rob_db)

spkr_db <- spkr_db %>%
  mutate(
    trmt_bin = if_else(treatmt == "dark", -1, 1)
  )

# filter calls to faint calls below -20 dB using the average dbfs column and drop no-species ID.
drop_ids<- c("hif", "mysp", "lof", "hilo")

rob_db_faint <- rob_db %>%
  filter(!is.na(avg_dbfs), avg_dbfs <= -20,!sp %in% drop_ids) # this removes dbfs NA and -20 dbfs observations. 

spkr_db_faint <- spkr_db %>% 
  filter(!is.na(avg_dbfs), avg_dbfs <= -20, !sp %in% drop_ids) # same filter as for robomoths. 


# now we do some kind of summary by observation site, noche, sp 
# the currently sampling unit I will use is site night and species. We also have to drop certain ID


rob_faint_counts <- rob_db_faint %>%
  count(
    site,
    noche,
    year,
    trmt_bin,
    sp,
    name = "n_calls" # this is the response variable 
  )

spkr_faint_counts <- spkr_db_faint %>% 
  count(
    site,
    noche,
    year,
    trmt_bin,
    sp,
    name = "n_calls"
  )

# I might have to add zeros, because as it is the data just represents hits. If I don't add zeros the questions might change to "when bats are detected, is there more or less calls/passes in lit sites"
# to add zeros we have to take in consideration the effort. we will do that using the amp all files created in the robomoth_build_v2.R script. 

#amp data We load amp data to build the effort table and add zeros to the data.

file_length_sec <- 15
drop_ids <- c("hif", "mysp", "lof", "hilo", "noise", "noid")
litsites <- c("iron01", "iron03", "iron05", "long01", "long03")

amp_robo_all<-read.csv(file = "data_for_analysis/robomoth_build/amp_only/amp_robo_all.csv") %>% 
  clean_names()
glimpse(amp_robo_all)

amp_spkr_all<-read.csv(file = "data_for_analysis/robomoth_build/amp_only/amp_spkr_all.csv") %>% 
  clean_names()
glimpse(amp_spkr_all)
# correct site names 

unique(amp_robo_all$site) # some vizc sites need to be relabelled. 

site_key <- c(
  "viz01" = "vizc01",
  "viz02" = "vizc02",
  "viz03" = "vizc03",
  "viz04" = "vizc04"
)

amp_robo_all <- amp_robo_all %>%
  mutate(
    site = recode(site, !!!site_key)
  )

amp_spkr_all <- amp_spkr_all %>%
  mutate(
    site = recode(site, !!!site_key)
  )

unique(amp_robo_all$site) # check that the site names are correct now.)
unique(amp_spkr_all$site) # check that the site names are correct now.)
# Build robomoth effort table from amplitude files

amp_robo_effort_raw <- amp_robo_all %>%
  mutate(
    site = str_to_lower(str_trim(site)),
    filename_norm = str_to_lower(str_trim(basename(filename_norm))),
    
    # Extract date and time from filename_norm
    # Works for patterns like iron01-20220804-220000.wav
    date_raw = str_extract(filename_norm, "20\\d{6}"),
    time_raw = str_extract(filename_norm, "(?<=20\\d{6}[-_])\\d{6}"),
    
    date_time = ymd_hms(
      paste(date_raw, time_raw),
      tz = "America/Denver"
    ),
    
    # Define monitoring night.
    # Files after midnight but before morning are assigned to previous night.
    noche = if_else(
      hour(date_time) < 9,
      as_date(date_time) - days(1),
      as_date(date_time)
    ),
    
    year = year(noche)
  )

# 2. Check date parsing
amp_robo_effort_raw %>%
  summarise(
    n_rows = n(),
    n_missing_date_raw = sum(is.na(date_raw)),
    n_missing_time_raw = sum(is.na(time_raw)),
    n_missing_date_time = sum(is.na(date_time)),
    n_missing_noche = sum(is.na(noche))
  )



# 3. Calculate effort
rob_effort <- amp_robo_effort_raw %>%
  distinct(site, noche, year, filename_norm) %>%
  group_by(site, noche, year) %>%
  summarise(
    n_files = n(),
    effort_sec = n_files * file_length_sec,
    effort_min = effort_sec / 60,
    effort_hours = effort_sec / 3600,
    .groups = "drop"
  ) %>%
  mutate(
    treatmt = if_else(site %in% litsites, "lit", "dark"),
    treatmt = factor(treatmt, levels = c("dark", "lit")),
    trmt_01 = if_else(treatmt == "lit", 1, 0),
    trmt_bin = if_else(treatmt == "lit", 1, -1)
  )

glimpse(rob_effort)

# 4. Count detected faint calls
rob_faint_counts <- rob_db %>%
  mutate(
    site = str_to_lower(str_trim(site)),
    noche = as.Date(noche),
    year = year(noche)
  ) %>%
  filter(
    year %in% c(2022, 2023),
    !is.na(avg_dbfs),
    avg_dbfs <= -20,
    !is.na(sp),
    !sp %in% drop_ids
  ) %>%
  count(site, noche, year, sp, name = "n_calls")

glimpse(rob_faint_counts)
summary(rob_faint_counts)


# 5. Define species set. Currently I have 10 sp but I might have to drop pahe and 
rob_species <- rob_db %>%
  filter(
    year %in% c(2022, 2023),
    !is.na(sp),
    !sp %in% drop_ids
  ) %>%
  distinct(sp) %>%
  arrange(sp) %>%
  pull(sp)

# 6 now we add zeros 
rob_faint_counts_zero <- rob_effort %>%
  tidyr::crossing(sp = rob_species) %>%
  left_join(
    rob_faint_counts,
    by = c("site", "noche", "year", "sp")
  ) %>%
  mutate(
    n_calls = replace_na(n_calls, 0)
  )


# 7. Check detected calls without effort
detected_without_effort <- rob_faint_counts %>%
  anti_join(
    rob_effort %>% distinct(site, noche, year),
    by = c("site", "noche", "year")
  )

detected_without_effort #now there are not detected calls without effort.

# 8 add covariates. we need to add them one by one as in the code before to avoid having NAs in the covariates. we take the rob_faint_counts_zero and add the covariates one by one.

# merge weather
# Merge crmo.wet.night into filtered_bm by matching dates
rob_db <- rob_faint_counts_zero %>%
  left_join(crmo.wet.night, by = c("noche" = "date"))
summary(rob_db)

# merge moon

rob_db <- rob_db %>%
  left_join(moon_daily_avg, by = "noche") 

summary(rob_db) #

# merge with insects
# I calculate week for insects to merge with robomot data. 
rob_db<- rob_db %>%
  mutate(wk = week(noche)) # add week and year to rob_db for merging

rob_db <- rob_db %>%
  left_join(ins_bm, by = c("site", "wk", "year"= "yr")) # merge

# check for NAs
summary(rob_db)  #lots of NAs rows for the t_leps column. there are many days where there's no data for insects.  

# we correct the nas from insects by adding the average year insect count by site 

rob_db <- rob_db %>%
  left_join(ins_bm_mean, 
            by = c("year"="yr", "site"),
            suffix = c("", "_mean")) %>%  # rename mean columns directly
  mutate(
    t_leps = coalesce(t_leps, t_leps_mean)
  ) %>%
  select( -t_leps_mean)

summary(rob_db) # 

# merge with elevation

rob_db<- rob_db %>%
  left_join(elev, by = "site")

# This is the file that should be exported now. It has the robomoth data filtered by -20 dBfs. It has covariates and is ready to analyze. 
# 
summary(rob_db) 



# speakr data -------------------------------------------------------------

# Build robomoth effort table from amplitude files

amp_spkr_effort_raw <- amp_spkr_all %>%
  mutate(
    site = str_to_lower(str_trim(site)),
    filename_norm = str_to_lower(str_trim(basename(filename_norm))),
    
    # Extract date and time from filename_norm
    # Works for patterns like iron01-20220804-220000.wav
    date_raw = str_extract(filename_norm, "20\\d{6}"),
    time_raw = str_extract(filename_norm, "(?<=20\\d{6}[-_])\\d{6}"),
    
    date_time = ymd_hms(
      paste(date_raw, time_raw),
      tz = "America/Denver"
    ),
    
    # Define monitoring night.
    # Files after midnight but before morning are assigned to previous night.
    noche = if_else(
      hour(date_time) < 9,
      as_date(date_time) - days(1),
      as_date(date_time)
    ),
    
    year = year(noche)
  )

# 2. Check date parsing
amp_spkr_effort_raw %>%
  summarise(
    n_rows = n(),
    n_missing_date_raw = sum(is.na(date_raw)),
    n_missing_time_raw = sum(is.na(time_raw)),
    n_missing_date_time = sum(is.na(date_time)),
    n_missing_noche = sum(is.na(noche))
  )



# 3. Calculate effort
spkr_effort <- amp_spkr_effort_raw %>%
  distinct(site, noche, year, filename_norm) %>%
  group_by(site, noche, year) %>%
  summarise(
    n_files = n(),
    effort_sec = n_files * file_length_sec,
    effort_min = effort_sec / 60,
    effort_hours = effort_sec / 3600,
    .groups = "drop"
  ) %>%
  mutate(
    treatmt = if_else(site %in% litsites, "lit", "dark"),
    treatmt = factor(treatmt, levels = c("dark", "lit")),
    trmt_01 = if_else(treatmt == "lit", 1, 0),
    trmt_bin = if_else(treatmt == "lit", 1, -1)
  )

glimpse(spkr_effort)

# 4. Count detected faint calls
spkr_faint_counts <- spkr_db %>%
  mutate(
    site = str_to_lower(str_trim(site)),
    noche = as.Date(noche),
    year = year(noche)
  ) %>%
  filter(
    year %in% c(2022, 2023),
    !is.na(avg_dbfs),
    avg_dbfs <= -20,
    !is.na(sp),
    !sp %in% drop_ids
  ) %>%
  count(site, noche, year, sp, name = "n_calls")

glimpse(spkr_faint_counts)
summary(spkr_faint_counts)


# 5. Define species set. Currently I have 10 sp but I might have to drop pahe and 
spkr_species <- spkr_db %>%
  filter(
    year %in% c(2022, 2023),
    !is.na(sp),
    !sp %in% drop_ids
  ) %>%
  distinct(sp) %>%
  arrange(sp) %>%
  pull(sp)

# 6 now we add zeros 
spkr_faint_counts_zero <- spkr_effort %>%
  tidyr::crossing(sp = spkr_species) %>%
  left_join(
    spkr_faint_counts,
    by = c("site", "noche", "year", "sp")
  ) %>%
  mutate(
    n_calls = replace_na(n_calls, 0)
  )


# 7. Check detected calls without effort
detected_without_effort <- spkr_faint_counts %>%
  anti_join(
    spkr_effort %>% distinct(site, noche, year),
    by = c("site", "noche", "year")
  )

detected_without_effort #now there are not detected calls without effort for spkr neither.


# now we add covariates

# 1. merge weather
spkr_db <- spkr_faint_counts_zero %>%
  left_join(crmo.wet.night, by = c("noche" = "date"))

summary(spkr_db)

# 2. merge moon
spkr_db <- spkr_db %>%
  left_join(moon_daily_avg, by = "noche")

summary(spkr_db)

# 3. add week for insect merge
spkr_db <- spkr_db %>%
  mutate(wk = week(noche))

# 4. merge insect data
spkr_db <- spkr_db %>%
  left_join(ins_bm, by = c("site", "wk", "year" = "yr"))

summary(spkr_db)

# 5. fill missing insect values with annual site mean
spkr_db <- spkr_db %>%
  left_join(ins_bm_mean, 
            by = c("year"="yr", "site"),
            suffix = c("", "_mean")) %>%  # rename mean columns directly
  mutate(
    t_leps = coalesce(t_leps, t_leps_mean)
  ) %>%
  select( -t_leps_mean)

summary(spkr_db) # 


# 6. merge elevation
spkr_db <- spkr_db %>%
  left_join(elev, by = "site")

summary(spkr_db) # this is the data that should be analyzed. 



# correlation  ------------------------------------------------------------
# Now we run a correlation for both the spkr and the robomoth data. then we need to assess if theres variables that shoud not be together in the same model

cor.rob<-run_correlation_analysis(rob_db, plot = TRUE) # for robomth wind and temp are correlated. 

cor.spkr <- run_correlation_analysis(spkr_db, plot = TRUE) #same for spkr wind and temp correlated. 



# standardize -------------------------------------------------------------

#calculate Julian day from 'noche' and scale variables
rob_db<-rob_db %>% 
  mutate(
    jday = lubridate::yday(noche),  # Calculate Julian day from 'noche'
  )

# check there's no NAs for jday
sum(is.na(rob_db$jday)) # no NAs for jday)

variables_to_scale <- c(
  "avg_moonlight",
  "avg_twilight",
  "avg_illumination",
  "nit_avg_temp_c",
  "nit_avg_wspm_s",
  "t_leps",
  "jday",
  "elev_max"
)

rob_db <- rob_db %>%
  scale_by_2sd_tidy(variables_to_scale)

summary(rob_db)

# year standardize. 

# make year between -1:1
rob_db <- rob_db %>%
  mutate(yr_s = case_when(
    year == 2022 ~ -1, # this will give NAs for the 2021 data but I am leaving those out of the analysis.
    year == 2023 ~ 1
  ))

summary(rob_db$yr_s)

# species names for graphs

species <- data.frame(
  sp = c("ANTPAL", "CORTOW", "EPTFUS", "EUDMAC", "LASCIN", "LASNOC",
         "MYOCAL", "MYOCIL", "MYOEVO", "MYOLUC", "MYOTHY", "MYOVOL",
         "MYOYUM", "PARHES"),
  
  species_name = c("Antrozous pallidus", "Corynorhinus townsendii", "Eptesicus fuscus", "Euderma maculatum",
                   "Lasiurus cinereus", "Lasiurus noctivagans", "Myotis californicus", "Myotis ciliolabrum",
                   "Myotis evotis", "Myotis lucifugus", "Myotis thysanodes", "Myotis volans",
                   "Myotis yumanensis", "Parastrellus hesperus")
)

# species to small caps 

species <- species %>%
  mutate(sp_s = tolower(sp)) # convert to sentence case

species <- species %>%
  mutate(
    genus = word(species_name, 1),
    species = word(species_name, 2),
    sp_label = paste0(substr(genus, 1, 1), ".", species)
  )

# 4 letter code

species <- species %>%
  mutate(
    sp_code = (tolower(substr(genus, 1, 2)) %>% paste0(substr(species, 1,2 )))
  )

# merge the rob_db and by sp_code

rob_db <- rob_db %>%
  left_join(species %>% select(sp_code, sp_label), by = c("sp" = "sp_code") )


glimpse(rob_db)


# now we standardize spkr data. 

#calculate Julian day from 'noche' and scale variables
spkr_db<-spkr_db %>% 
  mutate(
    jday = lubridate::yday(noche),  # Calculate Julian day from 'noche'
  )

# check there's no NAs in the jday column 

sum(is.na(spkr_db$jday)) # no NAs in jday column, good to go

variables_to_scale <- c(
  "avg_moonlight",
  "avg_twilight",
  "avg_illumination",
  "nit_avg_temp_c",
  "nit_avg_wspm_s",
  "t_leps",
  "jday",
  "elev_max"
)

spkr_db <- spkr_db %>%
  scale_by_2sd_tidy(variables_to_scale)

summary(spkr_db)

# year standardize. 

# make year between -1:1
spkr_db <- spkr_db %>%
  mutate(yr_s = case_when(
    year == 2022 ~ -1,
    year == 2023 ~ 1
  ))

summary(spkr_db)
# species names for graphs


# merge the ps_label with the spkr and by sp_code

spkr_db <- spkr_db %>%
  left_join(species %>% select(sp_code, sp_label), by = c("sp" = "sp_code") )


glimpse(spkr_db)
summary(spkr_db)

# outputs ----------------------------------------------------------------



# dir.create("data_for_analysis/prep_for_glm", showWarnings = FALSE) # just run if the dir is abscent




# Create a README file with information about the script
readme_content <- "Carlos Linares 6/17/2026 
this folde contains the database for robomoth and speaker with predictors ready to run models with. 
we updated the data and filterd calls to faint calls below -20 dB. we also added zeros for nights where no calls were detected. the data is filtered to include only 2022 and 2023 data because we don't have speaker data for 2021.
the columns are as follow:

noche: date of the night
site: site of the observation
year: year of the observation for robomoth and speaker data we have 2022 and 2023 data only. 
effort_sec: effort in seconds for that night and site. this is calculated from the number of files recorded that night and site multiplied by 15 seconds, which is the length of the files we are using for the analysis.
effort_min: effort in minutes for that night and site. this is calculated from the effort in seconds divided by 60.
effort_hours: effort in hours for that night and site. this is calculated from the effort_sec divided by 3600.
t_leps: total lepidoptera count summmarized as insects counts by week, site and year. when data for that week was missing we added the mean lepidoptera count for that site and year.
nit_avg: average night time temperature for that night
nit_avg_wspm: average night time wind speed for that night
treatment : treatment of the site, either lit or dark.
avg_illuminatio: both twilight and moon illumination
avg_moonlight: moonlight model prediction for that night. from the moonlit package (Kayba)



"
# Write the README content to a file
writeLines(readme_content, "data_for_analysis/rob_spkr_prep/README_v2.txt")


# write data
write.csv(rob_db, file = 'data_for_analysis/rob_spkr_prep/rob_db_v2.csv', row.names = F) 
write.csv(spkr_db, file = 'data_for_analysis/rob_spkr_prep/spkr_db_v2.csv', row.names = F)








