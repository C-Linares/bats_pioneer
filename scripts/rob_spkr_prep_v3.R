
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
  clean_names() %>% 
  select(-any_of("time"))

str(robomoth)
colSums(is.na(robomoth)) # there are 4577 NAs in the max-dbfs and avg_dbfs because we have not filter out 2021 data. 

# remove 2021 data 

# filter out 2021 from rob_db because we don't have speaker data for that year

robomoth <- robomoth %>% # we will do this in the prep script 
  filter(year != 2021)

colSums(is.na(robomoth)) # there's still 1661 NAs in the amplitude. some files are missing. I need to figure out why. we will proceed as this for now. 

summary(robomoth) # check the structure of the robomoth data.
# predictors --------------------------------------------------------------

# Craters weather (night)
crmo.wet.night <- read_csv("data_for_analysis/weather/craters_weater/craters_night.csv") %>% 
  clean_names() # load craters night weather
glimpse(crmo.wet.night)

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

summary(rob_db) # 

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
summary(rob_db)  #lots of NAs 13 k for the t_leps column. there are many days where there's no data for insects.  

# we remove the some of the nas from insects by adding the average year insect count by site 

rob_db <- rob_db %>%
  left_join(ins_bm_mean, 
            by = c("year"="yr", "site"),
            suffix = c("", "_mean")) %>%  # rename mean columns directly
  mutate(
    t_leps = coalesce(t_leps, t_leps_mean)
  ) %>%
  select( -t_leps_mean)

summary(rob_db)

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

speakr <- read_csv("data_for_analysis/robomoth_build/spkr_all_amp.csv") %>% 
  clean_names() %>% 
  select(-any_of("time"))

str(speakr)
colSums(is.na(speakr)) # no missing data for the amplitude columns in this data. 

unique(speakr$site) # check the unique sites in speakr

summary(speakr)

# merge -------------------------------------------------------------------
# when asking for the summary the time column gives an error. I am going to try to remove it. 
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
  filter(!is.na(avg_dbfs), avg_dbfs >= -40,!sp %in% drop_ids) # this removes dbfs NA and calls  lower than -40 dbfs. 

spkr_db_faint <- spkr_db %>% 
  filter(!is.na(avg_dbfs), avg_dbfs >= -40, !sp %in% drop_ids) # same filter as for robomoths. 


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
    avg_dbfs >= -40,
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

# This is the file that should be exported now. It has the robomoth data filtered by -40 dBfs. It has covariates and is ready to analyze. 
# 
summary(rob_db) 



# speakr data -------------------------------------------------------------

# Build spkr effort table from amplitude files

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

amp_spkr_effort_raw
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
    avg_dbfs >= -40,
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
# Now we run a correlation for both the spkr and the robomoth data. then we need to assess if there are variables that should not be together in the same model

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


# Now we divide the n_calls at noche site, year, sp by tocal calls for that sp at a noche, site, but in the sm3 data. 
# for this we are going to use the bm2 file from the prep for glmm V2 script. 
# I ill be creating a new variable that is ncalls in spkr and robomoth / total calls in the background. 

bm2<- read.csv("data_for_analysis/prep_for_glmm_v2/bm2.csv") %>% 
  clean_names() 

summary(bm2)

# bm2 uses 6 letter species codes while robomoth and speakr use 4 letter code. 

bm2_background <- bm2 %>%
  mutate(
    # Make sure site names are lower case and clean
    site = str_to_lower(str_trim(site)),
    
    # Make sure noche is a Date object
    noche = as.Date(noche),
    
    # Make year name consistent with rob_db and spkr_db
    year = yr,
    
    # Standardize bm2 species codes to match rob_db and spkr_db
    sp_lure_code = case_when(
      sp == "antpal" ~ "anpa",
      sp == "cortow" ~ "coto",
      sp == "eptfus" ~ "epfu",
      sp == "eudmac" ~ "euma",
      sp == "lascin" ~ "laci",
      sp == "lasnoc" ~ "lano",
      sp == "myocal" ~ "myca",
      sp == "myocil" ~ "myci",
      sp == "myoevo" ~ "myev",
      sp == "myoluc" ~ "mylu",
      sp == "myothy" ~ "myth",
      sp == "myovol" ~ "myvo",
      sp == "myoyum" ~ "myyu",
      sp == "parhes" ~ "pahe",
      TRUE ~ sp
    )
  ) %>%
  group_by(site, noche, year, sp = sp_lure_code) %>%
  summarise(
    # This is the total background acoustic activity for that species,
    # at that site, on that night.
    n_background = sum(n, na.rm = TRUE),
    
    # Optional: keep the background monitoring effort for checking.
    # If eff_hrs is identical within each site-night-species, first() is fine.
    eff_hrs_background = first(eff_hrs),
    
    .groups = "drop"
  )

# 2. Check that the background table looks correct
glimpse(bm2_background)

bm2_background %>%
  summarise(
    n_rows = n(),
    n_species = n_distinct(sp),
    n_site_nights = n_distinct(paste(site, noche)),
    min_background_calls = min(n_background, na.rm = TRUE),
    max_background_calls = max(n_background, na.rm = TRUE)
  )

# now we merge the backround data with te robomoth data. 


rob_db <- rob_db %>%
  mutate(
    site = str_to_lower(str_trim(site)),
    noche = as.Date(noche)
  ) %>%
  left_join(
    bm2_background,
    by = c("site", "noche", "year", "sp")
  ) %>%
  mutate(
    # If n_background is missing, there was no matching row in bm2.
    # This should be investigated before modeling.
    background_missing = is.na(n_background),
    
    # Raw ratio: lure-associated calls divided by total background calls.
    # If n_background == 0, the ratio is undefined, so we use NA_real_.
    n_calls_abund_controlled = if_else(
      !is.na(n_background) & n_background > 0,
      n_calls / n_background,
      NA_real_
    ),
    
    # Optional: calls per 100 background calls.
    # This is easier to interpret in tables and figures.
    n_calls_per_100_background = 100 * n_calls_abund_controlled
  )

# 4. Check robomoth merge
rob_db %>%
  summarise(
    n_rows = n(),
    missing_background = sum(is.na(n_background)),
    zero_background = sum(n_background == 0, na.rm = TRUE),
    positive_background = sum(n_background > 0, na.rm = TRUE),
    min_ratio = min(n_calls_abund_controlled, na.rm = TRUE),
    max_ratio = max(n_calls_abund_controlled, na.rm = TRUE)
  )

# 5. Show rows where background abundance did not match
rob_db %>%
  filter(is.na(n_background)) %>%
  distinct(site, noche, year, sp) %>%
  arrange(site, noche, sp)


rob_db %>%
  mutate(
    background_status = case_when(
      is.na(n_background) ~ "missing_background",
      n_background == 0 ~ "zero_background",
      n_background > 0 ~ "positive_background"
    )
  ) %>%
  group_by(background_status) %>%
  summarise(
    n_rows = n(),
    total_lure_calls = sum(n_calls, na.rm = TRUE),
    rows_with_lure_calls = sum(n_calls > 0, na.rm = TRUE),
    .groups = "drop"
  )


# we add background data to the speakr data --------------------------------------------------------

spkr_db <- spkr_db %>%
  mutate(
    site = str_to_lower(str_trim(site)),
    noche = as.Date(noche)
  ) %>%
  left_join(
    bm2_background,
    by = c("site", "noche", "year", "sp")
  ) %>%
  mutate(
    background_missing = is.na(n_background),
    
    n_calls_abund_controlled = if_else(
      !is.na(n_background) & n_background > 0,
      n_calls / n_background,
      NA_real_
    ),
    
    n_calls_per_100_background = 100 * n_calls_abund_controlled
  )

# 7. Check speaker merge
spkr_db %>%
  summarise(
    n_rows = n(),
    missing_background = sum(is.na(n_background)),
    zero_background = sum(n_background == 0, na.rm = TRUE),
    positive_background = sum(n_background > 0, na.rm = TRUE),
    min_ratio = min(n_calls_abund_controlled, na.rm = TRUE),
    max_ratio = max(n_calls_abund_controlled, na.rm = TRUE)
  )

# 8. Show rows where background abundance did not match
spkr_db %>%
  filter(is.na(n_background)) %>%
  distinct(site, noche, year, sp) %>%
  arrange(site, noche, sp) 



# now we have to handle the zeros as denominator 

# ------------------------------------------------------------------------------
# Impute zero or missing background activity for abundance-offset models
# ------------------------------------------------------------------------------
# Goal:
# We want to use background acoustic activity as an offset:
#
#   offset(log(n_background_for_offset))
#
# However, log(0) is undefined. Therefore, rows where n_background == 0
# cannot be used directly in an offset model.
#
# Here we create an imputed denominator:
#
#   n_background_for_offset
#
# Rules:
# 1. If observed background calls are positive, keep the observed value.
# 2. If observed background calls are zero, replace with the mean positive
#    background activity for the same site, year, and species.
# 3. If background is missing, also try to replace with the same site-year-species
#    mean, but we flag these rows separately because missing background means
#    there was no matching bm2 row.
#
# Important:
# This does not change the response variable n_calls.
# It only creates a denominator/exposure variable for the model offset.
# ------------------------------------------------------------------------------

# 1. Create site-year-species average background activity from bm2_background
#    We use only positive background counts to avoid calculating an average of zero
#    when the species was never detected.
bm2_background_means <- bm2_background %>%
  group_by(site, year, sp) %>%
  summarise(
    mean_positive_background = mean(n_background[n_background > 0], na.rm = TRUE),
    median_positive_background = median(n_background[n_background > 0], na.rm = TRUE),
    n_positive_nights = sum(n_background > 0, na.rm = TRUE),
    n_zero_nights = sum(n_background == 0, na.rm = TRUE),
    n_total_nights = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # If a species was never detected at a site-year combination,
    # mean_positive_background will be NaN. Convert that to NA.
    mean_positive_background = if_else(
      is.nan(mean_positive_background),
      NA_real_,
      mean_positive_background
    ),
    median_positive_background = if_else(
      is.nan(median_positive_background),
      NA_real_,
      median_positive_background
    )
  )

# Check the imputation table
bm2_background_means %>%
  summarise(
    n_rows = n(),
    missing_mean_background = sum(is.na(mean_positive_background)),
    min_mean_background = min(mean_positive_background, na.rm = TRUE),
    max_mean_background = max(mean_positive_background, na.rm = TRUE)
  )


# robomoth data. here we apply the imputation rules to the robomoth data.

# ------------------------------------------------------------------------------
# 2. Add imputed background denominator to robomoth data
# ------------------------------------------------------------------------------

rob_db <- rob_db %>%
  # Remove previous imputation columns if this section has been run before
  select(
    -any_of(c(
      "mean_positive_background",
      "median_positive_background",
      "n_positive_nights",
      "n_zero_nights",
      "n_total_nights",
      "background_status",
      "background_was_imputed",
      "n_background_for_offset",
      "log_background_offset"
    ))
  ) %>%
  left_join(
    bm2_background_means,
    by = c("site", "year", "sp")
  ) %>%
  mutate(
    # Classify each row based on observed background activity
    background_status = case_when(
      is.na(n_background) ~ "missing_background",
      n_background == 0 ~ "zero_background",
      n_background > 0 ~ "positive_background"
    ),
    
    # Create the denominator used in the offset model
    n_background_for_offset = case_when(
      n_background > 0 ~ as.numeric(n_background),
      n_background == 0 & !is.na(mean_positive_background) ~ mean_positive_background,
      is.na(n_background) & !is.na(mean_positive_background) ~ mean_positive_background,
      TRUE ~ NA_real_
    ),
    
    # Flag rows where the background value was not observed directly
    background_was_imputed = case_when(
      n_background > 0 ~ FALSE,
      n_background == 0 & !is.na(n_background_for_offset) ~ TRUE,
      is.na(n_background) & !is.na(n_background_for_offset) ~ TRUE,
      TRUE ~ NA
    ),
    
    # Create the log offset
    log_background_offset = log(n_background_for_offset)
  )

# 3. Check how many robomoth rows were kept or imputed
rob_db %>%
  count(background_status, background_was_imputed)

rob_db %>%
  summarise(
    n_rows = n(),
    observed_positive_background = sum(background_status == "positive_background", na.rm = TRUE),
    imputed_background = sum(background_was_imputed == TRUE, na.rm = TRUE),
    still_missing_offset = sum(is.na(log_background_offset)),
    min_background_for_offset = min(n_background_for_offset, na.rm = TRUE),
    max_background_for_offset = max(n_background_for_offset, na.rm = TRUE)
  )


# applythis to the spkr data. 

# ------------------------------------------------------------------------------
# 4. Add imputed background denominator to speaker data
# ------------------------------------------------------------------------------

spkr_db <- spkr_db %>%
  select(
    -any_of(c(
      "mean_positive_background",
      "median_positive_background",
      "n_positive_nights",
      "n_zero_nights",
      "n_total_nights",
      "background_status",
      "background_was_imputed",
      "n_background_for_offset",
      "log_background_offset"
    ))
  ) %>%
  left_join(
    bm2_background_means,
    by = c("site", "year", "sp")
  ) %>%
  mutate(
    background_status = case_when(
      is.na(n_background) ~ "missing_background",
      n_background == 0 ~ "zero_background",
      n_background > 0 ~ "positive_background"
    ),
    
    n_background_for_offset = case_when(
      n_background > 0 ~ as.numeric(n_background),
      n_background == 0 & !is.na(mean_positive_background) ~ mean_positive_background,
      is.na(n_background) & !is.na(mean_positive_background) ~ mean_positive_background,
      TRUE ~ NA_real_
    ),
    
    background_was_imputed = case_when(
      n_background > 0 ~ FALSE,
      n_background == 0 & !is.na(n_background_for_offset) ~ TRUE,
      is.na(n_background) & !is.na(n_background_for_offset) ~ TRUE,
      TRUE ~ NA
    ),
    
    log_background_offset = log(n_background_for_offset)
  )

# 5. Check how many speaker rows were kept or imputed
spkr_db %>%
  count(background_status, background_was_imputed)

spkr_db %>%
  summarise(
    n_rows = n(),
    observed_positive_background = sum(background_status == "positive_background", na.rm = TRUE),
    imputed_background = sum(background_was_imputed == TRUE, na.rm = TRUE),
    still_missing_offset = sum(is.na(log_background_offset)),
    min_background_for_offset = min(n_background_for_offset, na.rm = TRUE),
    max_background_for_offset = max(n_background_for_offset, na.rm = TRUE)
  )

# now we save this data again. 


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

the rob_db_v3.csv and spkr_db_v3.csv are the final databases for robomoth and speaker data respectively. I have added the sm3 calls as background data. I will try to analyze these with the sm3 activity as an offset. 


after adding the sm3 data (bm2), we have created a variable to use as an offset in the model.
n_background:
Observed background acoustic activity from the SM3 dataset. This is the total number of calls for the same species, site, and night in the broader acoustic monitoring dataset.

eff_hrs_background:
Background acoustic monitoring effort in hours from the SM3 dataset for the corresponding site and night.

background_missing:
Logical variable indicating whether no matching background acoustic value was found in the SM3 dataset for that species-site-night combination.

n_calls_abund_controlled:
Ratio of lure-associated calls to observed background acoustic calls, calculated as n_calls / n_background when n_background > 0. Rows with zero or missing background activity are assigned NA because division by zero is undefined.

n_calls_per_100_background:
Lure-associated calls per 100 background acoustic calls. Calculated as 100 * n_calls_abund_controlled. This variable is useful for summaries and figures but is not the main count response.

mean_positive_background:
Mean positive background acoustic activity for the same site, year, and species. This was calculated only from nights where n_background > 0 and can be used to replace zero or missing background values for offset-based sensitivity models.

median_positive_background:
Median positive background acoustic activity for the same site, year, and species. This provides a less extreme alternative to the mean when background activity has large outlier nights.

n_positive_nights:
Number of nights with positive background acoustic activity for the same site, year, and species.

n_zero_nights:
Number of nights with zero background acoustic activity for the same site, year, and species.

n_total_nights:
Total number of background-monitoring nights available for the same site, year, and species.

background_status:
Classification of the observed background value. Possible values are:

positive_background: observed background activity was greater than zero.

zero_background: a matching SM3 row existed, but background activity was zero.

missing_background: no matching SM3 row was found.

n_background_for_offset:
Final background denominator used to create the abundance offset. If n_background > 0, the observed value was used. If n_background was zero or missing, the mean positive background value for the same site, year, and species was used when available.

background_was_imputed:
Logical variable indicating whether n_background_for_offset was imputed. FALSE means the observed positive n_background value was used. TRUE means a zero or missing background value was replaced with the site-year-species mean positive background activity.

log_background_offset:
Log-transformed background denominator used as an offset in negative-binomial models. This variable was calculated as log(n_background_for_offset). It is on the log scale because the count models use a log link.

Notes:
The main response variable remains n_calls. Background acoustic activity can be incorporated in models as an offset using log_background_offset. If sampling effort differs among site-nights, effort_hours should also be included as an offset using log(effort_hours).


"
# Write the README content to a file
writeLines(readme_content, "data_for_analysis/rob_spkr_prep_v3/README_v3.txt")


# write data
write.csv(rob_db, file = 'data_for_analysis/rob_spkr_prep_v3/rob_db_v3.csv', row.names = F) 
write.csv(spkr_db, file = 'data_for_analysis/rob_spkr_prep_v3/spkr_db_v3.csv', row.names = F)








