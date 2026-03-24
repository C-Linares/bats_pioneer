
# =======================================================================
# Script Title:    rob_spkr_prep.R
# Project:         robomoth and speaker data analysis
# =======================================================================

# Description:
# this script helps prep the data for analysis, standardize, addpredictors, etc. 

# Author: Carlos Linares
# Date:   2026-03-3
# Contact: carlosgarcialina@u.boisestate.edu

# =======================================================================
# Inputs:
# - bat_robomoth.csv (home/r/bats_pioneer/data_for_analysis/bat_robomoth.csv)
# - bat_speakr.csv (home/r/bats_pioneer/data_for_analysis/bat_speakr.csv)

#  Outputs:
#  - ready to analyze data base ("data_for_analysis/rob_spkr_prep/rob_db.csv")
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


robomoth <- read_csv("data_for_analysis/robomoth_build/bat_robomoth.csv") %>% 
  clean_names()

str(robomoth)
colSums(is.na(robomoth)) # no missing data in robomoth

speakr <- read_csv("data_for_analysis/robomoth_build/spkr_all.csv") %>% 
  clean_names()

str(speakr)
colSums(is.na(speakr)) # some missing data in speakr, but not many.

unique(speakr$site) # check the unique sites in speakr
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
c_bugs <- read_csv("data_for_analysis/insect_wranglin/c_bugs.csv") %>%  # load insect data
  clean_names() %>%
  rename(yr = yrs) # safe rename
# add treatment. 

litsites<-c("iron01","iron03","iron05","long01","long03")


c_bugs$treatmt<-ifelse(c_bugs$site %in% litsites , "lit", "dark") # this makes a treatment variable.


# calculate mean by yr, trmt and site. I will use this to substitute the NA valeus from tha appeared when merging the bat data. 

c_bugs_mean <- c_bugs %>%
  group_by(yr,treatmt, site) %>%
  summarise(
    t_insect = mean(t_insect, na.rm = TRUE),
    t_lepidoptera = mean(t_lepidoptera, na.rm = TRUE)
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

summary(rob_db) #there is just one NA I can keep it that way. 

# merge with insects
# I calculate week for insects to merge with robomot data. 
rob_db<- rob_db %>%
  mutate(wk = week(noche)) # add week and year to rob_db for merging

rob_db <- rob_db %>%
  left_join(c_bugs, by = c("site", "wk", "year"= "yr")) # merge

# check for NAs
summary(rob_db)  #lots of NAs in t.insect and t.lepidoptera.

# I might need to check the insect data to see if the c_bugs is the best predictor to use for this analysis. 

# we remove the nas from insects by adding the average year insect count by site 

rob_db <- rob_db %>%
  left_join(c_bugs_mean, 
            by = c("year"="yr", "site"),
            suffix = c("", "_mean")) %>%  # rename mean columns directly
  mutate(
    t_insect = coalesce(t_insect, t_insect_mean),
    t_lepidoptera = coalesce(t_lepidoptera, t_lepidoptera_mean)
  ) %>%
  select(-t_insect_mean, -t_lepidoptera_mean)

summary(rob_db) # now there are no NAs in t.insect and t.lepidoptera.

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

# there's 3 diff treatment col keep just one.

rob_db <- rob_db %>%
  select( -treatmt, -treatmt.y) %>% 
  rename( trmt = treatmt.x) # rename the treatment column to be consistent with the speakr data.

summary(rob_db)

# remoe the zero_buzz column

rob_db <- rob_db %>%
  select(-zero_buzz)


# standardize -------------------------------------------------------------

#calculate Julian day from 'noche' and scale variables
rob_db<-rob_db %>% 
  mutate(
    jday = lubridate::yday(noche),  # Calculate Julian day from 'noche'
  )

variables_to_scale <- c(
  "avg_moonlight",
  "avg_twilight",
  "avg_illumination",
  "nit_avg_temp_c",
  "nit_avg_wspm_s",
  "t_lepidoptera",
  "t_insect",
  "jday"
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

# merge the ps_label with the rob_db and by sp_code

rob_db <- rob_db %>%
  left_join(species %>% select(sp_code, sp_label), by = c("sp" = "sp_code") )


glimpse(rob_db)



#  speakr data  -----------------------------------------------------------

# there's an error in the tags for sites on vizcaine. this needs to be fixed on robomoth_build.R script. (fixed)

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
  left_join(c_bugs, by = c("site", "wk", "year" = "yr"))

summary(spkr_db)

# 5. fill missing insect values with annual site mean
spkr_db <- spkr_db %>%
  left_join(
    c_bugs_mean,
    by = c("year" = "yr", "site"),
    suffix = c("", "_mean")
  ) %>%
  mutate(
    t_insect = coalesce(t_insect, t_insect_mean),
    t_lepidoptera = coalesce(t_lepidoptera, t_lepidoptera_mean)
  ) %>%
  select(-t_insect_mean, -t_lepidoptera_mean)

summary(spkr_db)

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

# remove unnecessary columns

names(spkr_db)

# there's 3 diff treatment col keep just one.

spkr_db <- spkr_db %>%
  select( -treatmt, -treatmt.y) %>% 
  rename( trmt = treatmt.x) # rename the treatment column to be consistent with the speakr data.

summary(spkr_db)



# standardize -------------------------------------------------------------

#calculate Julian day from 'noche' and scale variables
spkr_db<-spkr_db %>% 
  mutate(
    jday = lubridate::yday(noche),  # Calculate Julian day from 'noche'
  )

variables_to_scale <- c(
  "avg_moonlight",
  "avg_twilight",
  "avg_illumination",
  "nit_avg_temp_c",
  "nit_avg_wspm_s",
  "t_lepidoptera",
  "t_insect",
  "jday"
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






# outputs -----------------------------------------------------------------



# dir.create("data_for_analysis/prep_for_glm", showWarnings = FALSE) # just run if the dir is abscent




# Create a README file with information about the script
readme_content <- "Carlos Linares 3/23/2026 
this folde contains the database with predictors ready to run models with. 
the columns are as follow:

noche: date of the night
site: site of the observation
year: year of the observation
hi_f: high frequency bat call binomial not counts 
lo_f: low frequency bat call binomial not counts 
t_insect: total insect count summmarized insect counts by week, site and year. when data for that week was missing we added the mean insect count for that site and year.
t_lepidoptera: total lepidoptera count summmarized as insects counts by week, site and year. when data for that week was missing we added the mean lepidoptera count for that site and year.
nit_avg: average night time temperature for that night
nit_avg_wspm: average night time wind speed for that night
treatment : treatment of the site, either lit or dark.
c_buzz: corrected buzz counts taking in consideration both manually vetted buzz and auto_buzz counts. 
avg_illuminatio: both twilight and moon illumination
avf_moonlight: moonlight model prediction for that night. from the moonlit package (Kayba)


"
# Write the README content to a file
writeLines(readme_content, "data_for_analysis/rob_spkr_prep/README.txt")


# write data
write.csv(rob_db, file = 'data_for_analysis/rob_spkr_prep/rob_db.csv', row.names = F) 
write.csv(spkr_db, file = 'data_for_analysis/rob_spkr_prep/spkr_db', row.names = F)








