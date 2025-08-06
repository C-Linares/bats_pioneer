  #=============================================================================
  # Script Name:    glmm_v3.R
  # Purpose:        Fit GLMMs for bat paper using weather, and moon data, and light spectra. 
  #
  # Version:        v2 - Includes:
  #                   - 
  #                   - Moonlight estimates from the 'moonlit' package recalculated using midnight. 
  #
  # Author:         Carlos Linares  
  # Collaborator:   Jen Cruz
  # Created:        2025-02-11
  # Contact:        carlosgarcialina@u.boisestate.edu
  #
  # Notes:
  #   - Data prepped using prep_for_glmm.R
  #   - Outputs include marginal effect and random effects plots
  #
  # Inputs:
  #   - data_for_analysis/prep_for_glmm/bm.csv         : Bat call counts (daily summary)
  #   - data_for_analysis/...                          : Activity index and trait data
  #
  # Outputs:
  #   - Visualizations of model predictions
  #   - Diagnostic plots for model assumptions
  # =============================================================================

  

#  libraries-------------------------------------------------------------------------

  # List of packages
  packages <- c(
    "tidyverse",
    "magrittr",
    "lme4",
    "sjPlot",
    "ggeffects",
    "car",
    "glmmTMB",       # model
    "corrplot",
    "effects",
    "reshape2",
    "DHARMa",        # check assumptions 
    "marginaleffects", # plot models
    "MuMIn",         # models
    "performance",
    "viridis",
    "data.table",
    "janitor",       # to clean names
    "patchwork"
  )
  
  # Try loading each package one by one
  for (pkg in packages) {
    message("Attempting to load: ", pkg)
    tryCatch({
      library(pkg, character.only = TRUE)
      message("Successfully loaded: ", pkg)
    }, error = function(e) {
      message("Failed to load: ", pkg, "\nError: ", e$message)
    })
    # Optional: Pause between packages to observe behavior
    readline(prompt = "Press [Enter] to continue...")
  }
  

# Data --------------------------------------------------------------------

  # Load the data
  bm <- read_csv("data_for_analysis/prep_for_glmm/bm.csv") %>%
    clean_names()
  
  # Check the structure of the data
  str(bm)
  
  # Check for missing values
  colSums(is.na(bm))
  
  filtered_bm <- bm %>%
    rename(sp = auto_id) %>%
    filter(!sp %in% c("Noise", "NoID")) # remove noise and noid from the analysis. 


  # miller activity index
  
  bm_ai <- read_csv("data_for_analysis/prep_for_glmm/bm.miller.day.csv") %>% 
    clean_names() %>%
    filter(!sp %in% c("Noise", "NoID"))
  
  
# predictors --------------------------------------------------------------


  # Creaters weather (night).
  
  crmo.wet.night <- read_csv("data_for_analysis/weather/craters_weater/craters_night.csv") %>% 
    clean_names()# load creaters night weather. 
  
  
  # Moon
  
  moon.int<- read_csv('data_for_analysis/moon_pred/moon.int.csv') %>%
    clean_names()
  summary(moon.int) # check the structure of the moon data)
  
# convert date times from UTC to Denver/America

  moon.int$denver_time <- with_tz(moon.int$date, tzone = "America/Denver") # Convert to Denver time zone

  attr(moon.int$denver_time, "tzone") # check the timezomne is correct
  
  # summarize moonlight by date but conditional moon_alt_degrees>0

  # Step 1: Filter data where the moon is above the horizon
  moon_filtered <- moon.int %>%
    filter(moon_alt_degrees > 0)
  
  # Step 2: Create a new 'noche' variable (just the date part of the timestamp) but also makes nights any time stamps that are less than 9 am.
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
  
  c_bugs <- read_csv("data_for_analysis/insect_wranglin/c_bugs.csv") %>%  # load insect data.
    clean_names() %>%
    rename(yr = yrs) # safe rename
  
  
  # light
  
  light <- read_csv("data_for_analysis/lights/lightspectra_pioneer.csv") %>%
    clean_names() %>%
    filter(vert_horiz == "Horizontal") %>%                          # Keep only horizontal measures.
    mutate(mwatts = rowMeans(across(c(watts_m1, watts_m2, watts_m3)), na.rm = TRUE)) %>%  # Mean of watts
    select(site, lux, yr, mwatts)    # Select relevant columns  
  
  
# Merge datasets ------------------------------------------------------------

  
  
  

    