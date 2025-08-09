# ---------------------------
##
## Script name:  delta_sunset.R
##
## Purpose
# this script will is to answer the question: Do bat vocal activity patterns varies between experimental sites?
# 
## Author: Carlos Linares, 
## Date Created: 8/8/2025
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes: Bat data product of the prep_for_glmm.R script
##        sunrise and sunset data: from sunrise scripts/sunrise_sunset_cal.R from the Bird pioneer project. 
##   
##
## ---------------------------
## # inputs ------------------------------------------------------------------
#' - bat_combined, file = 'data_for_analysis/prep_for_glmm/bat_combined.csv'
#

# outputs ----------------------

# plots comparing activity between sites 


# libraries  --------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  "tidyverse",
  "lubridate",
  "here",# for reproducible file paths
  "janitor"
  )



# data --------------------------------------------------------------------

bat_combine<- read_csv("data_for_analysis/prep_for_glmm/bat_combined.csv") %>% 
  clean_names()
str(bat_combine)

# remove effort and jday. 

bat_combine<- bat_combine %>% 
  select(-c(eff_hrs, jday, trmt_bin))

# time zone
attr(bat_combine$date_time, "tzone") #the time is utc might need to change it to Americe/Denver. 

bat_combine<- bat_combine %>% 
  mutate(
    dtime= force_tz(date_time, tzone = "America/Denver") # keeps the time but makes it mtd
  )

# remove noise and NoID rows

bat_combine<- bat_combine %>%
  filter(!sp %in% c("Noise", "NoID"))


# Prepare sunset times data -----------------------------------------------

# Note: You'll need to create or load sunset times for your sites and dates
# This is a placeholder - replace with your actual sunset data loading


sunset_data<- read_csv('data_for_analysis/sunrise_sunset/sunrise_sunset_.csv') %>% 
  clean_names()
str(sunset_data)

sunset_data <- sunset_data %>%
  mutate(date = as.Date(date))

# Merge bat data with sunset times ----------------------------------------

# Join bat data with sunset times
bat_with_sunset <- bat_combine %>%
  left_join(sunset_data, by = c("site", "date_time"="date"))

missing_sunset_data <- bat_combine %>%
  anti_join(sunset_data, by = c("site", "noche" = "date"))

# Calculate time to sunset ------------------------------------------------

cat("Calculating time differences from sunset...\n")

bat_with_sunset <- bat_with_sunset %>%
  mutate(
    # Ensure date_time is POSIXct
    date_time = as.POSIXct(date_time, tz = "UTC"),
    
    # Ensure sunset is POSIXct  
    sunset = as.POSIXct(sunset, tz = "UTC"),
    
    # Calculate time difference from sunset in minutes
    # Negative values = before sunset, positive = after sunset
    time_to_sunset = as.numeric(difftime(date_time, sunset, units = "mins"))
  )

# Data quality checks -----------------------------------------------------

cat("Performing data quality checks...\n")

# Summary statistics
sunset_summary <- bat_with_sunset %>%
  filter(!is.na(time_to_sunset)) %>%
  summarise(
    n_records = n(),
    min_time = min(time_to_sunset),
    max_time = max(time_to_sunset),
    mean_time = mean(time_to_sunset),
    median_time = median(time_to_sunset),
    records_before_sunset = sum(time_to_sunset < 0),
    records_after_sunset = sum(time_to_sunset > 0)
  )

print(sunset_summary)

# Check for extreme values (might indicate data issues)
extreme_values <- bat_with_sunset %>%
  filter(!is.na(time_to_sunset)) %>%
  filter(abs(time_to_sunset) > 12 * 60) %>%  # More than 12 hours from sunset
  nrow()

if (extreme_values > 0) {
  warning("Found ", extreme_values, " records more than 12 hours from sunset. Check data quality.")
}

# Create final dataset ----------------------------------------------------

# Select final columns (modify as needed)
bat_final <- bat_with_sunset %>%
  select(
    sp, PULSES, site, noche, date_time, yr, treatmt, trmt_bin, 
    jday, eff.hrs, sunset, time_to_sunset
  ) %>%
  # Optional: filter out records with missing sunset data
  filter(!is.na(time_to_sunset))

cat("Final dataset contains", nrow(bat_final), "records with sunset calculations.\n")

# Visualization -----------------------------------------------------------

# Quick visualization to check the data
p1 <- bat_final %>%
  filter(abs(time_to_sunset) <= 480) %>%  # Focus on ±8 hours from sunset
  ggplot(aes(x = time_to_sunset)) +
  geom_histogram(bins = 50, alpha = 0.7, fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  labs(
    title = "Distribution of Bat Calls Relative to Sunset",
    x = "Time to Sunset (minutes)",
    y = "Number of Calls",
    subtitle = "Red line = sunset time"
  ) +
  theme_minimal()

print(p1)

# By treatment if available
if ("treatmt" %in% names(bat_final)) {
  p2 <- bat_final %>%
    filter(abs(time_to_sunset) <= 480) %>%
    ggplot(aes(x = time_to_sunset, fill = treatmt)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
    labs(
      title = "Bat Calls Relative to Sunset by Treatment",
      x = "Time to Sunset (minutes)",
      y = "Number of Calls"
    ) +
    theme_minimal() +
    facet_wrap(~treatmt)
  
  print(p2)
}

# Export results ----------------------------------------------------------

# Optional: save the processed dataset
# write_csv(bat_final, "data/processed/bat_combine_with_sunset.csv")

cat("Analysis complete! The 'bat_final' dataset now contains the 'time_to_sunset' column.\n")

# Return the final dataset
bat_final