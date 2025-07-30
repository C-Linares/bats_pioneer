# ---------------------------
##
## Script name:  map_bat_activity.R
##
## Purpose of script: plot the bat activity as a percentage or a ratio across sites. 
##
## Author: Carlos Linares
##
## Date Created: 4/9/2025
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes: 
##   
##
## ---------------------------
## # inputs ------------------------------------------------------------------
#   
# bat data matrix with the activity index for the bat species.
# bm.ai<- read_csv('data_for_analysis/') # activity index for the bat species.


# outputs ----------------------

# map of change of activity across years. 



# libraries  --------------------------------------------------------------

library(pacman)
pacman::p_load(
  "tidyverse",
  "magrittr",
  "lme4",
  "sjPlot",
  "ggeffects",
  "car",
  "glmmTMB",#model
  "corrplot",
  "effects",
  "reshape2",
  "DHARMa",#check assumptions 
  "marginaleffects", #plot models
  "MuMIn", #models
  "performance",
  "viridis",
  "data.table"
)



# data --------------------------------------------------------------------

bat.ai<-read_csv('data_for_analysis/prep_for_glmm/bm.miller.day.csv', col_names = T) # activity index for the bat species.

# remove the first column that is unecessary
bat.ai<-bat.ai[,-1]

#date as date
bat.ai$noche<-ymd(bat.ai$noche)
# make year column 
bat.ai$year<-year(bat.ai$noche)

# fiter out all observatons for noise and NoID in the species column
bat.ai<-bat.ai %>% filter(sp != "NoID" & sp != "Noise")

summary(bat.ai)
str(bat.ai)

# site coordinates

site.coords<-read_csv('data_for_analysis/sites_coordinates.csv', col_names = T) # activity index for the bat species.


# data wrangling -----------------------------------------------------------
# I need to calculate total bat activity by year and then the percentage of that total that each species was active. 

# First, calculate the total activity minutes by year, site, and species
bat_summary <- bat.ai %>%
  group_by(year, site, sp) %>%
  summarize(total_activity_sp = sum(activity_min)) %>%
  ungroup()

# Then calculate the total activity minutes by year and site (for all species)
site_totals <- bat.ai %>%
  group_by(year, site) %>%
  summarize(total_activity_site = sum(activity_min)) %>%
  ungroup()

# Join the two data sets
bat_final <- bat_summary %>%
  left_join(site_totals, by = c("year", "site"))


# calculate the percentage of activity by species by site and year. 

bat_final <- bat_final %>%
  mutate(percentage = (total_activity_sp / total_activity_site) * 100)


# coordinates merge with bat final

bat_final<-bat_final %>%
  left_join(site.coords, by = c("site" = "site"))

bat_change<-bat_final %>% 
  arrange(site, year) %>%
  group_by(site) %>%
mutate(
  pct_change = (total_activity_site - lag(total_activity_site)) / lag(total_activity_site) * 100
) %>%
  summarise(
    mean_annual_pct_change = mean(pct_change, na.rm = TRUE),
    lat = first(lat),
    lon = first(lon)
  )


# convert data to spatial
library(sf)
site_sf <- st_as_sf(bat_final, coords = c("lon", "lat"), crs = 4326)

# plot the data on a map

library(ggplot2)

ggplot(site_sf) +
  geom_sf(aes(color = percentage, size = abs(percentage))) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(color = "Call Decline Rate", size = "Decline Magnitude")

ggplot(site_sf) +
  geom_sf(aes(color = percentage, size = abs(percentage)), alpha = 0.8) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Change in Bat Calls"
  ) +
  scale_size(range = c(2, 8), name = "Magnitude of Change") +
  theme_minimal() +
  labs(
    title = "Bat Call Decline Across Sites (2021–2023)",
    subtitle = "Negative values indicate a decline in bat calls",
    caption = "Simulated data")



# predictions from the modesl ----------------------------------------------

preds<- read_csv('data_for_analysis/glmm_v2/preds.csv', col_names = T) # activity index for the bat species.

# join preds with site coords

preds<-preds %>%
  left_join(site.coords, by = c("site" = "site"))

site_sf <- st_as_sf(preds, coords = c("lon", "lat"), crs = 4326)


# graphs ------------------------------------------------------------------


ggplot(site_sf) +
  geom_sf(aes(color = estimate, size = estimate), alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "Predicted Activity") +
  scale_size(range = c(1, 12), name = "Predicted Activity") +
  facet_wrap(~ yr_s) +
  theme_minimal() +
  labs(
    title = "Predicted Bat Activity by Site and Year",
    subtitle = "Model-based estimates of bat calls (GLMM)",
    caption = "Each point represents a site. Color and size reflect predicted activity."
  )

ggplot(site_sf) +
  geom_sf(aes(color = estimate, size = estimate), alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "Predicted Activity") +
  scale_size(range = c(1, 12), name = "Predicted Activity") +
  facet_wrap(~ yr_s) +
  scale_x_continuous(breaks = c(-113.5, -113.0, -112.5)) +  # ← Only 3 tick marks
  theme_minimal() +
  labs(
    title = "Predicted Bat Activity by Site and Year",
    subtitle = "Model-based estimates of bat calls (GLMM)",
    caption = "Each point represents a site. Color and size reflect predicted activity."
  )

unique(site_sf$year)


library(gganimate)

site_sf <- site_sf %>% 
  filter(!is.na(year)) %>%
  mutate(year = as.numeric(as.character(year)))  # or as.factor(year) if categorical

ggplot(site_sf) +
  geom_sf(aes(color = estimate, size = estimate), alpha = 0.8) +
  scale_color_viridis_c(option = "magma", name = "Predicted Activity") +
  scale_size(range = c(2, 8), name = "Predicted Activity") +
  theme_minimal() +
  labs(
    title = "Predicted Bat Activity by Site",
    subtitle = "Year: {frame_time}",
    caption = "Animated from GLMM predictions"
  ) +
  transition_time(year) +
  ease_aes("linear")

site_change <- preds %>%
  group_by(site) %>%
  summarise(
    estimate_2021 = estimate[year == 2021],
    estimate_2023 = estimate[year == 2023],
    pct_change = 100 * (estimate_2023 - estimate_2021) / estimate_2021,
    lat = first(lat),
    lon = first(lon)
  ) %>%
  filter(!is.na(pct_change))

# Convert to sf
site_change_sf <- st_as_sf(site_change, coords = c("lon", "lat"), crs = 4326)

# Plot % change on a map
ggplot(site_change_sf) +
  geom_sf(aes(color = pct_change, size = abs(pct_change)), alpha = 0.8) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "% Change (2021–2023)") +
  theme_minimal() +
  labs(
    title = "Percent Change in Predicted Bat Activity (2021–2023)",
    subtitle = "Negative values indicate a decline",
    caption = "From GLMM model predictions"
  )


# bargraph of the percentage of activity by species by yera

bat_final %>%
  ggplot(aes(x = year, y = percentage, fill = sp)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Bat Activity Percentage by Species and Year Across Sites",
       x = "Year",
       y = "Percentage of Activity") +
  theme_minimal() +
  scale_fill_viridis_d() +
  # facet_wrap(~year) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(sf)

# Convert to spatial object
site_sf <- st_as_sf(bat_change, coords = c("lon", "lat"), crs = 4326)

# Plot
ggplot(site_sf) +
  geom_sf(aes(color = mean_annual_pct_change, size = abs(mean_annual_pct_change)), alpha = 2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "Mean Annual % Change") +
  scale_size(range = c(2, 8), name = "Change Magnitude") +
  theme_minimal() +
  labs(
    title = "Mean Annual Percent Change in Bat Activity (2021–2023)",
    subtitle = "Negative values indicate decline in activity",
    caption = "Source: Bat monitoring data"
  )
