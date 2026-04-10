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
  "janitor",
  "ovelap",
  "purrr",
  "patchwork"
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
  left_join(sunset_data, by = c("noche"="date"))

missing_sunset_data <- bat_combine %>%
  anti_join(sunset_data, by = c("noche" = "date"))

# Calculate time to sunset ------------------------------------------------

bat_with_sunset <- bat_with_sunset %>%
  mutate(
    # Ensure date_time is POSIXct in Denver time
    date_time = as.POSIXct(date_time, tz = "America/Denver"),
    
    # Ensure sunset is POSIXct in Denver time
    sunset = as.POSIXct(sunset, tz = "America/Denver"),
    
    # Calculate time difference from sunset in minutes
    time_to_sunset = as.numeric(difftime(date_time, sunset, units = "mins")),
    deltasunset = as.numeric(difftime(date_time, sunset, units = "hours"))
  )


# transfomr to radians

bat_with_sunset$deltasunset_rad <-
  ifelse(
    bat_with_sunset$deltasunset < 0,
    (24 + bat_with_sunset$deltasunset) / 24,
    bat_with_sunset$deltasunset / 24
  ) * 2 * pi



# see how many records there are by species.

bats<-read_csv('data_for_analysis/Species_bats.csv') %>% 
  clean_names()

bat_with_sunset<-left_join(bat_with_sunset, bats, by=c("sp"="six_letter_species_code"))

sp_site_treat <- table(bat_with_sunset$sp, bat_with_sunset$treatmt)
sp_site_treat

# species labels for plots


spp<-unique(bat_with_sunset$sp)
spplabs<-unique(bat_with_sunset$common_name)


for(i in seq_along(spp)) {
  sub_bm <- bat_with_sunset %>% filter(sp == spp[i])
  a_lit <- bat_with_sunset %>% filter(treatmt == "lit") %>% pull(deltasunset_rad)
  a_dark <- bat_with_sunset %>% filter(treatmt == "dark") %>% pull(deltasunset_rad)
  
  # Only plot if at least one group has more than 1 point
  if(sum(!is.na(a_lit)) > 1 & sum(!is.na(a_dark)) > 1) {
    ylim_max <- max(
      density(a_lit, na.rm = TRUE)$y,
      density(a_dark, na.rm = TRUE)$y
    )
    plot(density(a_lit, na.rm = TRUE), 
         main = spplabs[i], 
         xlab = "Time after sunrise (radians)",
         col = "orange", lwd = 2, ylim = c(0, ylim_max))
    rug(a_lit, col = "orange")
    lines(density(a_dark, na.rm = TRUE), col = "blue", lwd = 2)
    rug(a_dark, col = "blue")
    legend("topright", legend = c("Lit", "Dark"), col = c("orange", "blue"), lwd = 1)
  } else if(sum(!is.na(a_lit)) > 1) {
    plot(density(a_lit, na.rm = TRUE), 
         main = paste0(spplabs[i], " (Lit only)"), 
         xlab = "Time after sunrise (radians)",
         col = "orange", lwd = 2)
    rug(a_lit, col = "orange")
    legend("topright", legend = c("Lit"), col = c("orange"), lwd = 2)
  } else if(sum(!is.na(a_dark)) > 1) {
    plot(density(a_dark, na.rm = TRUE), 
         main = paste0(spplabs[i], " (Dark only)"), 
         xlab = "Time after sunrise (radians)",
         col = "blue", lwd = 2)
    rug(a_dark, col = "blue")
    legend("topright", legend = c("Dark"), col = c("blue"), lwd = 2)
  } else {
    plot.new()
    title(main = paste0(spplabs[i], "\nNot enough data for density"))
  }
}



# Density Plots ----------------------------------------------------------


bat_with_sunset %>%
  filter(sp %in% spp, !is.na(deltasunset_rad)) %>%
  ggplot(aes(x = deltasunset_rad, color = treatmt, fill = treatmt)) +
  geom_density(alpha = 0.3, adjust = 1.2) +
  geom_rug(aes(color = treatmt), sides = "b") +
  facet_wrap(~ sp, labeller = labeller(sp = setNames(spplabs, spp))) +
  labs(
    x = "Time after sunrise (radians)",
    y = "Density",
    color = "Treatment",
    fill = "Treatment"
  ) +
  scale_color_manual(values = c("lit" = "orange", "dark" = "blue")) +
  scale_fill_manual(values = c("lit" = "orange", "dark" = "blue")) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

# in the graph below 
density_graph<-bat_with_sunset %>%
  filter(sp %in% spp, !is.na(deltasunset_rad)) %>%
  ggplot(aes(x = deltasunset_rad, color = treatmt, fill = treatmt)) +
  geom_density(alpha = 0.3, adjust = 1.2) +
  geom_rug(aes(color = treatmt), sides = "b") +
  facet_wrap(~ sp, labeller = labeller(sp = setNames(spplabs, spp))) +
  labs(
    x = "Delta sunset",
    y = "Density",
    color = "Treatment",
    fill = "Treatment"
  ) +
  scale_color_manual(values = c("lit" = "orange", "dark" = "blue")) +
  scale_fill_manual(values = c("lit" = "orange", "dark" = "blue")) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

density_graph

ggsave(density_graph, file="figures/bat_delta_sunset/desnity_graph_v1.tiff", width = 10, height = 9)



# improved graph just looking at certain hours after sunset
density_graph <- bat_with_sunset %>%
  filter(sp %in% spp, !is.na(deltasunset_rad)) %>%
  ggplot(aes(x = deltasunset_rad, color = treatmt, fill = treatmt)) +
  geom_density(alpha = 0.3, adjust = 1.2) +
  # Added jitter and transparency so the rug lines don't completely overlap
  # geom_rug(aes(color = treatmt), sides = "b", alpha = 0.5, position = "jitter") +
  # Added scales = "free_y" to let each bat species scale its own density height
  facet_wrap(~ sp, labeller = labeller(sp = setNames(spplabs, spp)), scales = "free_y") +
  # Zoom in on the 0-4 range without recalculating/chopping the density curves
  coord_cartesian(xlim = c(0, 4)) +
  labs(
    x = "Delta sunset",
    y = "Density",
    color = "Treatment",
    fill = "Treatment"
  ) +
  scale_color_manual(values = c("lit" = "orange", "dark" = "blue")) +
  scale_fill_manual(values = c("lit" = "orange", "dark" = "blue")) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    # Optional: removes the y-axis text since density values are relative and vary by facet now
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank()
  )

density_graph


ggsave(
  filename = "figures/bat_delta_sunset/delta_sunset_density.png",
  plot = density_graph,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)


# Overlap -----------------------------------------------------------------

library(dplyr)
library(purrr)
library(overlap)

# lets separate the same species calls in to vectors one for lit and one for dark sites
# then we use the function overlapEst to estimate the overlap. 

# Antrozous pallidus

ap_lit <- bat_with_sunset %>%
  filter(scientific_name == "Antrozous pallidus", treatmt == "lit") %>% pull(deltasunset_rad)

ap_dark <- bat_with_sunset %>%
  filter(scientific_name == "Antrozous pallidus", treatmt == "dark") %>% pull(deltasunset_rad)



# Overlap Estimation ------------------------------------------------------


par( mfrow = c(1,1), cex = 1.7, lwd = 2, bty = "l", pty = "m" )
overlapPlot( ap_lit, ap_dark, main = " Antrozous pallidus overlap", lty =c(1,20),
             col = c(1,4),xlab = "Time after sunset", rug =T )
overlapEst(ap_lit, ap_dark)


# Function to compute overlap for a given species
compute_overlap <- function(df, sp_name) {
  lit <- df %>%
    filter(scientific_name == sp_name, treatmt == "lit") %>%
    pull(deltasunset_rad)
  
  dark <- df %>%
    filter(scientific_name == sp_name, treatmt == "dark") %>%
    pull(deltasunset_rad)
  
  # only compute if both have data
  if (length(lit) > 5 & length(dark) > 5) {
    est <- overlap::overlapEst(lit, dark, type = "all")
    return(est)
  } else {
    return(NA_real_)
  }
}

# Apply to all species
overlap_results <- bat_with_sunset %>%
  distinct(scientific_name) %>%
  pull(scientific_name) %>%
  set_names() %>% 
  map(~ compute_overlap(bat_with_sunset, .x)) %>%
  bind_rows(.id = "species")

overlap_results


# Example: pivot longer for plotting
overlap_long <- overlap_results %>%
  pivot_longer(
    cols = everything(),
    names_to = "species",
    values_to = "overlap"
  ) %>%
  mutate(species = factor(species))

# Plot as barplot with error bars if multiple rows per species
p1<-ggplot(overlap_long, aes(x = species, y = overlap)) +
  geom_point()+
  # geom_violin()+
  # geom_boxplot(fill = "skyblue", alpha = 0.6) +   # summarizes the 3 rows per species
  # geom_jitter(width = 0.2, alpha = 0.7, size = 2) + # shows individual values
  coord_flip() +
  labs(
    title = "Activity Overlap between Lit and Dark Treatments",
    y = "Overlap Index (Δ)",
    x = "Species"
  ) +
  theme_minimal(base_size = 14)
p1

ggsave(p1, file = "figures/bat_delta_sunset/overlap_v1.tiff", width = 10, height = 6)


density_graph+p1



# new version with help from chatgpt

# delta_sunset_clean.R
#
# Purpose:
# Compare bat activity timing between lit and dark sites
# using time relative to sunset in HOURS (not radians).
# ---------------------------

# libraries --------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse,
  lubridate,
  here,
  janitor,
  patchwork
)

# data -------------------------------------------------------------------
bat_combine <- read_csv("data_for_analysis/prep_for_glmm/bat_combined.csv") %>%
  clean_names()

# remove variables not needed here
bat_combine <- bat_combine %>%
  select(-c(eff_hrs, jday, trmt_bin))

# remove noise and NoID rows
bat_combine <- bat_combine %>%
  filter(!sp %in% c("Noise", "NoID"))

# ------------------------------------------------------------------------
# CHANGE 1:
# Make sure date_time is in the correct timezone.
# Use with_tz() if the timestamp is already stored in UTC and you want the
# equivalent local time in America/Denver.
# Use force_tz() only if the clock time is already local but missing the label.
# ------------------------------------------------------------------------

bat_combine <- bat_combine %>%
  mutate(
    date_time = force_tz(date_time, tzone = "America/Denver")
  )

# sunset data ------------------------------------------------------------
sunset_data <- read_csv("data_for_analysis/sunrise_sunset/sunrise_sunset_.csv") %>%
  clean_names() %>%
  mutate(date = as.Date(date))

# merge ------------------------------------------------------------------
bat_with_sunset <- bat_combine %>%
  left_join(sunset_data, by = c("noche" = "date"))

missing_sunset_data <- bat_combine %>%
  anti_join(sunset_data, by = c("noche" = "date"))

# ------------------------------------------------------------------------
# CHANGE 2:
# Keep time relative to sunset in HOURS.
# This is the variable you should use in plots.
# Negative values = before sunset
# Positive values = after sunset
# ------------------------------------------------------------------------

bat_with_sunset <- bat_with_sunset %>%
  mutate(
    sunset = as.POSIXct(sunset, tz = "America/Denver"),
    hours_since_sunset = as.numeric(difftime(date_time, sunset, units = "hours")),
    minutes_since_sunset = as.numeric(difftime(date_time, sunset, units = "mins"))
  )

# ------------------------------------------------------------------------
# CHANGE 3:
# If you still want the circular version for some other analysis, keep it.
# But DO NOT use it for the main figure if you want interpretable axes.
# ------------------------------------------------------------------------

bat_with_sunset <- bat_with_sunset %>%
  mutate(
    deltasunset_rad = ifelse(
      hours_since_sunset < 0,
      (24 + hours_since_sunset) / 24,
      hours_since_sunset / 24
    ) * 2 * pi
  )

# add species labels -----------------------------------------------------
bats <- read_csv("data_for_analysis/Species_bats.csv") %>%
  clean_names()

bat_with_sunset <- left_join(
  bat_with_sunset,
  bats,
  by = c("sp" = "six_letter_species_code")
)

spp <- unique(bat_with_sunset$sp)
spplabs <- unique(bat_with_sunset$common_name)

# ------------------------------------------------------------------------
# CHANGE 4:
# Main plot uses HOURS, not radians
# Also lower density smoothing slightly to preserve peaks better
# ------------------------------------------------------------------------

density_graph_hours <- bat_with_sunset %>%
  filter(sp %in% spp, !is.na(hours_since_sunset)) %>%
  ggplot(aes(x = hours_since_sunset, color = treatmt, fill = treatmt)) +
  geom_density(alpha = 0.3, adjust = 0.8) +
  facet_wrap(
    ~ sp,
    labeller = labeller(sp = setNames(spplabs, spp)),
    scales = "free_y"
  ) +
  coord_cartesian(xlim = c(-2, 10)) +
  labs(
    x = "Hours since sunset",
    y = "Density",
    color = "Treatment",
    fill = "Treatment"
  ) +
  scale_x_continuous(
    breaks = seq(-2, 10, by = 2)
  ) +
  scale_color_manual(values = c("lit" = "orange", "dark" = "blue")) +
  scale_fill_manual(values = c("lit" = "orange", "dark" = "blue")) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

density_graph_hours

# save -------------------------------------------------------------------
ggsave(
  filename = "figures/bat_delta_sunset/density_graph_hours_v1.tiff",
  plot = density_graph_hours,
  width = 10,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------------------
# OPTIONAL VERSION:
# If you want more visible peaks, use a histogram instead of density
# ------------------------------------------------------------------------

hist_graph_hours <- bat_with_sunset %>%
  filter(sp %in% spp, !is.na(hours_since_sunset)) %>%
  ggplot(aes(x = hours_since_sunset, fill = treatmt)) +
  geom_histogram(
    position = "identity",
    alpha = 0.4,
    binwidth = 0.5
  ) +
  facet_wrap(
    ~ sp,
    labeller = labeller(sp = setNames(spplabs, spp)),
    scales = "free_y"
  ) +
  coord_cartesian(xlim = c(-2, 10)) +
  labs(
    x = "Hours since sunset",
    y = "Number of calls",
    fill = "Treatment"
  ) +
  scale_x_continuous(
    breaks = seq(-2, 10, by = 2)
  ) +
  scale_fill_manual(values = c("lit" = "orange", "dark" = "blue")) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

hist_graph_hours




library(dplyr)
library(purrr)
library(tidyr)
library(overlap)

compute_overlap <- function(df, sp_name) {
  
  lit <- df %>%
    filter(scientific_name == sp_name, treatmt == "lit") %>%
    pull(deltasunset_rad) %>%
    na.omit()
  
  dark <- df %>%
    filter(scientific_name == sp_name, treatmt == "dark") %>%
    pull(deltasunset_rad) %>%
    na.omit()
  
  n_lit  <- length(lit)
  n_dark <- length(dark)
  
  # only compute if both groups have enough data
  if (n_lit > 5 && n_dark > 5) {
    
    # choose estimator depending on sample size
    estimator <- if (min(n_lit, n_dark) < 75) "Dhat1" else "Dhat4"
    
    overlap_est <- overlapEst(lit, dark, type = estimator)
    
    tibble(
      species = sp_name,
      n_lit = n_lit,
      n_dark = n_dark,
      estimator = estimator,
      overlap = overlap_est
    )
    
  } else {
    
    tibble(
      species = sp_name,
      n_lit = n_lit,
      n_dark = n_dark,
      estimator = NA_character_,
      overlap = NA_real_
    )
  }
}

overlap_results <- bat_with_sunset %>%
  distinct(scientific_name) %>%
  pull(scientific_name) %>%
  map_dfr(~ compute_overlap(bat_with_sunset, .x))

overlap_results

compute_overlap_ci <- function(df, sp_name, nboot = 1000) {
  
  lit <- df %>%
    filter(scientific_name == sp_name, treatmt == "lit") %>%
    pull(deltasunset_rad) %>%
    na.omit()
  
  dark <- df %>%
    filter(scientific_name == sp_name, treatmt == "dark") %>%
    pull(deltasunset_rad) %>%
    na.omit()
  
  n_lit  <- length(lit)
  n_dark <- length(dark)
  
  if (n_lit > 5 && n_dark > 5) {
    
    estimator <- if (min(n_lit, n_dark) < 75) "Dhat1" else "Dhat4"
    est <- overlapEst(lit, dark, type = estimator)
    
    # bootstrap
    boot_vals <- bootstrap(lit, dark, nb = nboot, type = estimator)
    ci <- quantile(boot_vals, probs = c(0.025, 0.975), na.rm = TRUE)
    
    tibble(
      species = sp_name,
      n_lit = n_lit,
      n_dark = n_dark,
      estimator = estimator,
      overlap = est,
      conf.low = ci[[1]],
      conf.high = ci[[2]]
    )
    
  } else {
    
    tibble(
      species = sp_name,
      n_lit = n_lit,
      n_dark = n_dark,
      estimator = NA_character_,
      overlap = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_
    )
  }
}

overlap_results <- bat_with_sunset %>%
  distinct(scientific_name) %>%
  pull(scientific_name) %>%
  map_dfr(~ compute_overlap_ci(bat_with_sunset, .x, nboot = 1000))

overlap_results


library(ggplot2)

p_overlap <- overlap_results %>%
  filter(!is.na(overlap)) %>%
  mutate(species = reorder(species, overlap)) %>%
  ggplot(aes(x = species, y = overlap)) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.2,
    color = "grey40"
  ) +
  geom_point(size = 2, color = "black") +
  coord_flip() +
  labs(
    title = "Activity overlap between lit and dark treatments",
    x = "Species",
    y = "Overlap coefficient (Δ)"
  ) +
  theme_minimal(base_size = 14)

p_overlap
