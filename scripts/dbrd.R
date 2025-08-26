# ---------------------------
##
## Script name:  bat_diversity_v2.R
##
## Purpose
#  We want to assess if the bat community changes with the light treatment. 
# 
## Author: Carlos Linares, 
## Date Created:8/25/2025
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes: Bat data product of the prep_for_glmm.R script
##        analysis following the methods from https://esajournals-onlinelibrary-wiley-com.libproxy.boisestate.edu/share/KSJBJ6ZJGCCGVC4TBJPR?target=10.1002/ecy.70128
##   
##
## ---------------------------
## # inputs ------------------------------------------------------------------
# - bm<- read_csv("data_for_analysis/prep_for_glmm/bm.csv")
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

# bat data 
bm<- read_csv("data_for_analysis/prep_for_glmm/bm.csv") %>% 
  clean_names()

bat_combined<-read_csv("data_for_analysis/prep_for_glmm/bat_combined.csv")

#effort

effort_days <- bat_combined %>%
  group_by(site, yr) %>%
  summarise(
    stard = min(noche),
    endd = max(noche),
    eff.days = as.numeric(difftime(max(noche), min(noche), units = "days"))
  )

total_effort <- effort_days %>%
  group_by(site) %>%
  summarise(total_effort = sum(eff.days, na.rm = TRUE))


# light data 

light<- read_csv("data_for_analysis/lights/lightspectra_pioneer.csv") %>% 
  clean_names() %>% 
  filter(vert_horiz == "Horizontal") %>%  # select just horizontal measures
  select("site","vert_horiz","watts_m1","watts_m2","watts_m3","lux","yr") %>% # keep this cols
  mutate(mwatts = rowMeans(across(c(watts_m1, watts_m2, watts_m3)), na.rm = TRUE)) # calculate mean watts. 

# calculate mean light by site

mlight <- light %>%
  group_by(site) %>%
  summarise(mean_mwatts = mean(mwatts, na.rm = TRUE))


# We need a site × species matrix with abundances.

bat_comm <- bm %>%
  group_by(site, auto_id) %>%     # site = sampling unit, auto_id = species
  summarise(abundance = sum(n), .groups = "drop") %>%
  pivot_wider(names_from = auto_id, values_from = abundance, values_fill = 0) %>% 
  select(-c(NoID, Noise)) # remove Noise and NoID

# Check matrix
head(bat_comm)

comm_matrix <- bat_comm %>% # make column site into rownames 
  arrange(site) %>%
  column_to_rownames("site")

# standardize and divide by number of days sampled after 3 years. 

comm_matrix_std <- comm_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("site") %>%
  left_join(total_effort, by = "site") %>%
  mutate(across(-c(site, total_effort), ~ .x / total_effort)) %>%
  select(-total_effort) %>%
  tibble::column_to_rownames("site")

comm_matrix_std

# Create Bray-Curtis distance matrix
bray_dist <- vegdist(comm_matrix_std, method = "bray")


# Run partial dbRDA---------------------------------------------------
# Example: test effect of treatment + year while controlling for site clusters


# Full model with predictors
dbrda_full <- dbrda(bray_dist ~ mean_mwatts, data = mlight)

# Partial dbRDA (controlling for site, if that’s a nuisance variable)
dbrda_partial <- dbrda(bray_dist ~ mean_mwatts + Condition(site), data = mlight)


#Assess model significance ---------------------------------------------------
anova(dbrda_full, by = "margin", permutations = 9999)     # marginal effects
anova(dbrda_full, permutations = 9999)                    # overall test

plot(dbrda_full, display = c("sites", "species"))
ordihull(dbrda_full, mlight$mean_mwatts, col = c("blue", "orange"), lwd = 2)
