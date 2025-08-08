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

# libraries-------------------------------------------------------------------------

# List of packages
packages <- c(
  "tidyverse", "magrittr", "lme4", "sjPlot", "ggeffects", "car",
  "glmmTMB", "corrplot", "effects", "reshape2", "DHARMa",
  "marginaleffects", "MuMIn", "performance", "viridis",
  "data.table", "janitor", "patchwork"
)

# Load all packages in one go
invisible(lapply(packages, library, character.only = TRUE))


load(file = "working_env/glmm_v3.RData")

# funtions ----------------------------------------------------------------

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


# Data --------------------------------------------------------------------

# Load the data
bm <- read_csv("data_for_analysis/prep_for_glmm/bm.csv") %>%
  clean_names()

# Check the structure of the data
str(bm)

# Check for missing values
colSums(is.na(bm))

# calculate week with the week function 

bm <- bm %>%
  mutate(
    wk = lubridate::week(noche)         # Extract week from noche
  )

filtered_bm <- bm %>%
  rename(sp = auto_id) %>%
  filter(!sp %in% c("Noise", "NoID")) # remove noise and noid from the analysis.

# miller activity index
bm_ai <- read_csv("data_for_analysis/prep_for_glmm/bm.miller.day.csv") %>% 
  clean_names() %>%
  filter(!sp %in% c("Noise", "NoID"))

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

# summarize moonlight by date but conditional moon_alt_degrees > 0

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



# light
light <- read_csv("data_for_analysis/lights/lightspectra_pioneer.csv") %>%
  clean_names() %>%
  filter(vert_horiz == "Horizontal") %>%                          # Keep only horizontal measures
  mutate(mwatts = rowMeans(across(c(watts_m1, watts_m2, watts_m3)), na.rm = TRUE)) %>%  # Mean of watts
  select(site, lux, yr, mwatts)    # Select relevant columns  

# Merge datasets ------------------------------------------------------------

# merge weather
# Merge crmo.wet.night into filtered_bm by matching dates
bm2 <- filtered_bm %>%
  left_join(crmo.wet.night, by = c("noche" = "date"))
summary(bm2)

# merge with moon
bm2 <- bm2 %>%
  left_join(moon_daily_avg, by = "noche") 

summary(bm2) #there is just one NA I can keep it that way. 

# merge with insects

bm2 <- bm2 %>%
  left_join(c_bugs, by = c("site", "wk", "yr")) # merge

# check for NAs
summary(bm2)  #lots of NAs in t.insect and t.lepidoptera.

# replace NAs in t.insect and t.lepidoptera with the mean values from c_bugs_mean


bm2 <- bm2 %>%
  left_join(c_bugs_mean, 
            by = c("yr", "treatmt.x" = "treatmt", "site"),
            suffix = c("", "_mean")) %>%  # rename mean columns directly
  mutate(
    t_insect = coalesce(t_insect, t_insect_mean),
    t_lepidoptera = coalesce(t_lepidoptera, t_lepidoptera_mean)
  ) %>%
  select(-t_insect_mean, -t_lepidoptera_mean)


# merge light 

bm2<- bm2 %>%
  left_join(light, by = c("site", "yr"))


# correlation -------------------------------------------------------------


# 1. Select numeric columns, optionally drop unique id/group columns
numeric_data <- bm2 %>% 
  select(where(is.numeric)) %>% 
  select(-any_of(c("yr", "trmt_bin"))) # add/remove columns as needed

# 2. Remove zero-variance columns
numeric_data <- numeric_data %>% select(where(~sd(., na.rm = TRUE) > 0))

# 3. Compute correlation matrix
cor_mat <- cor(numeric_data, use = "pairwise.complete.obs")

# 4. Visualize correlation matrix
corrplot(cor_mat, order = 'AOE')


# standardize -------------------------------------------------------------

#calculate Julian day from 'noche' and scale variables
bm2<-bm2 %>% 
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
  "lux",
  "mwatts",
  "jday"
)

bm2 <- bm2 %>%
  scale_by_2sd_tidy(variables_to_scale)

summary(bm2)

# year standardize. 

# make year between -1:1
bm2 <- bm2 %>%
  mutate(yr_s = case_when(
    yr == 2021 ~ -1,
    yr == 2022 ~ 0,
    TRUE ~ 1
  ))


# sopecies names for graphs

species <- data.frame(
  sp = c("ANTPAL", "CORTOW", "EPTFUS", "EUDMAC", "LASCIN", "LASNOC",
         "MYOCAL", "MYOCIL", "MYOEVO", "MYOLUC", "MYOTHY", "MYOVOL",
         "MYOYUM", "PARHES"),
  species_name = c("Antrozous pallidus", "Corynorhinus townsendii", "Eptesicus fuscus", "Euderma maculatum",
                   "Lasiurus cinereus", "Lasiurus noctivagans", "Myotis californicus", "Myotis ciliolabrum",
                   "Myotis evotis", "Myotis lucifugus", "Myotis thysanodes", "Myotis volans",
                   "Myotis yumanensis", "Parastrellus hesperus")
)

species <- species %>%
  mutate(
    genus = word(species_name, 1),
    species = word(species_name, 2),
    sp_label = paste0(substr(genus, 1, 1), ".", species)
  )



# modeling ----------------------------------------------------------------

# model with light spectra instead of binomial 

m1.1nb <- glmmTMB(
  #fixed effects
  n ~ lux_s + jday_s + I(jday_s^2) + avg_moonlight_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + lux_s + jday_s + I(jday_s^2) | sp) +
    #interactions
    jday_s * lux_s + I(jday_s^2) * lux_s + yr_s * lux_s + lux_s*t_lepidoptera_s + lux_s*avg_moonlight_s,
  data = bm2,
  family = nbinom2(link = "log")
)
summary(m1.1nb)
check_singularity(m1.1nb)
check_zeroinflation(m1.1nb)
check_overdispersion(m1.1nb)
r2(m1.1nb)
DHARMa::simulateResiduals(m1.1nb, plot = TRUE, quantreg = TRUE)


# model with watts

m1.2nb <- glmmTMB(
  #fixed effects
  n ~ mwatts_s + jday_s + I(jday_s^2) + avg_moonlight_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + mwatts_s + jday_s + I(jday_s^2) | sp) +
    #interactions
    jday_s * mwatts_s + I(jday_s^2) * mwatts_s + yr_s * mwatts_s + mwatts_s*t_lepidoptera_s + mwatts_s*avg_moonlight_s,
  data = bm2,
  family = nbinom2(link = "log")
)

summary(m1.2nb)
check_singularity(m1.2nb)
check_zeroinflation(m1.2nb)
check_overdispersion(m1.2nb)
r2(m1.2nb)
DHARMa::simulateResiduals(m1.2nb, plot = TRUE, quantreg = TRUE)
anova(m1.1nb,m1.2nb) # the model with light spectra is not better than the model with the insects and ear/arm ratio.)


m1.3 <- glmmTMB(
  #fixed effects
  n ~ mwatts_s + jday_s + I(jday_s^2) + avg_moonlight_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + mwatts_s + jday_s + I(jday_s^2) | sp) +
    #interactions
    jday_s * mwatts_s + I(jday_s^2) * mwatts_s + yr_s * mwatts_s + mwatts_s*t_lepidoptera_s + mwatts_s*avg_moonlight_s,
  data = bm2,
  family = poisson(link = "log")
)

summary(m1.3)
check_singularity(m1.3)
check_zeroinflation(m1.3)
check_overdispersion(m1.3)
r2(m1.3)

# marginal effects --------------------------------------------------------

# ------------------------------
# Marginal Effects of Light (Watts) on Bat Calls ------
# ------------------------------

# Standardization parameters (Gelman style: x / (2 * sd))
sdwatts <- sd(bm2$mwatts, na.rm = TRUE)
mwatts_mean <- mean(bm2$mwatts, na.rm = TRUE)

# ------------------------------
# 1. COMMUNITY-LEVEL PREDICTIONS (random effects removed)
# ------------------------------
pred_watts <- predictions(
  m1.2nb,
  newdata = datagrid(
    sp = NA,                     
    site = NA,                    
    mwatts_s = seq(min(bm2$mwatts_s, na.rm = TRUE),
                   max(bm2$mwatts_s, na.rm = TRUE),
                   length.out = 100)
  ),
  re.form = NA
) %>%
  mutate(
    mwatts = mwatts_s * (2 * sdwatts) + mwatts_mean
  )

# Plot community-level effect
p_comm <- ggplot(pred_watts, aes(x = mwatts, y = estimate)) +
  geom_line(linewidth = 1, color = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "grey70") +
  labs(
    x = "Light (watts)",
    y = "Predicted bat calls (n)",
    title = "Community-level effect of light on bat calls"
  ) +
  theme_minimal()

# ------------------------------
# 2. SPECIES-SPECIFIC PREDICTIONS (random slopes for sp retained)
# ------------------------------
pred_species <- predictions(
  m1.2nb,
  newdata = datagrid(
    sp = unique(bm2$sp),
    site = NA,
    mwatts_s = seq(min(bm2$mwatts_s, na.rm = TRUE),
                   max(bm2$mwatts_s, na.rm = TRUE),
                   length.out = 100)
  ),
  re.form = NULL
) %>%
  mutate(
    mwatts = mwatts_s * (2 * sdwatts) + mwatts_mean
  )

# Plot species-specific effects (facet per species)
p_species <- ggplot(pred_species, aes(x = mwatts, y = estimate)) +
  geom_line(color = "steelblue") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") +
  facet_wrap(~ sp, scales = "free_y") +
  labs(
    x = "Light (watts)",
    y = "Predicted bat calls (n)",
    title = "Species-specific effect of light on bat calls"
  ) +
  theme_minimal()

# ------------------------------
# 3. Display both plots
# ------------------------------
p_comm / p_species


# ------------------------------
# Marginal effect mwatts and year---------------------
# ------------------------------
pred_w.y <- predictions(
  m1.2nb,
  newdata = datagrid(
    yr_s = c(-1, 0, 1),  # standardized values for 2021, 2022, 2023
    mwatts_s = seq(
      min(bm2$mwatts_s, na.rm = TRUE),
      max(bm2$mwatts_s, na.rm = TRUE),
      length.out = 100
    )
  ),
  re.form = NA  # community-level predictions
) %>%
  mutate(
    year = case_when(
      yr_s == -1 ~ 2021,
      yr_s == 0  ~ 2022,
      yr_s == 1  ~ 2023
    )
  )

# Plot
ggplot(pred_w.y, aes(x = mwatts_s, y = estimate, color = factor(year))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = factor(year)),
              alpha = 0.2, color = NA) +
  labs(
    x = "Standardized watts",
    y = "Predicted bat calls (n)",
    color = "Year",
    fill  = "Year",
    title = "Interaction effect: watts × year"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ------------------------------
# Marginal effect Jday-----
# ------------------------------
sdjday <- sd(bm2$jday, na.rm = TRUE)
jday_mean <- mean(bm2$jday, na.rm = TRUE)

# ------------------------------
# COMMUNITY-LEVEL PREDICTIONS for Jday ----------------------
# ------------------------------
pred_jday <- predictions(
  m1.2nb,
  newdata = datagrid(
    sp = NA,            # remove random slope variation
    site = NA,
    jday_s = seq(min(bm2$jday_s, na.rm = TRUE),
                 max(bm2$jday_s, na.rm = TRUE),
                 length.out = 100),
    mwatts_s = mean(bm2$mwatts_s, na.rm = TRUE),
    avg_moonlight_s = mean(bm2$avg_moonlight_s, na.rm = TRUE),
    nit_avg_wspm_s_s = mean(bm2$nit_avg_wspm_s_s, na.rm = TRUE),
    yr_s = mean(bm2$yr_s, na.rm = TRUE),
    t_lepidoptera_s = mean(bm2$t_lepidoptera_s, na.rm = TRUE)
  ),
  re.form = NA
) %>%
  mutate(
    jday = jday_s * (2 * sdjday) + jday_mean
  )

# ------------------------------
# Plot
# ------------------------------
ggplot(pred_jday, aes(x = jday, y = estimate)) +
  geom_line(linewidth = 1, color = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "grey70", alpha = 0.3) +
  labs(
    x = "Julian day",
    y = "Predicted bat calls (n)",
    title = "Community-level effect of Julian day on bat calls"
  ) +
  theme_minimal()



# ------------------------------
# SPECIES-SPECIFIC PREDICTIONS for Jday ---------------------
# ------------------------------
pred_jday_sp <- predictions(
  m1.2nb,
  newdata = datagrid(
    sp = unique(bm2$sp),  # keep each species
    site = NA,            # no site-level variation
    jday_s = seq(min(bm2$jday_s, na.rm = TRUE),
                 max(bm2$jday_s, na.rm = TRUE),
                 length.out = 100),
    mwatts_s = mean(bm2$mwatts_s, na.rm = TRUE),
    avg_moonlight_s = mean(bm2$avg_moonlight_s, na.rm = TRUE),
    nit_avg_wspm_s_s = mean(bm2$nit_avg_wspm_s_s, na.rm = TRUE),
    yr_s = mean(bm2$yr_s, na.rm = TRUE),
    t_lepidoptera_s = mean(bm2$t_lepidoptera_s, na.rm = TRUE)
  ),
  re.form = NULL  # keep species-specific random effects
) %>%
  mutate(
    jday = jday_s * (2 * sdjday) + jday_mean
  )

# ------------------------------
# Plot
# ------------------------------
ggplot(pred_jday_sp, aes(x = jday, y = estimate, color = sp, fill = sp)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  labs(
    x = "Julian day",
    y = "Predicted bat calls (n)",
    title = "Species-specific phenology curves"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggplot(pred_jday_sp, aes(x = jday, y = estimate, color = sp, fill = sp)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  labs(
    x = "Julian day",
    y = "Predicted bat calls (n)",
    title = "Species-specific phenology curves"
  ) +
  facet_wrap(~ sp, scales = "free_y")  +
  theme_minimal() +
  theme(
    legend.position = "none",   # often clearer in facet plots
    strip.text.x = element_text(size = 8) # smaller labels if many species
  )


# ------------------------------
# Marginal effect for moonlight-------------------------
# ------------------------------
sdmoonlight <- sd(bm2$jday, na.rm = TRUE)
moonlight_mean <- mean(bm2$jday, na.rm = TRUE)

# ------------------------------
# COMMUNITY-LEVEL PREDICTIONS for moonlight---------------------------
# ------------------------------

pred_moon<- predictions(
  m1.2nb,
  newdata = datagrid(
    sp=NA,
    site= NA, 
    avg_moonlight_s = seq(
      min(bm2$avg_moonlight_s, na.rm=T),
      max(bm2$avg_moonlight_s, na.rm=T),
      length.out= 200
    )
    ),
    re.form= NA # this specify population=community of sp level predictions
  ) %>% 
    mutate(
      moon= avg_moonlight_s * (2*sdmoonlight) + moonlight_mean
    )

# Plot
ggplot(pred_moon, aes(x = moon, y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(
    x = "Average moonlight",
    y = "Predicted bat calls (n)",
    title = "Community-level effect of moonlight"
  ) +
  theme_minimal()


# ------------------------------
# Marginal effect for nit_avg_wspm_s_s-------------------------
# ------------------------------

sdwind <- sd(bm2$nit_avg_wspm_s, na.rm = TRUE)
wind_mean <- mean(bm2$nit_avg_wspm_s, na.rm = TRUE)


# COMMUNITY-LEVEL PREDICTIONS for moonlight---------------------------

pred_wind<-predictions(
  m1.2nb,
  newdata = datagrid(
    nit_avg_wspm_s_s = seq(
      min(bm2$nit_avg_wspm_s_s, na.rm=T),
      max(bm2$nit_avg_wspm_s_s, na.rm = T),
      length.out= 200
    )
  ),
  re.form=NA
) %>% 
  mutate(
    wind= nit_avg_wspm_s_s * (2*sdwind) + wind_mean

  )


ggplot(pred_wind, aes(x = wind, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .5, fill = "skyblue", color = NA) +
  geom_line(color = "blue", linewidth = 1) +
  labs(
    x = "Wind speed (m/s)",
    y = "Predicted bat calls (n)",
    title = "Community-level effect of wind speed"
  ) +
  theme_minimal()


# save working environment ------------------------------------------------

save.image("working_env/glmm_v3.RData")
