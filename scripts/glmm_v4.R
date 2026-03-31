#=============================================================================
# Script Name:    glmm_v4.R
# Purpose:        Fit GLMMs for bat paper using weather, and moon data, and light spectra.
#
# Version:        v4 - Includes: we are going back and build the model using light as -1/1 
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


# load(file = "working_env/glmm_v3.RData")

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


# Define a function to calculate and print c-hat
calculate_c_hat <- function(model) {
  residual_deviance <- deviance(model)
  residual_df <- df.residual(model)
  c_hat_deviance <- residual_deviance / residual_df
  print(c_hat_deviance)
}


# Data --------------------------------------------------------------------

# Load the data
bm <- read_csv("data_for_analysis/prep_for_glmm/bm.csv") %>%
  clean_names()

glimpse(bm) # check the structure of the data

# Check for missing values
colSums(is.na(bm))

# calculate week with the week function 

bm <- bm %>%
  mutate(
    wk = lubridate::week(noche)         # Extract week from noche
  )

# here we remove noise and noids. 

filtered_bm <- bm %>%
  rename(sp = sp) %>%
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

species <- species %>%
  mutate(
    genus = word(species_name, 1),
    species = word(species_name, 2),
    sp_label = paste0(substr(genus, 1, 1), ".", species)
  )

# we make treatment bin -1 for dark and 1 for lit. 
bm2 <- bm2 %>%
  mutate(trmt_bin = if_else(treatmt.x == "lit", 1, -1)) %>%
  left_join(species %>% select(sp, sp_label), by = "sp") # add species labels for plotting

glimpse(bm2)

# now we merge the bm_ai with the bm2 data to have both the Miller activity index and the echolocation count data.

bm2<-bm_ai %>%
  select(site,noche, sp, activity_min) %>%
  left_join(bm2, by = c("noche", "sp", "site")) # merge by noche and spm

# modeling ----------------------------------------------------------------

# model with treatment as binomial 

# first try poisson and just the data summary. 

m0<-glmmTMB(n ~ trmt_bin + (1 | site),
           data = bm2, family = poisson(link = "log"))

# check model
check_singularity(m0)
check_zeroinflation(m0)
calculate_c_hat(m0)  # indicates the is overdispersion we need to try negative binomial. 
performance_mae(m0)
range(bm2$n)
hist(bm2$n, breaks = 500)
performance::r2(m0)
DHARMa::simulateResiduals(m0, n = 1000, plot = TRUE)
summary(m0)


# model with a poisson distribution is overdispersed c-hat=81 so we try negative binomial.

m1<- glmmTMB(
  #fixed effects
  n ~ trmt_bin + (1 | site),
  data = bm2,
  family = nbinom2(link = "log"))

# check model
check_singularity(m1)
check_zeroinflation(m1)
calculate_c_hat(m1)  # indicates the is overdispersion we need to try negative binomial. 
performance_mae(m1)
range(bm2$n)
hist(bm2$n, breaks = 500)
performance::r2(m1)
DHARMa::simulateResiduals(m1, n = 1000, plot = TRUE)
summary(m1)

anova(m0, m1) # the model with negative binomial is better than the poisson.

# model with a negative binomial has a better c-hat value of 1.2 although still overdisperssed
# now we try adding seasonality as fixed effects. R conidtional is around 0.24 and marginal 0.155

# now we add seasonality to the model 

m1.1 <- glmmTMB(
  #fixed effects
  n ~ trmt_bin + jday_s + I(jday_s^2) +
    #random effects
    (1 | site),  
    #interactions
  data = bm2,
  family = nbinom2(link = "log")
)


check_singularity(m1.1)
m1.1$sdr$pdHess
m1.1$fit$message
check_zeroinflation(m1.1)
calculate_c_hat(m1.1)  # indicates the is overdispersion we need to try negative binomial. 
performance_mae(m1.1)
range(bm2$n)
hist(bm2$n, breaks = 500)
performance::r2(m1.1)
DHARMa::simulateResiduals(m1.1, n = 1000, plot = TRUE)
anova(m1, m1.1) 
summary(m1.1)


# it seems like the model with seasonality is better thatn the without it.
# now we try the environmental predictors.

m1.2<-glmmTMB(
  #fixed effects
  n ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site),
  data = bm2,
  family = nbinom2(link = "log")
)

check_singularity(m1.2)
m1.2$sdr$pdHess
m1.2$fit$message
check_zeroinflation(m1.2)
calculate_c_hat(m1.2)  # indicates the is overdispersion we need to try negative binomial. 
performance_mae(m1.2)
range(bm2$n)
hist(bm2$n, breaks = 500)
performance::r2(m1.2)
DHARMa::simulateResiduals(m1.2, n = 1000, plot = TRUE)
anova(m1, m1.1, m1.2) # I have an issue because there's a row that has a missing value for moon data. 
summary(m1.2)


# now lets try the model treatment inside random slope 

m1.3<- glmmTMB(
  #fixed effects
  n ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin | sp),  
  data = bm2,
  family = nbinom2(link = "log")
)

check_singularity(m1.3)
m1.3$sdr$pdHess
m1.3$fit$message
check_zeroinflation(m1.3)
calculate_c_hat(m1.3) 
performance_mae(m1.3)
range(bm2$n)
hist(bm2$n, breaks = 500)
performance::r2(m1.3)
DHARMa::simulateResiduals(m1.3, n = 1000, plot = TRUE)
anova( m1.2, m1.3)
summary(m1.3)


#now let's try seasonality inside the random slopes

m1.4<-glmmTMB(
  #fixed effects
  n ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp),  
  data = bm2,
  family = nbinom2(link = "log")
)

check_singularity(m1.4)
m1.4$sdr$pdHess
m1.4$fit$message
check_zeroinflation(m1.4)
calculate_c_hat(m1.4) 
performance_mae(m1.4)
range(bm2$n)
hist(bm2$n, breaks = 500)
performance::r2(m1.4)
DHARMa::simulateResiduals(m1.4, n = 1000, plot = TRUE)
anova( m1.2, m1.3,m1.4)
summary(m1.4)


# now let's try for interactions of treatment with seasonalty 

m1.5<-glmmTMB(
  #fixed effects
  n ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp) +
    #interactions
    jday_s * trmt_bin + I(jday_s^2) * trmt_bin,  
  data = bm2,
  family = nbinom2(link = "log")
)


check_singularity(m1.5)
m1.5$sdr$pdHess
m1.5$fit$message
check_zeroinflation(m1.5)
calculate_c_hat(m1.5) 
performance_mae(m1.5)
range(bm2$n)
hist(bm2$n, breaks = 500)
performance::r2(m1.5)
DHARMa::simulateResiduals(m1.5, n = 1000, plot = TRUE)
performance::check_collinearity(m1.5)
anova( m1.2, m1.3,m1.4,m1.5)
summary(m1.5)

# adding the treatment and light interaction we improve the fit but not as much as adding the seasonality inside the random slopes. 

# now we try to add the moon and treatment interaction 

m1.6<-glmmTMB(
  #fixed effects
  n ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp) +
    #interactions
    jday_s * trmt_bin + I(jday_s^2) * trmt_bin + trmt_bin*avg_moonlight_s,  
  data = bm2,
  family = nbinom2(link = "log")
)


check_singularity(m1.6)
m1.6$sdr$pdHess
m1.6$fit$message
check_zeroinflation(m1.6)
calculate_c_hat(m1.6) 
performance_mae(m1.6)
range(bm2$n)
hist(bm2$n, breaks = 500)
performance::r2(m1.6)
DHARMa::simulateResiduals(m1.6, n = 1000, plot = TRUE)
performance::check_collinearity(m1.6)
performance::check_collinearity(m1.6)
anova( m1.2, m1.3,m1.4,m1.5,m1.6)
summary(m1.6)


#now we try the interactrion with year

m1.8<-glmmTMB(
  #fixed effects
  n ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp) +
    #interactions
    jday_s * trmt_bin + I(jday_s^2) * trmt_bin + trmt_bin*avg_moonlight_s + trmt_bin*yr_s,  
  data = bm2,
  family = nbinom2(link = "log")
)

check_singularity(m1.8)
m1.8$sdr$pdHess
m1.8$fit$message
check_zeroinflation(m1.8)
calculate_c_hat(m1.8) 
performance_mae(m1.8)
range(bm2$n)
hist(bm2$n, breaks = 500)
performance::r2(m1.8)
DHARMa::simulateResiduals(m1.8, n = 1000, plot = TRUE)
performance::check_collinearity(m1.8)
anova( m1.2, m1.3,m1.4,m1.5,m1.6, m1.8)
summary(m1.8)

gtsummary::tbl_regression(m1.8, exponentiate = TRUE)

# this model seems strong and it show ecologically interesting effects light the treatment cahnging throuhg the years. plus this effect improves the model chisqu 30.6651   pr(chisqu) 3.066e-08 ***


# lets try interaction with lepidptera

m1.9<-glmmTMB(
  #fixed effects
  n ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp) +
    #interactions
    jday_s * trmt_bin + I(jday_s^2) * trmt_bin + trmt_bin*avg_moonlight_s + trmt_bin*yr_s + trmt_bin*t_lepidoptera_s,  
  data = bm2,
  family = nbinom2(link = "log")
)

check_singularity(m1.9)
m1.9$sdr$pdHess
m1.9$fit$message
check_zeroinflation(m1.9)
calculate_c_hat(m1.9) 
performance_mae(m1.9)
range(bm2$n)
hist(bm2$n, breaks = 500)
performance::r2(m1.9)
DHARMa::simulateResiduals(m1.9, n = 1000, plot = TRUE)
performance::check_collinearity(m1.9)
anova( m1.2, m1.3,m1.4,m1.5,m1.6, m1.8, m1.9)
summary(m1.9)

# including the lepidoptera does not improve the model so we keep m1.8


# runa a model with all the same predictors but with the miller activity index as the response variable.

m1.10<- glmmTMB(
  #fixed effects
  activity_min ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp) +
    #interactions
    jday_s * trmt_bin + I(jday_s^2) * trmt_bin + trmt_bin*avg_moonlight_s + trmt_bin*yr_s,  
  data = bm2,
  family = nbinom2(link = "log")
)

check_singularity(m1.10)
m1.10$sdr$pdHess
m1.10$fit$message
check_zeroinflation(m1.10)
calculate_c_hat(m1.10) 
performance_mae(m1.10)
range(bm2$n)
hist(bm2$n, breaks = 500)
performance::r2(m1.10)
DHARMa::simulateResiduals(m1.10, n = 1000, plot = TRUE)
performance::check_collinearity(m1.10)
anova( m1.2, m1.3,m1.4,m1.5,m1.6, m1.8, m1.10)
summary(m1.10)

# marginal effects m1.8 ---------------------------------------------------


typical <- list(
  jday_s = 0,
  avg_moonlight_s = 0,
  nit_avg_temp_c_s = 0,
  nit_avg_wspm_s_s = 0,
  yr_s = 0,
  t_lepidoptera_s = 0
)


pred_trmt <- predictions(
  m1.8,
  newdata = datagrid(
    trmt_bin = c(-1, 1),
    sp = unique(bm2$sp),
    jday_s = typical$jday_s,
    avg_moonlight_s = typical$avg_moonlight_s,
    nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
    nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
    yr_s = typical$yr_s,
    t_lepidoptera_s = typical$t_lepidoptera_s
  ),
  type = "response"   # VERY IMPORTANT → gives counts, not log scale
)



pred_trmt <- pred_trmt %>%
  mutate(
    sp_label = str_replace(sp, "([A-Z])[A-Z]+([A-Z]+)", "\\1. \\2"),
    trmt = ifelse(trmt_bin == -1, "Dark", "Lit")
  )

pred_trmt <- pred_trmt %>%
  group_by(sp) %>%
  mutate(
    change = estimate[trmt_bin == -1] - estimate[trmt_bin == 1]
  )

library(ggplot2)

p_trmt <- ggplot(pred_trmt, aes(x = trmt_bin, y = estimate, group = sp)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.1
  ) +
  facet_wrap(~ sp, scales = "free_y") +
  labs(
    x = "Treatment",
    y = "Predicted bat calls",
    title = "Species-specific response of bat activity to artificial light"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_trmt

bm_plot <- bm2 %>%
  mutate(
    trmt = ifelse(trmt_bin == -1, "Dark", "Lit"),
    sp_label = str_replace(sp, "([A-Z])[A-Z]+([A-Z]+)", "\\1. \\2")
  )

# improved graph with raw data in the bakground and log y x to better visualize the differences between species.

p_trmt <- ggplot() +
  
  # 🔹 RAW DATA (background)
  geom_jitter(
    data = bm_plot,
    aes(x = trmt_bin, y = n),
    width = 0.15,
    alpha = 0.15,
    size = 0.6,
    color = "grey40"
  ) +
  
  # 🔹 MODEL PREDICTIONS (foreground)
  geom_point(
    data = pred_trmt,
    aes(x = trmt_bin, y = estimate),
    size = 2
  ) +
  
  geom_line(
    data = pred_trmt,
    aes(x = trmt_bin, y = estimate, group = sp),
    linewidth = 0.8
  ) +
  
  facet_wrap(~ sp, scales = "free_y") +
  
  labs(
    x = "Treatment",
    y = "Predicted bat calls",
    title = "Species-specific response of bat activity to artificial light"
  ) +
  scale_y_continuous(trans = "log1p")+
  scale_x_continuous(
    breaks = c(-1, 1),
    labels = c("Dark", "Lit")
  ) +
  
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_trmt

ggsave(
  filename = "figures/glmm_v4/trmt_v1.png",
  plot = p_trmt,
  width = 10,
  height = 8,
  dpi = 300
)




# jday 


jday_seq <- seq(
  min(bm2$jday_s, na.rm = TRUE),
  max(bm2$jday_s, na.rm = TRUE),
  length.out = 100
)

pred_jday <- predictions(
  m1.8,
  newdata = datagrid(
    jday_s = jday_seq,
    trmt_bin = c(-1, 1),
    sp = unique(bm2$sp),
    avg_moonlight_s = 0,
    nit_avg_temp_c_s = 0,
    nit_avg_wspm_s_s = 0,
    yr_s = 0,
    t_lepidoptera_s = 0
  ),
  type = "response"
)

jday_mean <- mean(bm2$jday, na.rm = TRUE)
jday_sd   <- sd(bm2$jday, na.rm = TRUE)

pred_jday <- pred_jday %>%
  mutate(
    jday = jday_s * (2 * jday_sd) + jday_mean,
    trmt = ifelse(trmt_bin == -1, "Dark", "Lit"),
    sp_label = str_replace(sp, "([A-Z])[A-Z]+([A-Z]+)", "\\1. \\2")
  )

bm_jday_plot <- bm2 %>%
  mutate(
    trmt = ifelse(trmt_bin == -1, "Dark", "Lit"),
    sp_label = str_replace(sp, "([A-Z])[A-Z]+([A-Z]+)", "\\1. \\2")
  )


p_jday <- ggplot() +
  geom_point(
    data = bm_jday_plot,
    aes(x = jday, y = n, color = trmt),
    alpha = 0.10,
    size = 0.7
  ) +
  geom_line(
    data = pred_jday,
    aes(x = jday, y = estimate, color = trmt),
    linewidth = 1.2
  ) +
  facet_wrap(~ sp, scales = "free_y") +
  labs(
    x = "Julian day",
    y = "Predicted bat calls",
    color = "Treatment",
    title = "Species-specific seasonal patterns in bat activity"
  ) +
  scale_color_manual(
    values = c(
      "Lit" = "grey70",  # light grey
      "Dark"  = "grey20"   # darker grey
    )
  )+
  scale_y_continuous(trans = "log1p")+
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "top"
  )

p_jday


ggsave(
  filename = "figures/glmm_v4/jday_v1.png",
  plot = p_jday,
  width = 10,
  height = 8,
  dpi = 300
)


# marginal effects of light and moonlight interactions.


moon_seq <- seq(
  min(bm2$avg_moonlight_s, na.rm = TRUE),
  max(bm2$avg_moonlight_s, na.rm = TRUE),
  length.out = 100
)

pred_moon <- predictions(
  m1.8,
  newdata = datagrid(
    avg_moonlight_s = moon_seq,
    trmt_bin = c(-1, 1),
    sp = unique(bm2$sp),
    jday_s = 0,
    nit_avg_temp_c_s = 0,
    nit_avg_wspm_s_s = 0,
    yr_s = 0,
    t_lepidoptera_s = 0
  ),
  type = "response"
) %>%
  mutate(
    trmt = ifelse(trmt_bin == -1, "Dark", "Lit"),
    sp_label = str_replace(sp, "([A-Z])[A-Z]+([A-Z]+)", "\\1. \\2")
  )


p_moon <- ggplot() +
  # geom_point(
  #   data = bm2,
  #   aes(x = avg_moonlight_s, y = n, color = treatmt.x),
  #   alpha = 0.08,
  #   size = 0.6
  # ) +
  geom_line(
    data = pred_moon,
    aes(x = avg_moonlight_s, y = estimate, color = trmt),
    linewidth = 1
  ) +
  facet_wrap(~ sp, scales = "free_y") +
  labs(
    x = "Moonlight (scaled)",
    y = "Predicted bat calls",
    color = "Treatment"
  ) +
  scale_color_manual(values = c("Lit" = "grey70", "Dark" = "grey20")) +
  scale_y_continuous(trans = "log1p") +
  theme_minimal()

p_moon


ggsave(
  filename = "figures/glmm_v4/moon_v1.png",
  plot = p_moon,
  width = 10,
  height = 8,
  dpi = 300
)

# now the year interaction 


pred_year <- predictions(
  m1.8,
  newdata = datagrid(
    yr_s = c(-1, 0, 1),   # ← include ALL years
    trmt_bin = c(-1, 1),
    sp = unique(bm2$sp),
    jday_s = 0,
    avg_moonlight_s = 0,
    nit_avg_temp_c_s = 0,
    nit_avg_wspm_s_s = 0,
    t_lepidoptera_s = 0
  ),
  type = "response"
)

pred_year <- pred_year %>%
  mutate(
    trmt = ifelse(trmt_bin == -1, "Dark", "Lit"),
    year = case_when(
      yr_s == -1 ~ "2021",
      yr_s ==  0 ~ "2022",
      yr_s ==  1 ~ "2023"
    ),
    sp_label = str_replace(sp, "([A-Z])[A-Z]+([A-Z]+)", "\\1. \\2")
  )


p_year <- ggplot(pred_year,
                 aes(x = trmt, y = estimate, group = year, color = year)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  facet_wrap(~ sp, scales = "free_y") +
  labs(
    x = "Treatment",
    y = "Predicted bat calls",
    color = "Year",
    title = "Year-specific treatment effects on bat activity"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "top"
  )

p_year
# these panels might need to be organized the colors need to change and the olrder of the panesl too. I want to present the panels in three groups first the ones that constantly declined like antpal or edumac then the ones that fip their response like cortow dna eptfus 
# and then the ones that show increases in lit sites constantly like lasnoc or myluc should I add data points?


p_year_alt <- ggplot(pred_year,
                     aes(x = year, y = estimate, color = trmt, group = trmt)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  facet_wrap(~ sp, scales = "free_y") +
  labs(
    x = "Year",
    y = "Predicted bat calls",
    color = "Treatment",
    title = "Treatment × Year interaction in bat activity"
  ) +
  scale_color_manual(values = c("Lit" = "grey70", "Dark" = "grey20")) +
  theme_minimal()

p_year_alt


bm_year_plot <- bm2 %>%
  mutate(
    year = as.character(yr),
    trmt = ifelse(trmt_bin == -1, "Dark", "Lit"),
    sp_label = str_replace(sp, "([A-Z])[A-Z]+([A-Z]+)", "\\1. \\2")
  )

p_year_alt <- ggplot() +
  
  # 🔹 raw data (background)
  # geom_jitter(
  #   data = bm_year_plot,
  #   aes(x = year, y = n, color = trmt),
  #   width = 0.15,
  #   alpha = 0.08,
  #   size = 0.6
  # ) +
  
  # 🔹 model predictions
  geom_point(
    data = pred_year,
    aes(x = year, y = estimate, color = trmt),
    size = 2
  ) +
  geom_line(
    data = pred_year,
    aes(x = year, y = estimate, color = trmt, group = trmt),
    linewidth = 1
  ) +
  # geom_errorbar(
  #   data = pred_year,
  #   aes(x = year, ymin = conf.low, ymax = conf.high, color = trmt),
  #   width = 0.1,
  #   position = position_dodge(width = 0.2)
  # )+
  
  facet_wrap(~ sp, scales = "free_y") +
  
  labs(
    x = "Year",
    y = "Predicted bat calls",
    color = "Treatment",
    title = "Treatment × Year interaction in bat activity"
  ) +
  
  scale_color_manual(
    values = c("Lit" = "grey70", "Dark" = "grey20")
  ) +
  
  scale_y_continuous(trans = "log1p") +
  
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "top"
  )

p_year_alt
ggsave(
  filename = "figures/glmm_v4/year_v1.png",
  plot = p_year_alt,
  width = 10,
  height = 8,
  dpi = 300
)

summary(m1.8)

# not graphs for the temperature 


temp_seq <- seq(
  min(bm2$nit_avg_temp_c_s, na.rm = TRUE),
  max(bm2$nit_avg_temp_c_s, na.rm = TRUE),
  length.out = 100
)

pred_temp_comm <- predictions(
  m1.8,
  newdata = datagrid(
    nit_avg_temp_c_s = temp_seq,
    trmt_bin = 0,
    jday_s = 0,
    avg_moonlight_s = 0,
    nit_avg_wspm_s_s = 0,
    yr_s = 0,
    t_lepidoptera_s = 0
  ),
  type = "response",
  re.formula = NA
)

temp_mean <- mean(bm2$nit_avg_temp_c, na.rm = TRUE)
temp_sd   <- sd(bm2$nit_avg_temp_c, na.rm = TRUE)

pred_temp_comm <- pred_temp_comm %>%
  mutate(
    temp_c = nit_avg_temp_c_s * (2 * temp_sd) + temp_mean
  )


p_temp_comm <- ggplot(pred_temp_comm, aes(x = temp_c, y = estimate)) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2,
    fill = "grey70",
    color = NA
  ) +
  geom_line(
    linewidth = 1,
    color = "grey20"
  ) +
  labs(
    x = "Nightly average temperature (°C)",
    y = "Predicted bat call counts",
    title = "Community-level effect of temperature on bat activity"
  ) +
  theme_minimal()

p_temp_comm

ggsave(
  filename = "figures/glmm_v4/temp_comm_v1.png",
  plot = p_temp_comm,
  width = 6,
  height = 4,
  dpi = 300
)

# not lets do the wind 

wind_seq <- seq(
  min(bm2$nit_avg_wspm_s_s, na.rm = TRUE),
  max(bm2$nit_avg_wspm_s_s, na.rm = TRUE),
  length.out = 100
)

pred_wind_comm <- predictions(
  m1.8,
  newdata = datagrid(
    nit_avg_wspm_s_s = wind_seq,
    trmt_bin = 0,
    jday_s = 0,
    avg_moonlight_s = 0,
    nit_avg_wspm_s_s = 0,
    yr_s = 0,
    t_lepidoptera_s = 0
  ),
  type = "response",
  re.formula = NA
)

wind_mean <- mean(bm2$nit_avg_wspm_s, na.rm = TRUE)
wind_sd   <- sd(bm2$nit_avg_wspm_s, na.rm = TRUE)

pred_wind_comm <- pred_wind_comm %>%
  mutate(
    windms = nit_avg_wspm_s_s * (2 * wind_sd) + wind_mean
  )


p_wind_comm <- ggplot(pred_wind_comm, aes(x = windms, y = estimate)) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2,
    fill = "grey70",
    color = NA
  ) +
  geom_line(
    linewidth = 1,
    color = "grey20"
  ) +
  labs(
    x = "Nightly average wind (m/s)",
    y = "Predicted bat call counts",
    title = "Community-level effect of temperature on bat activity"
  ) +
  theme_minimal()

p_wind_comm

ggsave(
  filename = "figures/glmm_v4/wind_comm_v1.png",
  plot = p_wind_comm,
  width = 6,
  height = 4,
  dpi = 300
)




# marginal effects m1.10 acoustic index -----------------------------------

typical <- list(
  jday_s = 0,
  avg_moonlight_s = 0,
  nit_avg_temp_c_s = 0,
  nit_avg_wspm_s_s = 0,
  yr_s = 0,
  t_lepidoptera_s = 0
)


pred_trmt <- predictions(
  m1.10,
  newdata = datagrid(
    trmt_bin = c(-1, 1),
    sp = unique(bm2$sp),
    jday_s = typical$jday_s,
    avg_moonlight_s = typical$avg_moonlight_s,
    nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
    nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
    yr_s = typical$yr_s,
    t_lepidoptera_s = typical$t_lepidoptera_s
  ),
  type = "response"   # VERY IMPORTANT → gives counts, not log scale
)



pred_trmt <- pred_trmt %>%
  mutate(
    sp_label = str_replace(sp, "([A-Z])[A-Z]+([A-Z]+)", "\\1. \\2"),
    trmt = ifelse(trmt_bin == -1, "Dark", "Lit")
  )

pred_trmt <- pred_trmt %>%
  group_by(sp) %>%
  mutate(
    change = estimate[trmt_bin == -1] - estimate[trmt_bin == 1]
  )

library(ggplot2)

p_trmt <- ggplot(pred_trmt, aes(x = trmt_bin, y = estimate, group = sp)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.1
  ) +
  facet_wrap(~ sp, scales = "free_y") +
  labs(
    x = "Treatment",
    y = "Predicted bat calls",
    title = "Species-specific response of bat activity to artificial light"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_trmt

bm_plot <- bm2 %>%
  mutate(
    trmt = ifelse(trmt_bin == -1, "Dark", "Lit"),
    sp_label = str_replace(sp, "([A-Z])[A-Z]+([A-Z]+)", "\\1. \\2")
  )

# improved graph with raw data in the bakground and log y x to better visualize the differences between species.

p_trmt_m1.10 <- ggplot() +
  
  # 🔹 RAW DATA (background)
  geom_jitter(
    data = bm_plot,
    aes(x = trmt_bin, y = n),
    width = 0.15,
    alpha = 0.15,
    size = 0.6,
    color = "grey40"
  ) +
  
  # 🔹 MODEL PREDICTIONS (foreground)
  geom_point(
    data = pred_trmt,
    aes(x = trmt_bin, y = estimate),
    size = 2
  ) +
  
  geom_line(
    data = pred_trmt,
    aes(x = trmt_bin, y = estimate, group = sp),
    linewidth = 0.8
  ) +
  
  facet_wrap(~ sp, scales = "free_y") +
  
  labs(
    x = "Treatment",
    y = "Predicted bat calls",
    title = "Species-specific response of bat activity to artificial light"
  ) +
  scale_y_continuous(trans = "log1p")+
  scale_x_continuous(
    breaks = c(-1, 1),
    labels = c("Dark", "Lit")
  ) +
  
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_trmt_m1.10

ggsave(
  filename = "figures/glmm_v4/trmt_v1_acoustic_index.png",
  plot = p_trmt_m1.10,
  width = 10,
  height = 8,
  dpi = 300
)







###
# -------------------------------------------------------------------------


# below are previous models 

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
 # old modedls 
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

# Get mean and sd of original mwatts
mean_mwatts <- mean(bm2$mwatts, na.rm = TRUE)
sd_mwatts   <- sd(bm2$mwatts, na.rm = TRUE)

# Unstandardize the main prediction grid too
pred_w.y <- pred_w.y %>%
  mutate(mwatts = (mwatts_s * sd_mwatts) + mean_mwatts)

# Plot
ggplot(pred_w.y, aes(x = mwatts, y = estimate, color = factor(year))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = factor(year)),
              alpha = 0.2, color = NA) +
  labs(
    x = "watts",
    y = "Predicted bat calls (n)",
    color = "Year",
    fill  = "Year",
    title = "Interaction effect: watts × year"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# code suggested by chatgpt for depicting year effect. 

# Predictions at mean mwatts (0 after standardization)
pred_year <- predictions(
  m1.2nb,
  newdata = datagrid(
    yr_s = c(-1, 0, 1),
    mwatts_s = 0
  ),
  re.form = NA
) %>%
  mutate(
    year = case_when(
      yr_s == -1 ~ 2021,
      yr_s == 0  ~ 2022,
      yr_s == 1  ~ 2023
    )
  )

# Add to your interaction plot
ggplot(pred_w.y, aes(x = mwatts_s, y = estimate, color = factor(year))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = factor(year)),
              alpha = 0.2, color = NA) +
  geom_point(data = pred_year, aes(x = 0, y = estimate, color = factor(year)),
             size = 3, shape = 21, fill = "white") +  # highlights year differences at mean treatment
  geom_errorbar(data = pred_year,
                aes(x = 0, ymin = conf.low, ymax = conf.high, color = factor(year)),
                width = 0.1) +
  labs(
    x = "Standardized watts",
    y = "Predicted bat calls (n)",
    color = "Year",
    fill  = "Year",
    title = "Effect of year across treatments (interaction: watts × year)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


library(patchwork)

# Plot 1: Year effect at mean mwatts
p1 <- ggplot(pred_year, aes(x = factor(year), y = estimate, color = factor(year), fill = factor(year))) +
  geom_col(alpha = 0.4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(x = "Year", y = "Predicted bat calls (n)", title = "Effect of year (at mean watts)") +
  theme_minimal() +
  theme(legend.position = "none")
p1
# Plot 2: Interaction effect
p2 <- ggplot(pred_w.y, aes(x = mwatts_s, y = estimate, color = factor(year))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = factor(year)),
              alpha = 0.2, color = NA) +
  labs(x = "Standardized watts", y = "Predicted bat calls (n)", 
       color = "Year", fill = "Year", 
       title = "Interaction: watts × year") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine
p1 + p2

ggsave(
  filename = "figures/glmm_v3/year_v1.tiff", 
  plot = p1 + p2,   # combined plot goes here
  width = 20, 
  height = 10, 
  dpi = 300  # optional, good for publication-quality
)
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
p_jday<-ggplot(pred_jday, aes(x = jday, y = estimate)) +
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
p_sp_jday<-ggplot(pred_jday_sp, aes(x = jday, y = estimate, color = sp, fill = sp)) +
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

p_sp_jday<-ggplot(pred_jday_sp, aes(x = jday, y = estimate, color = sp, fill = sp)) +
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

p_sp_jday+p_jday

ggsave(
  filename = "figures/glmm_v3/jday_v1.tiff", 
  plot = p_sp_jday + p_jday,   # combined plot goes here
  width = 20, 
  height = 10, 
  dpi = 300  # optional, good for publication-quality
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
p_moon<-ggplot(pred_moon, aes(x = moon, y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(
    x = "Average moonlight",
    y = "Predicted bat calls (n)",
    title = "Community-level effect of moonlight"
  ) +
  theme_minimal()
p_moon

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


p_wind<-ggplot(pred_wind, aes(x = wind, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .5, fill = "skyblue", color = NA) +
  geom_line(color = "blue", linewidth = 1) +
  labs(
    x = "Wind speed (m/s)",
    y = "Predicted bat calls (n)",
    title = "Community-level effect of wind speed"
  ) +
  theme_minimal()


p_wind+p_moon

ggsave(
  filename = "figures/glmm_v3/wind_moon_v1.tiff", 
  plot = p_wind+p_moon,   # combined plot goes here
  width = 20, 
  height = 10, 
  dpi = 300  # optional, good for publication-quality
)


# save working environment ------------------------------------------------

save.image("working_env/glmm_v3.RData")
