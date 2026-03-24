
# =======================================================================
# Script Title:    rob_spkr_model.R
# Project:         robomoth and speaker data analysis
# =======================================================================

# Description:
# this script helps analyze the data I created with the robomoth_build.R and the rob_spkr_prepr.R scripts 

# Author: Carlos Linares
# Date:   2026-03-3
# Contact: carlosgarcialina@u.boisestate.edu

# =======================================================================
# Inputs:
# - write.csv(rob_db, file = 'data_for_analysis/rob_spkr_prep/rob_db.csv', row.names = F) 
# - write.csv(spkr_db, file = 'data_for_analysis/rob_spkr_prep/spkr_db', row.names = F)

#  Outputs:
# - Marginal effects plots in (images/rob_spkr_model/)
# =======================================================================

#  Load Libraries ------------------------------------------------------


# libraries 

if(!require("pacman")) install.packages("pacman")

pacman::p_load(
  "data.table",
  "lubridate",
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


# Define a function to calculate and print c-hat
calculate_c_hat <- function(model) {
  residual_deviance <- deviance(model)
  residual_df <- df.residual(model)
  c_hat_deviance <- residual_deviance / residual_df
  print(c_hat_deviance)
}


# load data 

rob_db <- fread('data_for_analysis/rob_spkr_prep/rob_db.csv') %>% 
  clean_names()
spkr_db <- fread('data_for_analysis/rob_spkr_prep/spkr_db') %>%
  clean_names()

glimpse(rob_db)
glimpse(spkr_db)

# filter out 2021 from rob_db because we don't have speaker data for that year

rob_db <- rob_db %>% 
  filter(year != 2021)


# models -------------------------------------------------------------------------

# possible model: occurrence of feeding using a binomial variable for those events where there is feeding and when there is no feedin. asks Does light changes the probability of feeding?
# feeding_present = c_buzz > 0
# feeding_present ~ trmt_bin + covariates
#(family = binomial)

# model intensity of feding mixing occurrence of feeding with intensity of feeding 
# this model askes: given fedding occurred, does light change how much feeding happens?

# leaving zeros or including zeros

m0 <- glmmTMB(
  c_buzz ~ trmt_bin + (1 | site),
  family = poisson(link = "log"),
  data = rob_db %>% filter(c_buzz > 0) # leaving zeros out
)

# check model
check_singularity(m0)
check_zeroinflation(m0)
calculate_c_hat(m0)  # indicates the is overdispersion we need to try negative binomial. 
performance_mae(m0)
range(rob_db$c_buzz)
hist(rob_db$c_buzz, breaks = 500)
performance::r2(m0)
DHARMa::simulateResiduals(m0, n = 1000, plot = TRUE)

# including zeros 


m0 <- glmmTMB(
  c_buzz ~ trmt_bin + (1 | site),
  family = poisson(link = "log"),
  data = rob_db  # leaving zeros in
)

# check model
check_singularity(m0)
check_zeroinflation(m0)
calculate_c_hat(m0)  # indicates the is overdispersion we need to try negative binomial. 
performance_mae(m0)
range(rob_db$c_buzz)
hist(rob_db$c_buzz, breaks = 500)
performance::r2(m0)
DHARMa::simulateResiduals(m0, n = 1000, plot = TRUE)

summary(m0)

# given the c_hat is 1.2 we might see overdisperssion. let's try negative binomial. 

m1 <- glmmTMB(
  c_buzz ~ trmt_bin + (1 | site),
  family = nbinom2,
  data = rob_db
)

# check model

check_singularity(m1)
check_zeroinflation(m1)
calculate_c_hat(m1)  # indicates the is overdispersion we need to try negative binomial.  
performance_mae(m1)
range(rob_db$c_buzz)
hist(rob_db$c_buzz, breaks = 100 )
performance::r2(m1)
DHARMa::simulateResiduals(m1, n = 1000, plot = TRUE)
summary(m1)

# phase 2
# now we add seasonality. 

m2 <- glmmTMB(
  c_buzz ~ trmt_bin + jday_s + I(jday_s^2) + (1 | site),
  family = nbinom2,
  data = rob_db
)

# check model

check_singularity(m2)
check_zeroinflation(m2)
calculate_c_hat(m2)  # indicates the is overdispersion we need to try negative binomial.  
performance_mae(m2)
range(rob_db$c_buzz)
hist(rob_db$c_buzz, breaks = 100 )
performance::r2(m2)
DHARMa::simulateResiduals(m2, n = 1000, plot = TRUE)
summary(m2)


#phase 3
# now we add environmental variables. 

m3 <- glmmTMB(
  c_buzz ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s + nit_avg_temp_c_s + avg_moonlight_s + elev_max + yr_s +
    (1 | site),
  family = nbinom2,
  data = rob_db
)

# check model

check_singularity(m3)
check_zeroinflation(m3)
calculate_c_hat(m3)  # indicates the is overdispersion we need to try negative binomial.  
performance_mae(m3)
range(rob_db$c_buzz)
hist(rob_db$c_buzz, breaks = 100 )
performance::r2(m3)
DHARMa::simulateResiduals(m3, n = 1000, plot = TRUE)
summary(m3)

# phase 4 we add the insect data 

m4 <- glmmTMB(
  c_buzz ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s + nit_avg_temp_c_s + avg_moonlight_s + elev_max + yr_s +
     + t_lepidoptera_s +
    (1 | site),
  family = nbinom2,
  data = rob_db
)

check_singularity(m4)
check_zeroinflation(m4)
calculate_c_hat(m4)  # indicates the is overdispersion we need to try negative binomial.  
performance_mae(m4)
range(rob_db$c_buzz)
hist(rob_db$c_buzz, breaks = 100 )
performance::r2(m4)
DHARMa::simulateResiduals(m4, n = 1000, plot = TRUE)
summary(m4)

anova(m0, m1, m2, m3, m4)


# Phase 5 interactions
# 

m5 <- glmmTMB(
  c_buzz ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s + nit_avg_temp_c_s + avg_moonlight_s + elev_max + yr_s +
    + t_lepidoptera_s +
    (1 | site)+
  #interactions
  trmt_bin*t_lepidoptera_s,
  family = nbinom2,
  data = rob_db
)
check_singularity(m5)
check_zeroinflation(m5)
calculate_c_hat(m5)  # indicates the is overdispersion we need to try negative binomial.  
performance_mae(m5)
range(rob_db$c_buzz)
hist(rob_db$c_buzz, breaks = 100 )
performance::r2(m5)
DHARMa::simulateResiduals(m5, n = 1000, plot = TRUE)
summary(m5)

anova(m0, m1, m2, m3, m4, m5)


# Phase 6 lets add random slopes by species

m6 <- glmmTMB(
  c_buzz ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s + nit_avg_temp_c_s + avg_moonlight_s + elev_max + yr_s +
    + t_lepidoptera_s +
    (1 | site)+ (1 | sp) +
    #interactions
    trmt_bin*t_lepidoptera_s,
  family = nbinom2,
  data = rob_db
)

check_singularity(m6)
check_zeroinflation(m6)
calculate_c_hat(m6)  # indicates the is overdispersion we need to try negative binomial.  
performance_mae(m6)
range(rob_db$c_buzz)
hist(rob_db$c_buzz, breaks = 100 )
performance::r2(m6)
DHARMa::simulateResiduals(m6, n = 1000, plot = TRUE)
summary(m6)

anova(m0, m1, m2, m3, m4, m5, m6)


# phase 7 lets add random slopes by species and site
m7 <- glmmTMB(
  c_buzz ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s + nit_avg_temp_c_s + avg_moonlight_s + elev_max + yr_s +
    + t_lepidoptera_s +
    (1 | site)+ (1 + trmt_bin | sp) +
    #interactions
    trmt_bin*t_lepidoptera_s,
  family = nbinom2,
  data = rob_db
)


check_singularity(m7)
check_zeroinflation(m7)
calculate_c_hat(m7)  # indicates the is overdispersion we need to try negative binomial.  
performance_mae(m7)
range(rob_db$c_buzz)
hist(rob_db$c_buzz, breaks = 100 )
performance::r2(m7)
DHARMa::simulateResiduals(m7, n = 1000, plot = TRUE)
summary(m7)

anova(m0, m1, m2, m3, m4, m5, m6, m7)
