
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
m7$sdr$pdHess
m7$fit$message
check_zeroinflation(m7)
calculate_c_hat(m7)   
performance_mae(m7)
range(rob_db$c_buzz)
hist(rob_db$c_buzz, breaks = 100 )
performance::r2(m7)
DHARMa::simulateResiduals(m7, n = 1000, plot = TRUE)
performance::check_collinearity(m7)
summary(m7)


anova(m0, m1, m2, m3, m4, m5, m6, m7)


# seasonality inside the random slope 
# the model breaks if we add seasonality inside the random slopes this.

# m8<- glmmTMB(
#   c_buzz ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s + nit_avg_temp_c_s + avg_moonlight_s + elev_max + yr_s +
#     + t_lepidoptera_s +
#     (1 | site)+ (1 + trmt_bin + jday_s + I(jday_s^2)| sp) +
#     #interactions
#     trmt_bin*t_lepidoptera_s,
#   family = nbinom2,
#   data = rob_db
# )
# 
# check_singularity(m8)
# check_zeroinflation(m8)
# calculate_c_hat(m8)   
# performance_mae(m8)
# range(rob_db$c_buzz)
# hist(rob_db$c_buzz, breaks = 100 )
# performance::r2(m8)
# DHARMa::simulateResiduals(m8, n = 1000, plot = TRUE)
# summary(m8)

# anova(m0, m1, m2, m3, m4, m5, m6, m7, )

# let's check what happens if we check for interaction between treatment and seasonality

m9 <- glmmTMB(
  c_buzz ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s + nit_avg_temp_c_s + avg_moonlight_s + elev_max + yr_s +
    + t_lepidoptera_s +
    (1 | site)+ (1 + trmt_bin | sp) +
    #interactions
    trmt_bin*t_lepidoptera_s + trmt_bin*jday_s + trmt_bin*I(jday_s^2),
  family = nbinom2,
  data = rob_db
)


check_singularity(m9)
m9$sdr$pdHess
m9$fit$message
check_zeroinflation(m9)
calculate_c_hat(m9)   
performance_mae(m9)
range(rob_db$c_buzz)
hist(rob_db$c_buzz, breaks = 100 )
performance::r2(m9)
DHARMa::simulateResiduals(m9, n = 1000, plot = TRUE)
performance::check_collinearity(m9)
summary(m9)


anova( m1, m2, m3, m4, m5, m6, m7, m9)

#now we add the moon interaction with treatment 

m10<- glmmTMB(
  c_buzz ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s + nit_avg_temp_c_s + avg_moonlight_s + elev_max + yr_s +
    + t_lepidoptera_s +
    (1 | site)+ (1 + trmt_bin | sp) +
    #interactions
    trmt_bin*t_lepidoptera_s + trmt_bin*jday_s + trmt_bin*I(jday_s^2) + trmt_bin*avg_moonlight_s,
  family = nbinom2,
  data = rob_db
) 

check_singularity(m10)
m10$sdr$pdHess
m10$fit$message
check_zeroinflation(m10)
calculate_c_hat(m10)   
performance_mae(m10)
range(rob_db$c_buzz)
hist(rob_db$c_buzz, breaks = 100 )
performance::r2(m10)
DHARMa::simulateResiduals(m10, n = 1000, plot = TRUE)
performance::check_collinearity(m10)
summary(m10)


anova( m1, m2, m3, m4, m5, m6, m7, m9, m10)

# now we remove lepidopter and treatment interaction becasuse it is not significant, it increases collinearity and is not really central to the question. 
# 

m11 <- glmmTMB(
  c_buzz ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s + nit_avg_temp_c_s + avg_moonlight_s + elev_max + yr_s   +t_lepidoptera_s +
    (1 | site) + (1 + trmt_bin | sp) +
    #interactions
    trmt_bin * jday_s + trmt_bin * I(jday_s^2) + trmt_bin * avg_moonlight_s,
  family = nbinom2,
  data = rob_db
)


check_singularity(m11)
m11$sdr$pdHess
m11$fit$message
check_zeroinflation(m11)
calculate_c_hat(m11)   
performance_mae(m11)
range(rob_db$c_buzz)
hist(rob_db$c_buzz, breaks = 100 )
performance::r2(m11)
DHARMa::simulateResiduals(m11, n = 1000, plot = TRUE)
performance::check_collinearity(m11)
summary(m11)
anova( m1, m2, m3, m4, m5, m6, m7, m9, m10, m11) # results indicate that the model m11 is the best model for now 


# Marginal effects  -------------------------------------------------------


# model mm11 --------------------------------------------------------------


# Lookup table: species code -> pretty label
sp_lookup <- rob_db %>%
  distinct(sp, sp_label)


typical <- list(
  jday_s = 0,
  nit_avg_wspm_s = 0,
  nit_avg_temp_c_s = 0,
  avg_moonlight_s = 0,
  t_lepidoptera_s = 0,
  elev_max = mean(rob_db$elev_max, na.rm = TRUE),
  yr_s = mean(rob_db$yr_s, na.rm = TRUE)
)

pred_comm <- predictions(
  m11,
  newdata = datagrid(
    trmt_bin = c(0, 1),
    jday_s = typical$jday_s,
    nit_avg_wspm_s = typical$nit_avg_wspm_s,
    nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
    avg_moonlight_s = typical$avg_moonlight_s,
    t_lepidoptera_s = typical$t_lepidoptera_s,
    elev_max = typical$elev_max,
    yr_s = typical$yr_s
  ),
  type = "response",
  re.form = NA
) %>%
  mutate(
    trmt = factor(trmt_bin, levels = c(0, 1), labels = c("Dark", "Lit")),
    panel = "Community"
  )

pred_sp <- predictions(
  m11,
  newdata = datagrid(
    sp = sort(unique(rob_db$sp)),
    trmt_bin = c(0, 1),
    jday_s = typical$jday_s,
    nit_avg_wspm_s = typical$nit_avg_wspm_s,
    nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
    avg_moonlight_s = typical$avg_moonlight_s,
    t_lepidoptera_s = typical$t_lepidoptera_s,
    elev_max = typical$elev_max,
    yr_s = typical$yr_s
  ),
  type = "response",
  re.form = NULL
) %>%
  left_join(sp_lookup, by = "sp") %>%
  mutate(
    trmt = factor(trmt_bin, levels = c(0, 1), labels = c("Dark", "Lit")),
    panel = if_else(is.na(sp_label), sp, sp_label)
  )

pred_all <- bind_rows(
  pred_comm %>% select(panel, trmt, estimate, conf.low, conf.high),
  pred_sp   %>% select(panel, trmt, estimate, conf.low, conf.high)
)

raw_sp <- rob_db %>%
  mutate(
    trmt = factor(trmt_bin, levels = c(0, 1), labels = c("Dark", "Lit")),
    panel = if_else(is.na(sp_label), sp, sp_label)
  )

raw_comm <- rob_db %>%
  mutate(
    trmt = factor(trmt_bin, levels = c(0, 1), labels = c("Dark", "Lit")),
    panel = "Community"
  )

raw_all <- bind_rows(
  raw_comm %>% select(panel, trmt, c_buzz),
  raw_sp   %>% select(panel, trmt, c_buzz)
)

panel_levels <- c("Community", sort(unique(raw_sp$panel)))

pred_all <- pred_all %>%
  mutate(panel = factor(panel, levels = panel_levels))

raw_all <- raw_all %>%
  mutate(panel = factor(panel, levels = panel_levels))


p_trmt_m11 <- ggplot() +
  geom_jitter(
    data = raw_all,
    aes(x = trmt, y = c_buzz),
    width = 0.12, height = 0,
    alpha = 0.18, size = 0.9
  ) +
  geom_point(
    data = pred_all,
    aes(x = trmt, y = estimate),
    size = 2
  ) +
  geom_errorbar(
    data = pred_all,
    aes(x = trmt, y = estimate, ymin = conf.low, ymax = conf.high),
    width = 0.08,
    linewidth = 0.5
  ) +
  geom_line(
    data = pred_all,
    aes(x = trmt, y = estimate, group = 1),
    linewidth = 0.6
  ) +
  facet_wrap(~ panel, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = "Treatment",
    y = "Feeding buzzes",
    title = "Marginal effect of treatment on feeding buzzes"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

p_trmt_m11


ggsave(
  p_trmt_m11,
  filename = "figures/rob_spkr_model/marginal_effects_trmt_m11.png",
  width = 8,
  height = 6,
  dpi = 300
)


# marginal effects for the jday

# sequence across observed range
j_seq <- seq(
  min(rob_db$jday_s, na.rm = TRUE),
  max(rob_db$jday_s, na.rm = TRUE),
  length.out = 100
)

# typical values for other predictors
typical <- list(
  nit_avg_wspm_s = 0,
  nit_avg_temp_c_s = 0,
  avg_moonlight_s = 0,
  t_lepidoptera_s = 0,
  elev_max = mean(rob_db$elev_max, na.rm = TRUE),
  yr_s = mean(rob_db$yr_s, na.rm = TRUE)
)

pred_jday <- predictions(
  m11,
  newdata = datagrid(
    jday_s = j_seq,
    trmt_bin = c(0, 1),
    nit_avg_wspm_s = typical$nit_avg_wspm_s,
    nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
    avg_moonlight_s = typical$avg_moonlight_s,
    t_lepidoptera_s = typical$t_lepidoptera_s,
    elev_max = typical$elev_max,
    yr_s = typical$yr_s
  ),
  type = "response",
  re.form = NA   # <-- community-level
) %>%
  mutate(
    trmt = factor(trmt_bin, levels = c(0, 1), labels = c("Dark", "Lit"))
  )

p_jday <- ggplot(pred_jday, aes(x = jday_s, y = estimate, color = trmt)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = trmt),
    alpha = 0.2,
    color = NA
  ) +
  labs(
    x = "Season (scaled Julian day)",
    y = "Predicted feeding buzzes",
    color = "Treatment",
    fill = "Treatment",
    title = "Seasonal effect on feeding buzzes by light treatment"
  ) +
  theme_minimal(base_size = 13)

p_jday


ggsave(
  "figures/rob_spkr_model/p_jday.tiff", # this graph needs to be fixed show average increas (misleading),
  width = 22,
  height = 15, 
  dpi = 600,
  compression = "lzw",
  units = "cm"
)


# effect of night average temp

# sequence across observed range
wspm_seq <- seq(
  min(rob_db$nit_avg_wspm_s, na.rm = TRUE),
  max(rob_db$nit_avg_wspm_s, na.rm = TRUE),
  length.out = 100
)

typical <- list(
  avg_moonlight_s = 0,
  t_lepidoptera_s = 0,
  elev_max = mean(rob_db$elev_max, na.rm = TRUE),
  yr_s = mean(rob_db$yr_s, na.rm = TRUE)
)

pred_wspm <- predictions(
  m11,
  newdata = datagrid(
    nit_avg_wspm_s = wspm_seq,
    trmt_bin = mean(c(0, 1)),  # average treatment effect
    jday_s = 0,
    nit_avg_temp_c_s = 0,
    avg_moonlight_s = typical$avg_moonlight_s,
    t_lepidoptera_s = typical$t_lepidoptera_s,
    elev_max = typical$elev_max,
    yr_s = typical$yr_s
  ),
  type = "response",
  re.form = NA
) %>%
  mutate(
    trmt = factor(trmt_bin, levels = c(0, 1), labels = c("Dark", "Lit"))
  )


ggplot(pred_wspm, aes(x = nit_avg_wspm_s, y = estimate)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  labs(
    x = "Wind speed (scaled)",
    y = "Predicted feeding buzzes",
    color = "Treatment",
    fill = "Treatment",
    title = "Marginal effect of wind speed on feeding buzzes"
  ) +
  theme_minimal(base_size = 13)


# moon marginal effect 

moon_seq <- seq(
  min(rob_db$avg_moonlight_s, na.rm = TRUE),
  max(rob_db$avg_moonlight_s, na.rm = TRUE),
  length.out = 100
)

typical <- list(
  jday_s = 0,
  nit_avg_wspm_s = 0,
  nit_avg_temp_c_s = 0,
  t_lepidoptera_s = 0,
  elev_max = mean(rob_db$elev_max, na.rm = TRUE),
  yr_s = mean(rob_db$yr_s, na.rm = TRUE)
)

pred_moon <- predictions(
  m11,
  newdata = datagrid(
    avg_moonlight_s = moon_seq,
    trmt_bin = c(0, 1),
    jday_s = typical$jday_s,
    nit_avg_wspm_s = typical$nit_avg_wspm_s,
    nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
    t_lepidoptera_s = typical$t_lepidoptera_s,
    elev_max = typical$elev_max,
    yr_s = typical$yr_s
  ),
  type = "response",
  re.form = NA
) %>%
  mutate(
    trmt = factor(trmt_bin, levels = c(0, 1), labels = c("Dark", "Lit"))
  )

p_moon <- ggplot(pred_moon, aes(x = avg_moonlight_s, y = estimate, color = trmt, fill = trmt)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  labs(
    x = "Moonlight (scaled)",
    y = "Predicted feeding buzzes",
    color = "Treatment",
    fill = "Treatment",
    title = "Marginal effect of moonlight on feeding buzzes"
  ) +
  theme_minimal(base_size = 13)

p_moon


ggsave(
  "figures/rob_spkr_model/p_moon.tiff", # this graph needs to be fixed show average increas (misleading),
  width = 22,
  height = 15, 
  dpi = 600,
  compression = "lzw",
  units = "cm"
)


# Marginal effects plot for year

yr_seq <- seq(
  min(rob_db$yr_s, na.rm = TRUE),
  max(rob_db$yr_s, na.rm = TRUE),
  length.out = 100
)

typical <- list(
  jday_s = 0,
  nit_avg_wspm_s = 0,
  nit_avg_temp_c_s = 0,
  avg_moonlight_s = 0,
  t_lepidoptera_s = 0,
  elev_max = mean(rob_db$elev_max, na.rm = TRUE)
)

pred_year_disc <- predictions(
  m11,
  newdata = datagrid(
    yr_s = c(-1, 1),
    trmt_bin = 0.5,
    jday_s = 0,
    nit_avg_wspm_s = 0,
    nit_avg_temp_c_s = 0,
    avg_moonlight_s = 0,
    t_lepidoptera_s = 0,
    elev_max = mean(rob_db$elev_max, na.rm = TRUE)
  ),
  type = "response",
  re.form = NA
) %>%
  mutate(
    year = factor(yr_s, levels = c(-1, 0, 1),
                  labels = c("2021", "2022", "2023"))
  )

ggplot(pred_year_disc, aes(x = year, y = estimate)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  geom_line(aes(group = 1)) +
  labs(
    x = "Year",
    y = "Predicted feeding buzzes",
    title = "Marginal effect of year on feeding buzzes"
  ) +
  theme_minimal(base_size = 13)

raw_year <- rob_db %>%
  mutate(
    year = factor(yr_s, levels = c(-1, 0, 1),
                  labels = c("2021", "2022", "2023"))
  )

ggplot() +
  geom_jitter(
    data = raw_year,
    aes(x = year, y = c_buzz),
    alpha = 0.05,
    width = 0.1
  ) +
  geom_point(
    data = pred_year_disc,
    aes(x = year, y = estimate),
    size = 2
  ) +
  geom_errorbar(
    data = pred_year_disc,
    aes(x = year, ymin = conf.low, ymax = conf.high),
    width = 0.1
  ) +
  labs(
    x = "Year",
    y = "Feeding buzzes",
    title = "Year effect on feeding buzzes"
  ) +
  scale_y_continuous(trans = "log1p") +
  theme_minimal(base_size = 13)
