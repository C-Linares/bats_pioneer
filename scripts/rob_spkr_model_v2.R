
# =======================================================================
# Script Title:    rob_spkr_model.R
# Project:         robomoth and speaker data analysis
# =======================================================================

# Description:
# this script helps analyze the data I created with the robomoth_build_v2.R and the updated version of  rob_spkr_prepr_v.R scripts
# now we filter calls using an amplitude threshold and analyze passes and not feeding buzzes. And we have filtered calls >= -20 dBfs

# Author: Carlos Linares
# Date:   2026-03-3
# Contact: carlosgarcialina@u.boisestate.edu
# version : v3 updated for the models using the v3 version of databases. 

# =======================================================================
# Inputs:
# - robomoth data 'data_for_analysis/rob_spkr_prep/rob_db_v3.csv' 
# - speakr data 'data_for_analysis/rob_spkr_prep/spkr_db_v3.csv'

#  Outputs:
# - Marginal effects plots for the robomoth and speaker data showing the effect of light on bat calls/passes under an amplitude threshold. 

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
  "janitor",
  "flextable"
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

rob_db <- fread('data_for_analysis/rob_spkr_prep_v3/rob_db_v3.csv') %>% 
  clean_names()



spkr_db <- fread('data_for_analysis/rob_spkr_prep_v3/spkr_db_v3.csv') %>% 
  clean_names()

glimpse(rob_db)
glimpse(spkr_db)




# models -------------------------------------------------------------------------

# possible model: occurrence of feeding using a binomial variable for those events where there is feeding and when there is no feeding. asks Does light changes the probability of feeding?
# feeding_present = c_buzz > 0
# feeding_present ~ trmt_bin + covariates
#(family = binomial)

# model intensity of feeding mixing occurrence of feeding with intensity of feeding 
# this model askes: given fedding occurred, does light change how much feeding happens?

# leaving zeros or including zeros

m0 <- glmmTMB(
  n_calls ~ trmt_bin + (1 | site),
  family = poisson(link = "log"),
  data = rob_db
)

# check model
check_singularity(m0)
check_zeroinflation(m0) # probably zero inflated. 
calculate_c_hat(m0)  # indicates the is overdispersion we need to try negative binomial. 
performance_mae(m0)
range(rob_db$n_calls)
hist(rob_db$n_calls, breaks = 50)
performance::r2(m0)
DHARMa::simulateResiduals(m0, n = 1000, plot = TRUE)

# including zeros 



# given the c_hat is big 40.025 it indicates overdissperssion. let's try negative binomial. 

m1 <- glmmTMB(
  n_calls ~ trmt_bin + (1 | site),
  family = nbinom2,
  data = rob_db
)

# check model

check_singularity(m1)
m1$sdr$pdHess
m1$fit$message
check_zeroinflation(m1) # not zero inflated.
calculate_c_hat(m1)    # no overdispersion, c_hat is 1.08 which is good.
performance_mae(m1)
performance::r2(m1)
DHARMa::simulateResiduals(m1, n = 1000, plot = TRUE)
summary(m1)

# phase 2
# now we add seasonality. 

m2 <- glmmTMB(
  n_calls ~ trmt_bin + jday_s + I(jday_s^2) + (1 | site),
  family = nbinom2,
  data = rob_db
)

# check model

check_singularity(m2)
m2$sdr$pdHess
m2$fit$message
check_zeroinflation(m2)
calculate_c_hat(m2)    
performance_mae(m2)
performance::r2(m2)
DHARMa::simulateResiduals(m2, n = 1000, plot = TRUE)
summary(m2)

anova(m0, m1, m2) # the model with seasonality is better than the model without seasonality.

#phase 3
# now we add environmental variables. 

m3 <- glmmTMB(
  n_calls ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s_s + avg_moonlight_s + elev_max_s + yr_s +
    (1 | site),
  family = nbinom2,
  data = rob_db
)

# check model

check_singularity(m3)
m3$sdr$pdHess
m3$fit$message
check_zeroinflation(m3)
calculate_c_hat(m3)    
performance_mae(m3)
performance::r2(m3)
DHARMa::simulateResiduals(m3, n = 1000, plot = TRUE)
summary(m3)

# phase 4 we add the insect data 

m4 <- glmmTMB(
  n_calls ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s_s  + avg_moonlight_s + elev_max_s + yr_s +
     + t_leps_s +
    (1 | site),
  family = nbinom2,
  data = rob_db
)

check_singularity(m4)
m4$sdr$pdHess
m4$fit$message
check_zeroinflation(m4)
calculate_c_hat(m4)    
performance_mae(m4)
performance::r2(m4)
DHARMa::simulateResiduals(m4, n = 1000, plot = TRUE)
summary(m4)

anova(m0, m1, m2, m3, m4)

# m4 does not seem to improve significantly when adding insects. compared to m3 


# Phase 5 interactions aleps and treatment 
# 

m5 <- glmmTMB(
  n_calls ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s_s  + avg_moonlight_s + elev_max_s + yr_s 
    + t_leps_s +
    (1 | site)+
  #interactions
  trmt_bin*t_leps_s,
  family = nbinom2,
  data = rob_db
)
check_singularity(m5)
m5$sdr$pdHess
m5$fit$message
check_zeroinflation(m5)
calculate_c_hat(m5)  #  
performance_mae(m5)
performance::r2(m5)
DHARMa::simulateResiduals(m5, n = 1000, plot = TRUE)
summary(m5)

anova(m0, m1, m2, m3, m4, m5)

# m5 does improve significantly compared to m3,4

# Phase 6 lets add random slopes by species

m6 <- glmmTMB(
  n_calls ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s_s  + avg_moonlight_s + elev_max_s + yr_s + t_leps_s +
    (1 | site)+ (1 | sp) +
    #interactions
    trmt_bin*t_leps_s,
  family = nbinom2,
  data = rob_db
)

check_singularity(m6)
m6$sdr$pdHess
m6$fit$message
check_zeroinflation(m6)
calculate_c_hat(m6)  #  
performance_mae(m6)
performance::r2(m6)
DHARMa::simulateResiduals(m6, n = 1000, plot = TRUE)
summary(m6)

anova(m0, m1, m2, m3, m4, m5, m6)

#m6 does improve the model fit. 

# phase 7 lets add random slopes by species and treatment. but this time we add the offset for background calls. 

m7_bg <- glmmTMB(
  n_calls ~ trmt_bin + jday_s + I(jday_s^2) +
    nit_avg_wspm_s_s + avg_moonlight_s + elev_max_s +
    yr_s + t_leps_s +
    offset(log_background_offset) +
    (1 | site) +
    (1 + trmt_bin | sp) +
    trmt_bin * t_leps_s,
  family = nbinom2,
  data = rob_db
)


check_singularity(m7_bg)
m7_bg$sdr$pdHess
m7_bg$fit$message
check_zeroinflation(m7_bg)
calculate_c_hat(m7_bg)   
performance_mae(m7_bg)
performance::r2(m7_bg)
DHARMa::simulateResiduals(m7_bg, n = 1000, plot = TRUE)
performance::check_collinearity(m7_bg)
summary(m7_bg)




# seasonality inside the random slope
# the model breaks if we add seasonality inside the random slopes this.

# m8<- glmmTMB(
#   n_calls ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s_s + nit_avg_temp_c_s + avg_moonlight_s + elev_max_s + yr_s +
#     + t_leps_s +
#     (1 | site)+ (1 + trmt_bin + jday_s + I(jday_s^2)| sp) +
#     #interactions
#     trmt_bin*t_leps_s,
#   family = nbinom2,
#   data = rob_db
# )
# 
# # #
# check_singularity(m8)
# m8$sdr$pdHess
# m8$fit$message
# check_zeroinflation(m8)
# calculate_c_hat(m8)
# performance_mae(m8)
# range(rob_db$c_buzz)
# hist(rob_db$c_buzz, breaks = 100 )
# performance::r2(m8)
# DHARMa::simulateResiduals(m8, n = 1000, plot = TRUE)
# summary(m8)
# 
# anova(m0, m1, m2, m3, m4, m5, m6, m7, m8)

# let's check what happens if we check for interaction between treatment and seasonality. It seems model 9 is not better than m7

m9_bg <- glmmTMB(
  n_calls ~ trmt_bin + jday_s + I(jday_s^2) +
    nit_avg_wspm_s_s + avg_moonlight_s + elev_max_s +
    yr_s + t_leps_s +
    offset(log_background_offset) +
    (1 | site) +
    (1 + trmt_bin | sp) +
    trmt_bin * t_leps_s +
    trmt_bin * jday_s +
    trmt_bin * I(jday_s^2),
  family = nbinom2,
  data = rob_db
)


check_singularity(m9_bg)
m9_bg$sdr$pdHess
m9_bg$fit$message
check_zeroinflation(m9_bg)
calculate_c_hat(m9_bg)   
performance_mae(m9_bg)
performance::r2(m9_bg)
DHARMa::simulateResiduals(m9_bg, n = 1000, plot = TRUE)
performance::check_collinearity(m9_bg)
summary(m9_bg)



anova( m7_bg, m9_bg)

# m9 is not better and than m7 


#now we add the moon interaction with treatment. but we remove the seasonality interaction as it does not improve the model. 

m10_bg <- glmmTMB(
  n_calls ~ trmt_bin + jday_s + I(jday_s^2) +
    nit_avg_wspm_s_s + avg_moonlight_s + elev_max_s +
    yr_s + t_leps_s +
    offset(log(effort_hours)) +
    offset(log_background_offset) +
    (1 | site) +
    (1 + trmt_bin | sp) +
    trmt_bin * t_leps_s +
    trmt_bin * avg_moonlight_s,
  family = nbinom2,
  data = rob_db
)

check_singularity(m10_bg)
m10_bg$sdr$pdHess
m10_bg$fit$message
check_zeroinflation(m10_bg)
calculate_c_hat(m10_bg)
performance_mae(m10_bg)
performance::r2(m10_bg)
DHARMa::simulateResiduals(m10_bg, n = 1000, plot = TRUE)
performance::check_collinearity(m10_bg)
summary(m10_bg)



anova( m7_bg, m9_bg, m10_bg) #the model doe improve when we add the treatment and moonlight interaction 


# now we check the model with the conservative data where we just keep data that has n_background > 0 indicating those where we have some background calls from sm3 and for those effort hours > 0 as well. 

rob_db_offset_conservative <- rob_db %>%
  filter(
    effort_hours > 0,
    n_background > 0,
    !is.na(n_calls),
    !is.na(trmt_bin),
    !is.na(jday_s),
    !is.na(nit_avg_wspm_s_s),
    !is.na(avg_moonlight_s),
    !is.na(elev_max_s),
    !is.na(yr_s),
    !is.na(t_leps_s),
    !is.na(site),
    !is.na(sp)
  )


# we run the same model as m10_bg but with the conservative data 

m10_bg_conservative <- update(
  m10_bg,
  data = rob_db_offset_conservative
)

check_singularity(m10_bg_conservative)
m10_bg_conservative$sdr$pdHess
m10_bg_conservative$fit$message
check_zeroinflation(m10_bg_conservative)
calculate_c_hat(m10_bg_conservative)
performance_mae(m10_bg_conservative)
performance::r2(m10_bg_conservative)
DHARMa::simulateResiduals(m10_bg_conservative, n = 1000, plot = TRUE)
performance::check_collinearity(m10_bg_conservative)


summary(m10_bg_conservative)




# it seems like the best model is m10_bg even when compared with the conservative model m10_bg_conservative. i might need to do a table for this one instead. 

# in the past models leps interaction with treatment was not significant but now it is. I won't use m11 because it removes the interaction. 


# m11 <- glmmTMB(
#   n_calls ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s_s   + avg_moonlight_s + elev_max_s + yr_s   + t_leps_s +
#     (1 | site) + (1 + trmt_bin | sp) +
#     #interactions
#     trmt_bin * jday_s + trmt_bin * I(jday_s^2) + trmt_bin * avg_moonlight_s,
#   family = nbinom2,
#   data = rob_db
# )

# 
# check_singularity(m11)
# m11$sdr$pdHess
# m11$fit$message
# check_zeroinflation(m11)
# calculate_c_hat(m11)   
# performance_mae(m11)
# range(rob_db$c_buzz)
# performance::r2(m11)
# DHARMa::simulateResiduals(m11, n = 1000, plot = TRUE)
# performance::check_collinearity(m11)
# summary(m11)
# anova( m1, m2, m3, m4, m5, m6, m7, m9, m10, m11) # after testing AIC shos m10 is could be better candidate but If I want to keep it simple and focus only in the answer. then the model m11 works as well


m10_table <- tbl_regression(
  m10_bg,
  exponentiate = FALSE
)

m10_table <- gtsummary::tbl_regression(
  m10_bg,
  exponentiate = TRUE,
  label = list(
    trmt_bin ~ "Treatment (Lit vs Dark)",
    jday_s ~ "Julian Day ",
    `I(jday_s^2)` ~ "Julian Day Squared",
    avg_moonlight_s ~ "Average Moonlight ",
    elev_max_s  ~ "Maximum Elevation",
    nit_avg_wspm_s_s ~ "Nigt Average Wind Speed",
    t_leps_s ~ "Total Lepidoptera",
    yr_s ~ "Year",
    `trmt_bin:avg_moonlight_s` ~ "Treatment × Moonlight"
  )
)

m10_table

m10_flex <- gtsummary::as_flex_table(m10_table) %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::fontsize(size = 14, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::autofit() %>% 
  # bold sig values
  flextable::bold(
    i = ~ p.value < 0.05,
    j = "p.value",
    bold = TRUE,
    part = "body"
  )

m10_flex

save_as_docx(
  "m10_bg_table" = m10_flex,
  path = "figures/rob_spkr_model/v3/tables/m10_bg_table.docx")


save_as_image(m10_flex, path = "figures/rob_spkr_model/v3/tables/m10_bg_table.png")



####################################
####################################
####################################
# Marginal effects robomoth -------------------------------------------------------

# we are using the m11 model to produce marginal effects plots for the robomoth data.

# Treatment colors used in your previous figures
trt_cols <- c(
  "Dark" = "grey10",
  "Lit"  = "grey"
)

# Species order
rob_species <- rob_db %>%
  filter(!is.na(sp)) %>%
  distinct(sp) %>%
  arrange(sp) %>%
  pull(sp)

# typical offsevalue

typical_log_background<- median(
  rob_db$log_background_offset,
  na.rm = TRUE
)

typical_background<- exp(typical_log_background)

typical_background


# typical effort hours 
typical_effort_hours <- median(
  rob_db$effort_hours,
  na.rm = TRUE
)

typical_effort_hours

# -------------------------------------------------------------------------
# Species-specific prediction grid
# -------------------------------------------------------------------------

pred_grid_species <- tidyr::expand_grid(
  sp = rob_species,
  trmt_bin = c(-1, 1)
) %>%
  mutate(
    site = NA,                 # removes site-specific random intercept
    jday_s = 0,
    nit_avg_wspm_s_s = 0,
    avg_moonlight_s = 0,
    elev_max_s = 0,
    yr_s = 0,
    t_leps_s = 0,
    
    # required ofsset 
    
    log_background_offset = typical_log_background, # is around 6
    # Required because the new model includes offset(log(effort_hours))
    effort_hours = typical_effort_hours
  )

is.factor(rob_db$site)
is.factor(rob_db$sp)

# -------------------------------------------------------------------------
# Species-specific predictions
# -------------------------------------------------------------------------

pred_species <- predictions(
  m10_bg,
  newdata = pred_grid_species,
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    ),
    panel = as.character(sp)
  )

# -------------------------------------------------------------------------
# Community-average predictions across species
# -------------------------------------------------------------------------

pred_community <- avg_predictions(
  m10_bg,
  newdata = pred_grid_species,
  by = "trmt_bin",
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    ),
    panel = "Community"
  )

# -------------------------------------------------------------------------
# Combine species-specific and community predictions
# -------------------------------------------------------------------------

pred_treatment_all <- bind_rows(
  pred_community,
  pred_species
) %>%
  mutate(
    panel = factor(panel, levels = c("Community", rob_species))
  )

# -------------------------------------------------------------------------
# Marginal predictions plot: treatment effect by species
# -------------------------------------------------------------------------

p_1 <- ggplot(
  pred_treatment_all,
  aes(x = treatment, y = estimate)
) +
  geom_line(
    aes(group = panel),
    color = "black",
    linewidth = 0.7
  ) +
  geom_pointrange(
    aes(
      ymin = conf.low,
      ymax = conf.high,
      color = treatment
    ),
    linewidth = 0.7,
    size = 1.5
  ) +
  facet_wrap(
    ~ panel,
    scales = "free_y"
  ) +
  scale_color_manual(values = trt_cols) +
  labs(
    title = "Predicted faint robomoth passes by treatment",
    subtitle = "Predictions from m11; covariates held at mean scaled values",
    x = "Treatment",
    y = "Predicted number of faint bat passes",
    color = "Treatment"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

p_1

library(tidyverse)
library(marginaleffects)

# Species in robomoth data
rob_species <- rob_db %>%
  filter(!is.na(sp)) %>%
  distinct(sp) %>%
  arrange(sp) %>%
  pull(sp)

# Prediction grid at mean covariate values
rob_trt_grid <- tidyr::expand_grid(
  sp = rob_species,
  trmt_bin = -1
) %>%
  mutate(
    site = NA,
    # Mean scaled covariate values
    jday_s = 0,
    nit_avg_wspm_s_s = 0,
    avg_moonlight_s = 0,
    elev_max_s = 0,
    yr_s = 0,
    t_leps_s = 0,
    
    # Required offset value
    log_background_offset = typical_log_background,
    # Required because the new model includes offset(log(effort_hours))
    effort_hours = typical_effort_hours
  )



# Species-specific lit/dark ratios
rob_species_ratio <- avg_comparisons(
  m10_bg,
  newdata = rob_trt_grid,
  variables = list(trmt_bin = c(-1, 1)),
  comparison = "lnratio",
  by = "sp",
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    percent_change_lit_vs_dark = 100 * (estimate - 1),
    percent_conf_low = 100 * (conf.low - 1),
    percent_conf_high = 100 * (conf.high - 1),
    supported_response = case_when(
      conf.low > 1 ~ "Higher in lit",
      conf.high < 1 ~ "Lower in lit",
      TRUE ~ "Uncertain"
    )
  ) %>%
  arrange(desc(percent_change_lit_vs_dark))

rob_species_ratio

rob_species_ratio_clean <- rob_species_ratio %>%
  transmute(
    species = sp,
    lit_dark_ratio = round(estimate, 2),
    ratio_CI = paste0(round(conf.low, 2), "–", round(conf.high, 2)),
    percent_change = round(percent_change_lit_vs_dark, 1),
    percent_CI = paste0(
      round(percent_conf_low, 1),
      " to ",
      round(percent_conf_high, 1),
      "%"
    ),
    p_value = signif(p.value, 3),
    supported_response
  )

rob_species_ratio_clean

ts4 <- flextable(rob_species_ratio_clean) %>%
  theme_booktabs() %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  bold(
    i = ~ supported_response != "Uncertain",
    bold = TRUE,
    part = "body"
  ) %>%
  autofit() %>%
  set_caption(
    "Table S4. Species-specific lit/dark ratios for faint robomoth passes from the abundance-offset robomoth model."
  )

ts4

save_as_docx(
  "table.s4" = ts4,
  path = "figures/rob_spkr_model/v3/tables/ts4_robomoth.docx")

save_as_image(ts4, path = "figures/rob_spkr_model/v3/tables/ts4_robomoth.png")

# change to log ratio for the comparaison 

rob_species_logratio <- avg_comparisons(
  m10_bg,
  newdata = rob_trt_grid,
  variables = list(trmt_bin = c(-1, 1)),
  comparison = "lnratio",
  by = "sp",
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    # Back-transform from log ratio to lit/dark ratio
    lit_dark_ratio = exp(estimate),
    ratio_conf_low = exp(conf.low),
    ratio_conf_high = exp(conf.high),
    
    # Convert ratio to percent change
    percent_change_lit_vs_dark = 100 * (lit_dark_ratio - 1),
    percent_conf_low = 100 * (ratio_conf_low - 1),
    percent_conf_high = 100 * (ratio_conf_high - 1),
    
    # On the log-ratio scale, 0 means no difference.
    supported_response = case_when(
      conf.low > 0 ~ "Higher in lit",
      conf.high < 0 ~ "Lower in lit",
      TRUE ~ "Uncertain"
    )
  ) %>%
  arrange(desc(percent_change_lit_vs_dark))

rob_species_logratio


# then we have a clean table 

rob_species_ratio_clean <- rob_species_logratio %>%
  transmute(
    species = sp,
    lit_dark_ratio = round(lit_dark_ratio, 2),
    ratio_CI = paste0(
      round(ratio_conf_low, 2),
      "–",
      round(ratio_conf_high, 2)
    ),
    percent_change = round(percent_change_lit_vs_dark, 1),
    percent_CI = paste0(
      round(percent_conf_low, 1),
      " to ",
      round(percent_conf_high, 1),
      "%"
    ),
    p_value = signif(p.value, 3),
    supported_response
  )

rob_species_ratio_clean

ts4 <- flextable(rob_species_ratio_clean) %>%
  theme_booktabs() %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  bold(
    i = ~ supported_response != "Uncertain",
    bold = TRUE,
    part = "body"
  ) %>%
  autofit() %>%
  set_caption(
    "Table S4. Species-specific lit/dark ratios for faint robomoth passes from the abundance-offset robomoth model."
  )

ts4

save_as_docx(
  "table.s4" = ts4,
  path = "figures/rob_spkr_model/v3/tables/ts4_robomoth.docx")

save_as_image(ts4, path = "figures/rob_spkr_model/v3/tables/ts4_robomoth.png")


# jday robomot marginal  --------------------------------------------------

# Treatment colors
trt_cols <- c(
  "Dark" = "grey20",
  "Lit"  = "grey"
)

# Get scaling values from your data
# This assumes jday_s was scaled as: (jday - mean(jday)) / (2 * sd(jday))
jday_mean <- mean(rob_db$jday, na.rm = TRUE)
jday_sd   <- sd(rob_db$jday, na.rm = TRUE)

# Raw Julian day sequence for plotting
jday_seq <- seq(
  min(rob_db$jday, na.rm = TRUE),
  max(rob_db$jday, na.rm = TRUE),
  length.out = 100
)

# Species included in model/data
rob_species <- rob_db %>%
  filter(!is.na(sp)) %>%
  distinct(sp) %>%
  arrange(sp) %>%
  pull(sp)

# Prediction grid
pred_grid_jday <- tidyr::expand_grid(
  jday = jday_seq,
  trmt_bin = c(-1, 1),
  sp = rob_species
) %>%
  mutate(
    jday_s = (jday - jday_mean) / (2 * jday_sd),
    site = NA,
    nit_avg_wspm_s_s = 0,
    avg_moonlight_s = 0,
    elev_max_s = 0,
    yr_s = 0,
    t_leps_s = 0,
    # Required offset value
    log_background_offset = typical_log_background,
    # Required because the new model includes offset(log(effort_hours))
    effort_hours = typical_effort_hours
  )


pred_jday_comm <- avg_predictions(
  m10_bg,
  newdata = pred_grid_jday,
  by = c("jday", "jday_s", "trmt_bin"),
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    )
  )


p_jday_comm <- ggplot(
  pred_jday_comm,
  aes(x = jday, y = estimate, color = treatment, fill = treatment)
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = trt_cols) +
  scale_fill_manual(values = trt_cols) +
  labs(
    title = "Seasonal pattern of robomoth bat passes (>= -40 dBfs) by treatment",
    subtitle = "Predictions from m10_bg; non-focal covariates held at mean scaled values",
    x = "Julian day",
    y = "Predicted number of bat passes",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

p_jday_comm


raw_jday_summary <- rob_db %>%
  group_by(treatmt, trmt_bin, jday) %>%
  summarise(
    mean_calls = mean(n_calls, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    )
  )



p_jday_comm_raw <- ggplot(
  pred_jday_comm,
  aes(x = jday, y = estimate, color = treatment, fill = treatment)
) +
  geom_point(
    data = raw_jday_summary,
    aes(x = jday, y = mean_calls, color = treatment),
    alpha = 0.25,
    size = 1.4,
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = trt_cols) +
  scale_fill_manual(values = trt_cols) +
  labs(
    # title = "Seasonal pattern of faint robomoth bat passes by treatment",
    # subtitle = "Points are observed treatment-level means; lines are model predictions from m11",
    # x = "Julian day",
    y = "Predicted number of faint bat passes",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

p_jday_comm_raw



# year robomoth --------------------------------------------------------------------



# Treatment colors
trt_cols <- c(
  "Dark" = "grey35",
  "Lit"  = "goldenrod2"
)

# Extract the real year-to-yr_s mapping from your model data
year_lookup <- rob_db %>%
  filter(!is.na(year), !is.na(yr_s)) %>%
  distinct(year, yr_s) %>%
  arrange(year)

year_lookup


# Species included in the model
rob_species <- rob_db %>%
  filter(!is.na(sp)) %>%
  distinct(sp) %>%
  arrange(sp) %>%
  pull(sp)

# Prediction grid
pred_grid_year <- tidyr::expand_grid(
  year = year_lookup$year,
  trmt_bin = c(-1, 1),
  sp = rob_species
) %>%
  left_join(year_lookup, by = "year") %>%
  mutate(
    site = NA,
    jday_s = 0,
    nit_avg_wspm_s_s = 0,
    avg_moonlight_s = 0,
    elev_max_s = 0,
    t_leps_s = 0,
    # Required offset value
    log_background_offset = typical_log_background,
    # Required because the new model includes offset(log(effort_hours))
    effort_hours = typical_effort_hours
  )


# -------------------------------------------------------------------------
# Community-level year predictions by treatment
# -------------------------------------------------------------------------

pred_year_trt <- avg_predictions(
  m10_bg,
  newdata = pred_grid_year,
  by = c("year", "trmt_bin"),
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    )
  )

pred_year_trt


p_year_trt <- ggplot(
  pred_year_trt,
  aes(x = year, y = estimate, color = treatment, fill = treatment)
) +
  geom_line(linewidth = 1.1) +
  geom_pointrange(
    aes(ymin = conf.low, ymax = conf.high),
    linewidth = 0.7,
    fatten = 2.3
  ) +
  scale_color_manual(values = trt_cols) +
  scale_fill_manual(values = trt_cols) +
  scale_x_continuous(
    breaks = year_lookup$year
  ) +
  labs(
    title = "Predicted faint robomoth bat passes by year and treatment",
    subtitle = "Predictions from m11; non-focal covariates held at mean scaled values",
    x = "Year",
    y = "Predicted number of faint bat passes",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

p_year_trt



# -------------------------------------------------------------------------
# Community-level year predictions averaged across treatment and species
# -------------------------------------------------------------------------

pred_year_comm <- avg_predictions(
  m10_bg,
  newdata = pred_grid_year,
  by = "year",
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble()

p_year_comm <- ggplot(
  pred_year_comm,
  aes(x = year, y = estimate)
) +
  geom_line(linewidth = 1.1, color = "black") +
  geom_pointrange(
    aes(ymin = conf.low, ymax = conf.high),
    linewidth = 0.7,
    fatten = 2.3,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = year_lookup$year
  ) +
  labs(
    title = "Predicted faint robomoth bat passes by year",
    subtitle = "Community-level predictions from m11 averaged across treatments and species",
    x = "Year",
    y = "Predicted number of faint bat passes"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

p_year_comm




# moon robmoth marginal  --------------------------------------------------


# Treatment colors
trt_cols <- c(
  "Dark" = "grey20",
  "Lit"  = "grey"
)

# Species included in the model
rob_species <- rob_db %>%
  filter(!is.na(sp)) %>%
  distinct(sp) %>%
  arrange(sp) %>%
  pull(sp)

# Sequence of observed moonlight values on the scaled scale
moon_seq <- seq(
  min(rob_db$avg_moonlight_s, na.rm = TRUE),
  max(rob_db$avg_moonlight_s, na.rm = TRUE),
  length.out = 100
)

# Prediction grid
pred_grid_moon <- tidyr::expand_grid(
  avg_moonlight_s = moon_seq,
  trmt_bin = c(-1, 1),
  sp = rob_species
) %>%
  mutate(
    site = NA,
    jday_s = 0,
    nit_avg_wspm_s_s = 0,
    elev_max_s = 0,
    yr_s = 0,
    t_leps_s = 0,
    effort_hours = typical_effort_hours,
    log_background_offset = typical_log_background
  )



pred_moon_comm <- avg_predictions(
  m10_bg,
  newdata = pred_grid_moon,
  by = c("avg_moonlight_s", "trmt_bin"),
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    )
  )


p_moon_comm <- ggplot(
  pred_moon_comm,
  aes(x = avg_moonlight_s, y = estimate, color = treatment, fill = treatment)
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = trt_cols) +
  scale_fill_manual(values = trt_cols) +
  labs(
    title = "Moonlight modifies the effect of artificial light on moth decoy passes",
    subtitle = "Predictions from m10_bg; non-focal covariates held at mean scaled values",
    x = "Moonlight intensity (scaled)",
    y = "Predicted number of faint bat passes",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

p_moon_comm


#low moon  = 10th percentile
#mean moon = 50th percentile / median
#high moon = 90th percentile


moon_levels <- quantile(
  rob_db$avg_moonlight_s,
  probs = c(0.10, 0.50, 0.90),
  na.rm = TRUE
)

moon_level_df <- tibble(
  moon_level = c("Low moonlight", "Average moonlight", "High moonlight"),
  avg_moonlight_s = as.numeric(moon_levels)
)

moon_level_df

pred_grid_moon_levels <- tidyr::expand_grid(
  moon_level = moon_level_df$moon_level,
  trmt_bin = c(-1, 1),
  sp = rob_species
) %>%
  left_join(moon_level_df, by = "moon_level") %>%
  mutate(
    moon_level = factor(
      moon_level,
      levels = c("Low moonlight", "Average moonlight", "High moonlight")
    ),
    site = NA,
    jday_s = 0,
    nit_avg_wspm_s_s = 0,
    elev_max_s = 0,
    yr_s = 0,
    t_leps_s = 0,
    effort_hours = typical_effort_hours,
    log_background_offset = typical_log_background
    
  )

pred_moon_levels <- avg_predictions(
  m10_bg,
  newdata = pred_grid_moon_levels,
  by = c("moon_level", "avg_moonlight_s", "trmt_bin"),
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    ),
    moon_level = factor(
      moon_level,
      levels = c("Low moonlight", "Average moonlight", "High moonlight")
    )
  )

p_moon_levels <- ggplot(
  pred_moon_levels,
  aes(x = moon_level, y = estimate, color = treatment, group = treatment)
) +
  geom_line(linewidth = 1.1, position = position_dodge(width = 0.25)) +
  geom_pointrange(
    aes(ymin = conf.low, ymax = conf.high),
    linewidth = 0.7,
    fatten = 2.3,
    position = position_dodge(width = 0.25)
  ) +
  scale_color_manual(values = trt_cols) +
  labs(
    title = "Predicted faint robomoth passes at low and high moonlight",
    subtitle = "Predictions from m11; moonlight set to 10th, 50th, and 90th percentiles",
    x = "Moonlight level",
    y = "Predicted number of faint bat passes",
    color = "Treatment"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

p_moon_levels


# lepidoptera robmoth marginal effect  ------------------------------------



# Sequence of observed Lepidoptera values on the scaled scale
leps_seq <- seq(
  min(rob_db$t_leps_s, na.rm = TRUE),
  max(rob_db$t_leps_s, na.rm = TRUE),
  length.out = 100
)

# Prediction grid
pred_grid_leps <- tidyr::expand_grid(
  t_leps_s = leps_seq,
  trmt_bin = c(-1, 1),
  sp = rob_species
) %>%
  mutate(
    site = NA,
    jday_s = 0,
    nit_avg_wspm_s_s = 0,
    avg_moonlight_s = 0,
    elev_max_s = 0,
    yr_s = 0,
    
    # Required offsets
    effort_hours = typical_effort_hours,
    log_background_offset = typical_log_background
  )

# Community-average predictions across species
pred_leps_comm <- avg_predictions(
  m10_bg,
  newdata = pred_grid_leps,
  by = c("t_leps_s", "trmt_bin"),
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    )
  )

# Plot
p_leps_comm <- ggplot(
  pred_leps_comm,
  aes(x = t_leps_s, y = estimate, color = treatment, fill = treatment)
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = trt_cols) +
  scale_fill_manual(values = trt_cols) +
  labs(
    title = "Lepidoptera abundance modifies the effect of artificial light on faint robomoth passes",
    subtitle = paste0(
      "Predictions from m10_bg; non-focal covariates held at mean scaled values; ",
      "effort held at ", round(typical_effort_hours, 2), " h; ",
      "background activity held at ", round(typical_background, 1), " calls"
    ),
    x = "Total Lepidoptera abundance (scaled)",
    y = "Predicted robomoth bat passes",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

p_leps_comm



# spkr models.  -----------------------------------------------------------

# leaving zeros or including zeros

m0_s <- glmmTMB(
  n_calls ~ trmt_bin + (1 | site),
  family = poisson(link = "log"),
  data = spkr_db
)

# check model
check_singularity(m0)
check_zeroinflation(m0) # probably zero inflated. 
calculate_c_hat(m0)  # indicates the is overdispersion we need to try negative binomial. 
performance_mae(m0)
range(rob_db$n_calls)
hist(rob_db$n_calls, breaks = 50)
performance::r2(m0)
DHARMa::simulateResiduals(m0, n = 1000, plot = TRUE)

# including zeros 



# we will go straignt tot he negative binomial model. 
# phase 1: treatment-only model with offsets ------------------------------

m1_s_bg <- glmmTMB(
  n_calls ~ trmt_bin +
    offset(log(effort_hours)) +
    offset(log_background_offset) +
    (1 | site),
  family = nbinom2,
  data = spkr_db
)

check_singularity(m1_s_bg)
m1_s_bg$sdr$pdHess
m1_s_bg$fit$message
check_zeroinflation(m1_s_bg)
calculate_c_hat(m1_s_bg)
performance_mae(m1_s_bg)
performance::r2(m1_s_bg)
DHARMa::simulateResiduals(m1_s_bg, n = 1000, plot = TRUE)
summary(m1_s_bg)

# phase 2
# now we add seasonality.

m2_s_bg <- glmmTMB(
  n_calls ~ trmt_bin +
    jday_s + I(jday_s^2) +
    offset(log(effort_hours)) +
    offset(log_background_offset) +
    (1 | site),
  family = nbinom2,
  data = spkr_db
)

check_singularity(m2_s_bg)
m2_s_bg$sdr$pdHess
m2_s_bg$fit$message
check_zeroinflation(m2_s_bg)
calculate_c_hat(m2_s_bg)
performance_mae(m2_s_bg)
performance::r2(m2_s_bg)
DHARMa::simulateResiduals(m2_s_bg, n = 1000, plot = TRUE)
summary(m2_s_bg)

anova(m1_s_bg, m2_s_bg) # when we add seasonality the model improves 


#phase 3
# now we add environmental variables. treatment is not significant even with environmental variables.but adding environmental predictors does improve the model fit. 

# phase 3: add environmental predictors ----------------------------------

m3_s_bg <- glmmTMB(
  n_calls ~ trmt_bin +
    jday_s + I(jday_s^2) +
    nit_avg_wspm_s_s +
    avg_moonlight_s +
    elev_max_s +
    yr_s +
    offset(log(effort_hours)) +
    offset(log_background_offset) +
    (1 | site),
  family = nbinom2,
  data = spkr_db
)

check_singularity(m3_s_bg)
m3_s_bg$sdr$pdHess
m3_s_bg$fit$message
check_zeroinflation(m3_s_bg)
calculate_c_hat(m3_s_bg)
performance_mae(m3_s_bg)
performance::r2(m3_s_bg)
DHARMa::simulateResiduals(m3_s_bg, n = 1000, plot = TRUE)
performance::check_collinearity(m3_s_bg)
summary(m3_s_bg)

anova(m1_s_bg, m2_s_bg, m3_s_bg)

# phase 4 we add the insect data. with the insect data added we see an effect of moon and year but no treatment. 

# phase 4: add Lepidoptera abundance -------------------------------------
  
  m4_s_bg <- glmmTMB(
    n_calls ~ trmt_bin +
      jday_s + I(jday_s^2) +
      nit_avg_wspm_s_s +
      avg_moonlight_s +
      elev_max_s +
      yr_s +
      t_leps_s +
      offset(log(effort_hours)) +
      offset(log_background_offset) +
      (1 | site),
    family = nbinom2,
    data = spkr_db
  )

check_singularity(m4_s_bg)
m4_s_bg$sdr$pdHess
m4_s_bg$fit$message
check_zeroinflation(m4_s_bg)
calculate_c_hat(m4_s_bg)
performance_mae(m4_s_bg)
performance::r2(m4_s_bg)
DHARMa::simulateResiduals(m4_s_bg, n = 1000, plot = TRUE)
performance::check_collinearity(m4_s_bg)
summary(m4_s_bg)

anova(m1_s_bg, m2_s_bg, m3_s_bg, m4_s_bg) # adding insect data does not seem to improve model fit. 

# Phase 5 interactions leps and treatment. the interaction is not significant and the model does not improve. 
# 
# phase 5: treatment x Lepidoptera ---------------------------------------

m5_s_bg <- glmmTMB(
  n_calls ~ trmt_bin +
    jday_s + I(jday_s^2) +
    nit_avg_wspm_s_s +
    avg_moonlight_s +
    elev_max_s +
    yr_s +
    t_leps_s +
    offset(log(effort_hours)) +
    offset(log_background_offset) +
    (1 | site) +
    trmt_bin * t_leps_s,
  family = nbinom2,
  data = spkr_db
)

check_singularity(m5_s_bg)
m5_s_bg$sdr$pdHess
m5_s_bg$fit$message
check_zeroinflation(m5_s_bg)
calculate_c_hat(m5_s_bg)
performance_mae(m5_s_bg)
performance::r2(m5_s_bg)
DHARMa::simulateResiduals(m5_s_bg, n = 1000, plot = TRUE)
performance::check_collinearity(m5_s_bg)
summary(m5_s_bg)

anova(m4_s_bg, m5_s_bg) # the interaction is not significant neither so we might need to remove it. 

# Phase 6 lets add random slopes by species. 

# phase 6: add species random intercept ----------------------------------

m6_s_bg <- glmmTMB(
  n_calls ~ trmt_bin +
    jday_s + I(jday_s^2) +
    nit_avg_wspm_s_s +
    avg_moonlight_s +
    elev_max_s +
    yr_s +
    t_leps_s +
    offset(log(effort_hours)) +
    offset(log_background_offset) +
    (1 | site) +
    (1 | sp) +
    trmt_bin * t_leps_s,
  family = nbinom2,
  data = spkr_db
)

check_singularity(m6_s_bg)
m6_s_bg$sdr$pdHess
m6_s_bg$fit$message
check_zeroinflation(m6_s_bg)
calculate_c_hat(m6_s_bg)
performance_mae(m6_s_bg)
performance::r2(m6_s_bg)
DHARMa::simulateResiduals(m6_s_bg, n = 1000, plot = TRUE)
performance::check_collinearity(m6_s_bg)
summary(m6_s_bg)

anova(m5_s_bg, m6_s_bg) # the model improves with the addition of the random intercept for species.

# phase 7 lets add random treatment inside the random slopes by specie. 

m7_s_bg <- glmmTMB(
  n_calls ~ trmt_bin +
    jday_s + I(jday_s^2) +
    nit_avg_wspm_s_s +
    avg_moonlight_s +
    elev_max_s +
    yr_s +
    t_leps_s +
    offset(log(effort_hours)) +
    offset(log_background_offset) +
    (1 | site) +
    (1 + trmt_bin | sp) +
    trmt_bin * t_leps_s,
  family = nbinom2,
  data = spkr_db
)

check_singularity(m7_s_bg)
m7_s_bg$sdr$pdHess
m7_s_bg$fit$message
check_zeroinflation(m7_s_bg)
calculate_c_hat(m7_s_bg)
performance_mae(m7_s_bg)
performance::r2(m7_s_bg)
DHARMa::simulateResiduals(m7_s_bg, n = 1000, plot = TRUE)
performance::check_collinearity(m7_s_bg)
summary(m7_s_bg)

anova(m6_s_bg, m7_s_bg) # the model improves with the addition of the treatment inside the random slopes for species.

# seasonality inside the random slope 
# the model breaks if we add seasonality inside the random slopes
# 
# m8_s<- glmmTMB(
#   n_calls ~ trmt_bin + jday_s + I(jday_s^2) + nit_avg_wspm_s_s + nit_avg_temp_c_s + avg_moonlight_s + elev_max_s + yr_s +
#     + t_leps_s +
#     (1 | site)+ (1 + trmt_bin + jday_s + I(jday_s^2)| sp) +
#     #interactions
#     trmt_bin*t_leps_s,
#   family = nbinom2,
#   data = spkr_db
# )
# # # #
# check_singularity(m8_s)
# m8_s$sdr$pdHess
# m8_s$fit$message
# check_zeroinflation(m8_s)
# calculate_c_hat(m8_s)
# performance_mae(m8_s)
# range(rob_db$c_buzz)
# hist(rob_db$c_buzz, breaks = 100 )
# performance::r2(m8_s)
# DHARMa::simulateResiduals(m8_s, n = 1000, plot = TRUE)
# summary(m8_s)
# 
# anova(m0_s, m2_s, m3_s, m4_s, m5_s, m6_s, m7_s, m8_s)

# let's check what happens if we check for interaction between treatment and seasonality. It seems model 9 is not better than m7_s

# phase 9: treatment x seasonality ---------------------------------------

m9_s_bg <- glmmTMB(
  n_calls ~ trmt_bin +
    jday_s + I(jday_s^2) +
    nit_avg_wspm_s_s +
    avg_moonlight_s +
    elev_max_s +
    yr_s +
    t_leps_s +
    offset(log(effort_hours)) +
    offset(log_background_offset) +
    (1 | site) +
    (1 + trmt_bin | sp) +
    trmt_bin * t_leps_s +
    trmt_bin * jday_s +
    trmt_bin * I(jday_s^2),
  family = nbinom2,
  data = spkr_db
)

check_singularity(m9_s_bg)
m9_s_bg$sdr$pdHess
m9_s_bg$fit$message
check_zeroinflation(m9_s_bg)
calculate_c_hat(m9_s_bg)
performance_mae(m9_s_bg)
performance::r2(m9_s_bg)
DHARMa::simulateResiduals(m9_s_bg, n = 1000, plot = TRUE)
performance::check_collinearity(m9_s_bg)
summary(m9_s_bg)

anova(m7_s_bg, m9_s_bg) # model does improve with the addition of seasonality interaction with treatment. 

#now we add the moon interaction with treatment. Again it seems we observed the model improve in m10_s. The interaction is significant and the model improves.a

# phase 10: treatment x Lepidoptera + seasonality + moonlight ------------- this seems to be the best model. 

m10_s_bg <- glmmTMB(
  n_calls ~ trmt_bin +
    jday_s + I(jday_s^2) +
    nit_avg_wspm_s_s +
    avg_moonlight_s +
    elev_max_s +
    yr_s +
    t_leps_s +
    offset(log(effort_hours)) +
    offset(log_background_offset) +
    (1 | site) +
    (1 + trmt_bin | sp) +
    trmt_bin * t_leps_s +
    trmt_bin * jday_s +
    trmt_bin * I(jday_s^2) +
    trmt_bin * avg_moonlight_s,
  family = nbinom2,
  data = spkr_db
)

check_singularity(m10_s_bg)
m10_s_bg$sdr$pdHess
m10_s_bg$fit$message
check_zeroinflation(m10_s_bg)
calculate_c_hat(m10_s_bg)
performance_mae(m10_s_bg)
performance::r2(m10_s_bg)
DHARMa::simulateResiduals(m10_s_bg, n = 1000, plot = TRUE)
performance::check_collinearity(m10_s_bg)
summary(m10_s_bg)

anova(m7_s_bg, m9_s_bg, m10_s_bg) # model does improve with the addition of moonlight interaction with treatment.


# now we create the table for the model. 

m10_s_bg_table <- gtsummary::tbl_regression(
  m10_s_bg,
  
  # FALSE reports coefficients on the model's log scale
  exponentiate = FALSE,
  conf.level = 0.95,
  
  label = list(
    trmt_bin ~ "Treatment (Lit vs. Dark)",
    jday_s ~ "Julian date",
    `I(jday_s^2)` ~ "Julian date²",
    nit_avg_wspm_s_s ~ "Average wind speed",
    avg_moonlight_s ~ "Average moonlight",
    elev_max_s ~ "Maximum elevation",
    yr_s ~ "Year",
    t_leps_s ~ "Lepidoptera abundance",
    `trmt_bin:t_leps_s` ~ "Treatment × Lepidoptera abundance",
    `trmt_bin:jday_s` ~ "Treatment × Julian date",
    `trmt_bin:I(jday_s^2)` ~ "Treatment × Julian date²",
    `trmt_bin:avg_moonlight_s` ~ "Treatment × moonlight"
  )
) %>%
  modify_header(
    label ~ "**Predictor**",
    estimate ~ "**β**",
    conf.low ~ "**95% CI**",
    p.value ~ "**p-value**"
  ) %>%
  modify_caption(
    "**Table S__. Negative binomial mixed-effects model of bat call activity during speaker-lure experiments.**"
  ) %>%
  bold_p(t = 0.05) %>%
  modify_footnote(
    estimate ~ "Coefficients are presented on the log scale. The model includes offsets for speaker sampling effort and species-specific background acoustic activity."
  )

m10_s_bg_table


# make ir a flex table 

m10_s_bg_flex <- m10_s_bg_table %>%
  gtsummary::as_flex_table() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::align(j = 1, align = "left", part = "all") %>%
  flextable::align(j = 2:4, align = "left", part = "all") %>%
  flextable::padding(padding = 4, part = "all") %>%
  flextable::autofit()

m10_s_bg_flex


flextable::save_as_image(
  x = m10_s_bg_flex,
  path = "figures/rob_spkr_model/v3/tables/m10_s_bg_table.png",
  zoom = 3,
  expand = 15
)

flextable::save_as_docx(
  "Table m10_s_bg" = m10_s_bg_flex,
  path = "figures/rob_spkr_model/v3/tables/m10_s_bg_table.docx"
)

# now we remove lepidopter and treatment interaction to see if it improves model fit. 
# I didn't show improvements so we keep m10_s_bg as the best model.

# phase 11: remove treatment x Lepidoptera -------------------------------

m11_s_bg <- glmmTMB(
  n_calls ~ trmt_bin +
    jday_s + I(jday_s^2) +
    nit_avg_wspm_s_s +
    avg_moonlight_s +
    elev_max_s +
    yr_s +
    t_leps_s +
    offset(log(effort_hours)) +
    offset(log_background_offset) +
    (1 | site) +
    (1 + trmt_bin | sp) +
    trmt_bin * jday_s +
    trmt_bin * I(jday_s^2) +
    trmt_bin * avg_moonlight_s,
  family = nbinom2,
  data = spkr_db
)

check_singularity(m11_s_bg)
m11_s_bg$sdr$pdHess
m11_s_bg$fit$message
check_zeroinflation(m11_s_bg)
calculate_c_hat(m11_s_bg)
performance_mae(m11_s_bg)
performance::r2(m11_s_bg)
DHARMa::simulateResiduals(m11_s_bg, n = 1000, plot = TRUE)
performance::check_collinearity(m11_s_bg)
summary(m11_s_bg)

anova( m10_s_bg, m11_s_bg) # it seems he best model might be m10_s_bg 

# m11_s_table <- gtsummary::tbl_regression(
#   m11_s,
#   
#   # FALSE reports coefficients on the model's log scale
#   exponentiate = FALSE,
#   
#   conf.level = 0.95,
#   
#   label = list(
#     trmt_bin ~ "Treatment (Lit vs. Dark)",
#     jday_s ~ "Julian date",
#     `I(jday_s^2)` ~ "Julian date²",
#     nit_avg_wspm_s_s ~ "Average wind speed",
#     avg_moonlight_s ~ "Average moonlight",
#     elev_max_s ~ "Maximum elevation",
#     yr_s ~ "Year",
#     t_leps_s ~ "Lepidoptera abundance",
#     `trmt_bin:jday_s` ~ "Treatment × Julian date",
#     `trmt_bin:I(jday_s^2)` ~ "Treatment × Julian date²",
#     `trmt_bin:avg_moonlight_s` ~ "Treatment × moonlight"
#   )
# ) %>%
#   
#   modify_header(
#     label ~ "**Predictor**",
#     estimate ~ "**β**",
#     ci ~ "**95% CI**",
#     p.value ~ "**p-value**"
#   ) %>%
#   
#   modify_caption(
#     "**Table S__. Negative binomial mixed-effects model of bat call activity during speaker experiments.**"
#   ) %>%
#   
#   bold_p(t = 0.05) %>%
#   
#   modify_footnote(
#     estimate ~
#       "Coefficients are presented on the log scale."
#   )
# 
# m11_s_table
# 
# 
# # convert into a flex table 
# 
# m11_s_flex <- m11_s_table %>%
#   gtsummary::as_flex_table() %>%
#   flextable::theme_booktabs() %>%
#   flextable::font(
#     fontname = "Times New Roman",
#     part = "all"
#   ) %>%
#   flextable::fontsize(
#     size = 10,
#     part = "all"
#   ) %>%
#   flextable::bold(
#     part = "header"
#   ) %>%
#   flextable::align(
#     j = 1,
#     align = "left",
#     part = "all"
#   ) %>%
#   flextable::align(
#     j = 2:ncol(m11_s_flex$body$dataset),
#     align = "center",
#     part = "all"
#   ) %>%
#   flextable::padding(
#     padding = 4,
#     part = "all"
#   ) %>%
#   flextable::autofit()
# 
# m11_s_flex
# 
# m11_s_flex <- m11_s_table %>%
#   gtsummary::as_flex_table() %>%
#   flextable::theme_booktabs() %>%
#   flextable::font(fontname = "Times New Roman", part = "all") %>%
#   flextable::fontsize(size = 10, part = "all") %>%
#   flextable::bold(part = "header") %>%
#   flextable::autofit()
# 
# flextable::save_as_image(
#   x = m11_s_flex,
#   path = "figures/rob_spkr_model/v2/tables/m11_s_table.png",
#   zoom = 3,
#   expand = 15
# )
# 
# 


# # we have been suggested to test a model with an offset lets see how it changes. 
# 
# m12_s <- glmmTMB(
#   n_calls ~ trmt_bin + jday_s + I(jday_s^2) +
#     nit_avg_wspm_s_s + avg_moonlight_s + elev_max_s +
#     yr_s + t_leps_s +
#     offset(log(effort_hours)) +
#     (1 | site) +
#     (1 + trmt_bin | sp) +
#     trmt_bin * jday_s +
#     trmt_bin * I(jday_s^2) +
#     trmt_bin * avg_moonlight_s,
#   family = nbinom2,
#   data = spkr_db
# )
# 
# check_singularity(m12_s)
# m12_s$sdr$pdHess
# m12_s$fit$message
# check_zeroinflation(m12_s)
# calculate_c_hat(m12_s)   
# performance_mae(m12_s)
# performance::r2(m12_s)
# DHARMa::simulateResiduals(m12_s, n = 1000, plot = TRUE)
# performance::check_collinearity(m12_s)
# summary(m12_s)
# anova( m10_s, m11_s, m12_s)
# 
# # now we generate marginal effects plot for the model m11_s for the sepeaker data. It seems we have a final model. 


# spkr marginal effect plots.  --------------------------------------------


# Treatment colors
trt_cols <- c(
  "Dark" = "grey20",
  "Lit"  = "grey"
)

# Species included in speaker model/data
spkr_species <- spkr_db %>%
  filter(!is.na(sp)) %>%
  distinct(sp) %>%
  arrange(sp) %>%
  pull(sp)

# typical values for offsets. 
typical_effort_hours <- median(
  spkr_db$effort_hours,
  na.rm = TRUE
)

typical_log_background <- median(
  spkr_db$log_background_offset,
  na.rm = TRUE
)

typical_background <- exp(typical_log_background)

typical_effort_hours
typical_background

# -------------------------------------------------------------------------
# Speaker treatment prediction grid
# -------------------------------------------------------------------------

pred_grid_spkr_trt <- tidyr::expand_grid(
  sp = spkr_species,
  trmt_bin = c(-1, 1)
) %>%
  mutate(
    site = NA,
    # Non-focal covariates held at mean scaled values
    jday_s = 0,
    nit_avg_wspm_s_s = 0,
    avg_moonlight_s = 0,
    elev_max_s = 0,
    yr_s = 0,
    t_leps_s = 0,
    
    # Required offsets for m10_s_bg
    effort_hours = typical_effort_hours,
    log_background_offset = typical_log_background
  )

# -------------------------------------------------------------------------
# Species-specific treatment predictions
# -------------------------------------------------------------------------

pred_spkr_species <- predictions(
  m10_s_bg,
  newdata = pred_grid_spkr_trt,
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    ),
    panel = as.character(sp)
  )

# -------------------------------------------------------------------------
# Community-average treatment predictions
# -------------------------------------------------------------------------

pred_spkr_comm <- avg_predictions(
  m10_s_bg,
  newdata = pred_grid_spkr_trt,
  by = "trmt_bin",
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    ),
    panel = "Community"
  )


# -------------------------------------------------------------------------
# Combine community and species predictions
# -------------------------------------------------------------------------

pred_spkr_treatment_all <- bind_rows(
  pred_spkr_comm,
  pred_spkr_species
) %>%
  mutate(
    panel = factor(panel, levels = c("Community", spkr_species))
  )



# -------------------------------------------------------------------------
# Speaker treatment marginal effects plot
# -------------------------------------------------------------------------

p_spkr_treatment <- ggplot(
  pred_spkr_treatment_all,
  aes(x = treatment, y = estimate)
) +
  geom_line(
    aes(group = panel),
    color = "black",
    linewidth = 0.75
  ) +
  geom_pointrange(
    aes(
      ymin = conf.low,
      ymax = conf.high,
      color = treatment
    ),
    linewidth = 0.7,
    fatten = 2.4
  ) +
  facet_wrap(
    ~ panel,
    scales = "free_y"
  ) +
  scale_color_manual(values = trt_cols) +
  labs(
    title = "Predicted faint bat passes at speaker lures by treatment",
    subtitle = "Predictions from m10_s_bg; non-focal covariates held at mean values",
    x = "Treatment",
    y = "Predicted number of faint bat passes",
    color = "Treatment"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

p_spkr_treatment



spkr_ranef_trt <- ranef(m10_s_bg)$cond$sp %>%
  as.data.frame() %>%
  rownames_to_column("sp") %>%
  rename(
    species_intercept_deviation = `(Intercept)`,
    species_trmt_deviation = trmt_bin
  )

fixed_trmt <- fixef(m10_s_bg)$cond["trmt_bin"]

spkr_species_trt_effects <- spkr_ranef_trt %>%
  mutate(
    fixed_trmt = as.numeric(fixed_trmt),
    species_trmt_half_contrast = fixed_trmt + species_trmt_deviation,
    lit_dark_log_ratio = 2 * species_trmt_half_contrast,
    lit_dark_ratio = exp(lit_dark_log_ratio),
    percent_change_lit_vs_dark = 100 * (lit_dark_ratio - 1)
  ) %>%
  arrange(percent_change_lit_vs_dark)

spkr_species_trt_effects



spkr_species_trt_effects %>%
  mutate(
    sp = factor(sp, levels = sp),
    response = if_else(percent_change_lit_vs_dark > 0,
                       "Higher in lit",
                       "Lower in lit")
  ) %>%
  ggplot(aes(x = percent_change_lit_vs_dark, y = sp)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7) +
  geom_point(size = 3) +
  labs(
    title = "Species-specific treatment response at speaker lures",
    subtitle = "Percent change in predicted faint passes from dark to lit sites",
    x = "Percent change in predicted faint passes: lit vs dark",
    y = "Species"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )


# uncertainty table 

library(tidyverse)
library(marginaleffects)

# Species in robomoth data
spkr_species <- spkr_db %>%
  filter(!is.na(sp)) %>%
  distinct(sp) %>%
  arrange(sp) %>%
  pull(sp)

# Prediction grid at mean covariate values
spkr_trt_grid <- tidyr::expand_grid(
  sp = spkr_species,
  trmt_bin = -1
) %>%
  mutate(
    site = NA,
    jday_s = 0,
    nit_avg_wspm_s_s = 0,
    avg_moonlight_s = 0,
    elev_max_s = 0,
    yr_s = 0,
    t_leps_s = 0,
    
    # Required offsets
    effort_hours = typical_effort_hours,
    log_background_offset = typical_log_background
  )



# Species-specific lit/dark ratios
spkr_species_ratio <- avg_comparisons(
  m10_s_bg,
  newdata = spkr_trt_grid,
  variables = list(trmt_bin = c(-1, 1)),
  comparison = "ratio",
  by = "sp",
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    percent_change_lit_vs_dark = 100 * (estimate - 1),
    percent_conf_low = 100 * (conf.low - 1),
    percent_conf_high = 100 * (conf.high - 1),
    supported_response = case_when(
      conf.low > 1 ~ "Higher in lit",
      conf.high < 1 ~ "Lower in lit",
      TRUE ~ "Uncertain"
    )
  ) %>%
  arrange(desc(percent_change_lit_vs_dark))

spkr_species_ratio

# clean uncertainty table. 

spkr_species_ratio_clean <- spkr_species_ratio %>%
  transmute(
    species = sp,
    lit_dark_ratio = round(estimate, 2),
    ratio_CI = paste0(round(conf.low, 2), "–", round(conf.high, 2)),
    percent_change = round(percent_change_lit_vs_dark, 1),
    percent_CI = paste0(
      round(percent_conf_low, 1),
      " to ",
      round(percent_conf_high, 1),
      "%"
    ),
    p_value = signif(p.value, 3),
    supported_response
  )

spkr_species_ratio_clean


ts5<-flextable(spkr_species_ratio_clean) %>%
  autofit() %>%
  set_caption(
    ""
  )

ts5

ts5 <- flextable(spkr_species_ratio_clean) %>%
  theme_booktabs() %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  bold(
    i = ~ supported_response != "Uncertain",
    bold = FALSE,
    part = "body"
  ) %>%
  autofit() %>%
  set_caption(
    "Table S__. Species-specific lit/dark ratios for faint bat passes at speaker lures from the effort- and background-offset speaker model."
  )

ts5

save_as_docx(
  "table.s2" = ts5,
  path = "figures/rob_spkr_model/v3/tables/speaker_uncertainty.docx")

save_as_image(ts5, path = "figures/rob_spkr_model/v3/tables/spkr_uncertainty.png")



# improved table 

# -------------------------------------------------------------------------
# Species-specific lit/dark log ratios with uncertainty
# -------------------------------------------------------------------------

spkr_species_logratio <- avg_comparisons(
  m10_s_bg,
  newdata = spkr_trt_grid,
  variables = list(trmt_bin = c(-1, 1)),
  comparison = "lnratio",
  by = "sp",
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    # Convert log ratio back to lit/dark ratio
    lit_dark_ratio = exp(estimate),
    ratio_conf_low = exp(conf.low),
    ratio_conf_high = exp(conf.high),
    
    # Convert ratio to percent change
    percent_change_lit_vs_dark = 100 * (lit_dark_ratio - 1),
    percent_conf_low = 100 * (ratio_conf_low - 1),
    percent_conf_high = 100 * (ratio_conf_high - 1),
    
    # On the log-ratio scale, 0 means no difference.
    # If the CI is entirely below 0, lit is lower than dark.
    # If the CI is entirely above 0, lit is higher than dark.
    supported_response = case_when(
      conf.low > 0 ~ "Higher in lit",
      conf.high < 0 ~ "Lower in lit",
      TRUE ~ "Uncertain"
    )
  ) %>%
  arrange(desc(percent_change_lit_vs_dark))

spkr_species_logratio


# clean table 

spkr_species_ratio_clean <- spkr_species_logratio %>%
  transmute(
    species = sp,
    lit_dark_ratio = round(lit_dark_ratio, 2),
    ratio_CI = paste0(
      round(ratio_conf_low, 2),
      "–",
      round(ratio_conf_high, 2)
    ),
    percent_change = round(percent_change_lit_vs_dark, 1),
    percent_CI = paste0(
      round(percent_conf_low, 1),
      " to ",
      round(percent_conf_high, 1),
      "%"
    ),
    p_value = signif(p.value, 3),
    supported_response
  )

spkr_species_ratio_clean

ts5<-flextable(spkr_species_ratio_clean) %>%
  autofit() %>%
  set_caption(
    ""
  )

ts5

ts5 <- flextable(spkr_species_ratio_clean) %>%
  theme_booktabs() %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  bold(
    i = ~ supported_response != "Uncertain",
    bold = TRUE,
    part = "body"
  ) %>%
  autofit() %>%
  set_caption(
    "Table S__. Species-specific lit/dark ratios for faint bat passes at speaker lures from the effort- and background-offset speaker model."
  )

ts5

save_as_docx(
  "table.s2" = ts5,
  path = "figures/rob_spkr_model/v3/tables/speaker_uncertainty.docx")

save_as_image(ts5, path = "figures/rob_spkr_model/v3/tables/spkr_uncertainty.png")




# speaker Jday marginal effects  ------------------------------------------

# Treatment colors
trt_cols <- c(
  "Dark" = "grey20",
  "Lit"  = "grey"
)

# Get scaling values from your data
# This assumes jday_s was scaled as: (jday - mean(jday)) / (2 * sd(jday))
jday_mean <- mean(spkr_db$jday, na.rm = TRUE)
jday_sd   <- sd(spkr_db$jday, na.rm = TRUE)

# Raw Julian day sequence for plotting
jday_seq <- seq(
  min(spkr_db$jday, na.rm = TRUE),
  max(spkr_db$jday, na.rm = TRUE),
  length.out = 100
)

# Species included in model/data
rob_species <- spkr_db %>%
  filter(!is.na(sp)) %>%
  distinct(sp) %>%
  arrange(sp) %>%
  pull(sp)

# Prediction grid
pred_grid_jday <- tidyr::expand_grid(
  jday = jday_seq,
  trmt_bin = c(-1, 1),
  sp = rob_species
) %>%
  mutate(
    jday_s = (jday - jday_mean) / (2 * jday_sd),
    site = NA,
    nit_avg_wspm_s_s = 0,
    avg_moonlight_s = 0,
    elev_max_s = 0,
    yr_s = 0,
    t_leps_s = 0,
    # Required offsets for m10_s_bg
    effort_hours = typical_effort_hours,
    log_background_offset = typical_log_background
    
  )


pred_jday_comm <- avg_predictions(
  m10_s_bg,
  newdata = pred_grid_jday,
  by = c("jday", "jday_s", "trmt_bin"),
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    )
  )


p_jday_comm <- ggplot(
  pred_jday_comm,
  aes(x = jday, y = estimate, color = treatment, fill = treatment)
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = trt_cols) +
  scale_fill_manual(values = trt_cols) +
  labs(
    title = "Seasonal pattern of faint robomoth bat passes by treatment",
    subtitle = "Predictions from m11; non-focal covariates held at mean scaled values",
    x = "Julian day",
    y = "Predicted number of faint bat passes",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

p_jday_comm


raw_jday_summary <- spkr_db %>%
  group_by(treatmt, trmt_bin, jday) %>%
  summarise(
    mean_calls = mean(n_calls, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    )
  )



p_jday_comm_raw_spkr <- ggplot(
  pred_jday_comm,
  aes(x = jday, y = estimate, color = treatment, fill = treatment)
) +
  geom_point(
    data = raw_jday_summary,
    aes(x = jday, y = mean_calls, color = treatment),
    alpha = 0.25,
    size = 1.4,
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = trt_cols) +
  scale_fill_manual(values = trt_cols) +
  labs(
    # title = "Seasonal pattern of faint acoustic lure bat passes by treatment",
    # subtitle = "Points are observed treatment-level means; lines are model predictions from m11_s",
    # x = "Julian day",
     y = "",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

p_jday_comm_raw_spkr

ggsave(
  p_jday_comm_raw_spkr,
  filename = "figures/rob_spkr_model/v3/spkr_jday_marginal_effects.png",
  width = 8,
  height = 6,
  dpi = 300)





# lepidoptera spkr marginal effects plot ----------------------------------


# Sequence of observed Lepidoptera values on the scaled scale
leps_seq <- seq(
  quantile(spkr_db$t_leps_s, 0.05, na.rm = TRUE),
  quantile(spkr_db$t_leps_s, 0.95, na.rm = TRUE),
  length.out = 100
)

# Prediction grid
pred_grid_leps <- tidyr::expand_grid(
  t_leps_s = leps_seq,
  trmt_bin = c(-1, 1),
  sp = spkr_species
) %>%
  mutate(
    site = NA,
    jday_s = 0,
    nit_avg_wspm_s_s = 0,
    avg_moonlight_s = 0,
    elev_max_s = 0,
    yr_s = 0,
    
    # Required offsets
    effort_hours = typical_effort_hours,
    log_background_offset = typical_log_background
  )

# Community-average predictions across species
pred_leps_comm <- avg_predictions(
  m10_s_bg,
  newdata = pred_grid_leps,
  by = c("t_leps_s", "trmt_bin"),
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    )
  )

# Plot
p_leps_comm_spkr <- ggplot(
  pred_leps_comm,
  aes(x = t_leps_s, y = estimate, color = treatment, fill = treatment)
) +
  geom_ribbon(
    aes(ymin = pmax(conf.low,0), ymax = conf.high),
    alpha = 0.08,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = trt_cols) +
  scale_fill_manual(values = trt_cols) +
  labs(
    title = "Lepidoptera abundance modifies the effect of artificial light on faint spkr passes",
    subtitle = paste0(
      "Predictions from m10_s_bg; non-focal covariates held at mean scaled values; ",
      "effort held at ", round(typical_effort_hours, 2), " h; ",
      "background activity held at ", round(typical_background, 1), " calls"
    ),
    x = "Total Lepidoptera abundance (scaled)",
    y = "Predicted number of bat acoustic lure passes",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

p_leps_comm_spkr





# both robomoth and speaker data ------------------------------------------

# color for the two decoys and lures 
assay_cols <- c(
  "Speaker"  = "grey45",  # teal
  "Robomoth" = "black"   # red/salmon
)

# ------------------------------------------------------------
# Helper function to build prediction grid
# ------------------------------------------------------------

make_trt_grid <- function(dat) {
  
  species <- dat %>%
    filter(!is.na(sp)) %>%
    distinct(sp) %>%
    arrange(sp) %>%
    pull(sp)
  
  typical_effort_hours <- median(
    dat$effort_hours,
    na.rm = TRUE
  )
  
  typical_log_background <- median(
    dat$log_background_offset,
    na.rm = TRUE
  )
  
  grid <- tidyr::expand_grid(
    sp = species,
    trmt_bin = c(-1, 1)
  ) %>%
    mutate(
      site = NA,
      
      # Non-focal scaled covariates held at mean/reference value
      jday_s = 0,
      nit_avg_wspm_s_s = 0,
      avg_moonlight_s = 0,
      elev_max_s = 0,
      yr_s = 0,
      t_leps_s = 0,
      
      # Required offsets for current models
      effort_hours = typical_effort_hours,
      log_background_offset = typical_log_background
    )
  
  if (is.factor(dat$sp)) {
    grid$sp <- factor(grid$sp, levels = levels(dat$sp))
  }
  
  if (is.factor(dat$site)) {
    grid$site <- factor(NA, levels = levels(dat$site))
  }
  
  return(grid)
}


# robomoth predictions 
# 
# 
rob_grid <- make_trt_grid(rob_db)

pred_rob_sp <- predictions(
  m10_bg,
  newdata = rob_grid,
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    assay = "Robomoth",
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    ),
    panel = as.character(sp)
  )

pred_rob_comm <- avg_predictions(
  m10_bg,
  newdata = rob_grid,
  by = "trmt_bin",
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    assay = "Robomoth",
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    ),
    panel = "Community"
  )



# spkr predictions

spkr_grid <- make_trt_grid(spkr_db)

pred_spkr_sp <- predictions(
  m10_s_bg,
  newdata = spkr_grid,
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    assay = "Speaker",
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    ),
    panel = as.character(sp)
  )

pred_spkr_comm <- avg_predictions(
  m10_s_bg,
  newdata = spkr_grid,
  by = "trmt_bin",
  type = "response",
  re.form = NULL,
  allow.new.levels = TRUE
) %>%
  as_tibble() %>%
  mutate(
    assay = "Speaker",
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    ),
    panel = "Community"
  )

#combine predictions 
pred_compare <- bind_rows(
  pred_rob_comm,
  pred_rob_sp,
  pred_spkr_comm,
  pred_spkr_sp
) %>%
  mutate(
    assay = factor(assay, levels = c("Speaker", "Robomoth")),
    treatment = factor(treatment, levels = c("Dark", "Lit")),
    panel = as.character(panel)
  )

# Keep only species present in both assays, plus Community
shared_panels <- intersect(
  unique(as.character(pred_rob_sp$panel)),
  unique(as.character(pred_spkr_sp$panel))
) %>% 
  sort()

panel_levels <- c("Community", shared_panels)

pred_compare <- pred_compare %>%
  filter(panel %in% panel_levels) %>%
  mutate(
    panel = factor(panel, levels = panel_levels)
  )


raw_rob <- rob_db %>%
  filter(
    !is.na(n_calls),
    !is.na(sp),
    !is.na(trmt_bin),
    !is.na(effort_hours),
    effort_hours > 0,
    !is.na(log_background_offset)
  ) %>%
  mutate(
    assay = "Robomoth",
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    ),
    panel = as.character(sp)
  )

raw_spkr <- spkr_db %>%
  filter(
    !is.na(n_calls),
    !is.na(sp),
    !is.na(trmt_bin),
    !is.na(effort_hours),
    effort_hours > 0,
    !is.na(log_background_offset)
  ) %>%
  mutate(
    assay = "Speaker",
    treatment = factor(
      if_else(trmt_bin == -1, "Dark", "Lit"),
      levels = c("Dark", "Lit")
    ),
    panel = as.character(sp)
  )

raw_species <- bind_rows(raw_rob, raw_spkr) %>%
  filter(panel %in% shared_panels)

# community level raw observations
raw_comm <- bind_rows(raw_rob, raw_spkr) %>%
  group_by(assay, treatment, trmt_bin, site, noche) %>%
  summarise(
    n_calls = sum(n_calls, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(panel = "Community")

raw_compare <- bind_rows(raw_comm, raw_species) %>%
  mutate(
    assay = factor(assay, levels = c("Speaker", "Robomoth")),
    treatment = factor(treatment, levels = c("Dark", "Lit")),
    panel = factor(panel, levels = panel_levels)
  )

# color for the two decoys and lures 
assay_cols <- c(
  "Speaker"  = "green2",  # teal
  "Robomoth" = "black"   # red/salmon
)

# now we plot both analysis. 

pd <- position_dodge(width = 0.18)

 
# raw observations
p_gomes_style_treatment <- ggplot() +
  geom_jitter(
    data = raw_compare,
    aes(x = treatment, y = n_calls, color = assay),
    width = 0.08,
    height = 0,
    alpha = 0.10,
    size = 0.7
  ) +
  # model estimates
  geom_errorbar(
    data = pred_compare,
    aes(
      x = treatment,
      ymin = conf.low,
      ymax = conf.high,
      color = assay
    ),
    width = 0.08,
    linewidth = 0.5,
    position = pd) +
  
  # geom_line(
  #   data = pred_compare,
  #   aes(
  #     x = treatment,
  #     y = estimate,
  #     group = assay,
  #     color = assay
  #   ),
  #   linewidth = 0.9,
  #   position = pd
  # ) +
  # model -estimated predicted mean
  geom_point(
    data = pred_compare,
    aes(
      x = treatment,
      y = estimate,
      color = assay
    ),
    size = 2.2,
    position = pd
  ) +
  # geom_ribbon(
  #   data = pred_compare,
  #   aes(
  #     x = treatment,
  #     ymin = conf.low,
  #     ymax = conf.high,
  #     fill = assay
  #   ),
  #   alpha = 0.2,
  #   position = pd
  # ) +
  facet_wrap(~ panel, scales = "free_y", ncol = 4) +
  scale_color_manual(values = assay_cols) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10)
  ) +
  labs(
    title = "Contrasting speaker and robomoth responses to artificial light",
    subtitle = paste(
      "Predictions from selected negative-binomial mixed models;",
      "models include offsets for sampling effort and species-specific background activity;",
      "points show observed faint passes"
    ),
    x = "Treatment",
    y = "Bat passes",
    color = "experiment"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_gomes_style_treatment


pd <- position_dodge(width = 0.25)

p_gomes_style_treatment <- ggplot() +
  geom_jitter(
    data = raw_compare,
    aes(x = treatment, y = n_calls, color = assay),
    width = 0.08,
    height = 0,
    alpha = .4,
    size = 1
  ) +
  geom_line(
    data = pred_compare,
    aes(
      x = treatment,
      y = estimate,
      group = assay,
      color = assay
    ),
    linewidth = 1,
    position = pd
  ) +
  # geom_pointrange(
  #   data = pred_compare,
  #   aes(
  #     x = treatment,
  #     y = estimate,
  #     ymin = conf.low,
  #     ymax = conf.high,
  #     color = assay
  #   ),
  #   linewidth = 1,
  #   fatten = 2.2,
  #   position = pd
  # ) +
  facet_wrap(~ panel, scales = "free_y", ncol = 4) +
  scale_color_manual(values = assay_cols) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10)
  ) +
  labs(
    title = "Contrasting speaker and robomoth responses to artificial light",
    subtitle = "Predictions from selected negative-binomial mixed models; points show observed faint passes",
    x = "Treatment",
    y = "Faint bat passes",
    color = "Assay"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_gomes_style_treatment

ggsave(
  filename = "figures/rob_spkr_model/v3/rob_spkr_treatment.tiff",)



# trial for improved graph 

plot_ratio_df <- bind_rows(
  rob_species_ratio_clean %>%
    mutate(assay = "Robomoth"),
  
  spkr_species_ratio_clean %>%
    mutate(assay = "Speaker")
) %>%
  mutate(
    assay = factor(assay, levels = c("Speaker", "Robomoth")),
    
    # Extract lower and upper percent confidence limits from text.
    # Example: "-96.1 to 110.7%" becomes -96.1 and 110.7
    percent_conf_low = map_dbl(
      str_extract_all(percent_CI, "-?\\d+\\.?\\d*"),
      ~ as.numeric(.x[1])
    ),
    percent_conf_high = map_dbl(
      str_extract_all(percent_CI, "-?\\d+\\.?\\d*"),
      ~ as.numeric(.x[2])
    ),
    
    support_status = if_else(
      supported_response == "Uncertain",
      "Uncertain",
      "Supported evidence"
    )
  )

plot_ratio_df <- plot_ratio_df %>%
  group_by(assay) %>%
  arrange(percent_change, .by_group = TRUE) %>%
  mutate(
    species_facet = paste0(assay, "__", species)
  ) %>%
  ungroup()

plot_ratio_df$species_facet <- factor(
  plot_ratio_df$species_facet,
  levels = unique(plot_ratio_df$species_facet)
)


# -------------------------------------------------------------------------

support_cols <- c(
  "Supported evidence" = "black",
  "Uncertain" = "grey70"
)

p_species_support <- ggplot(
  plot_ratio_df,
  aes(
    x = percent_change,
    y = species_facet,
    color = support_status
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey40"
  ) +
  geom_segment(
    aes(
      x = percent_conf_low,
      xend = percent_conf_high,
      y = species_facet,
      yend = species_facet
    ),
    linewidth = 0.8
  ) +
  geom_point(size = 3) +
  facet_wrap(
    ~ assay,
    scales = "free_y",
    ncol = 2
  ) +
  scale_y_discrete(
    labels = function(x) sub(".*__", "", x)
  ) +
  scale_color_manual(values = support_cols) +
  labs(
    title = "Species-specific responses to artificial light",
    subtitle = "Black points indicate species with confidence intervals excluding no change",
    x = "Percent change in lit relative to dark sites",
    y = "Species",
    color = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_species_support

ggsave(
  filename = "figures/rob_spkr_model/v3/support_spkr_rob.png",
  plot = p_species_support,
  width = 9,
  height = 6,
  dpi = 300
)

#Jday with patchwork for jday for both speaker and robomoth data. 

p_jday_robomoth_spkr <-
  (p_jday_comm_raw | p_jday_comm_raw_spkr) +
  plot_annotation(
    title = "Faint bat passes at moth decoys (A) and speaker lures (B) by treatment and Julian day",
    subtitle = paste(
      "Lines and points show model call model predictions",
      "shaded areas show 95% confidence intervals"
    ),
    tag_levels = "A"
  ) &
  theme(
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold")
  )

p_jday_robomoth_spkr


ggsave(
  plot = p_jday_robomoth_spkr,
  filename = "figures/rob_spkr_model/v3/jday_spkr_robomoth.tiff",
  width = 12,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw", # this compression just works with tiff files. 
  bg = "white"
)


# lepidoptea with patwork

leps_robomoth_spkr <-
  (
    p_leps_comm +
      labs(title = NULL, subtitle = NULL) +
      labs(y = "Predicted number of faint bat passes")
    |
      p_leps_comm_spkr +
      labs(title = NULL, subtitle = NULL) +
      labs(y = NULL)
  ) +
  plot_annotation(
    title = "Bat calls at moth decoys (A) and speaker lures (B) by treatment and Lepidoptera abundance",
    subtitle = paste(
      "Lines show model predictions;",
      "shaded areas show 95% confidence intervals"
    ),
    tag_levels = "A"
  ) &
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    plot.tag = element_text(face = "bold")
  )

leps_robomoth_spkr

ggsave(
  plot = leps_robomoth_spkr,
  filename = "figures/rob_spkr_model/v3/leps_spkr_robomoth.tiff",
  width = 12,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw", # this compression just works with tiff files. 
  bg = "white"
)

# moon marginal effects for both robomoth and speakr moon
moon_predictions <- function(model, data, experiment_label) {
  
  moon_seq <- seq(
    min(data$avg_moonlight_s, na.rm = TRUE),
    max(data$avg_moonlight_s, na.rm = TRUE),
    length.out = 100
  )
  
  species <- data %>%
    filter(!is.na(sp)) %>%
    distinct(sp) %>%
    arrange(sp) %>%
    pull(sp)
  
  typical_effort_hours <- median(
    data$effort_hours,
    na.rm = TRUE
  )
  
  typical_log_background <- median(
    data$log_background_offset,
    na.rm = TRUE
  )
  
  typical_background <- exp(typical_log_background)
  
  new_data <- tidyr::crossing(
    avg_moonlight_s = moon_seq,
    trmt_bin = c(-1, 1),
    sp = species
  ) %>%
    mutate(
      site = NA,
      jday_s = 0,
      nit_avg_wspm_s_s = 0,
      elev_max_s = 0,
      yr_s = 0,
      t_leps_s = 0,
      effort_hours = typical_effort_hours,
      log_background_offset = typical_log_background
    )
  
  if (is.factor(data$sp)) {
    new_data$sp <- factor(new_data$sp, levels = levels(data$sp))
  }
  
  if (is.factor(data$site)) {
    new_data$site <- factor(NA, levels = levels(data$site))
  }
  
  marginaleffects::avg_predictions(
    model,
    newdata = new_data,
    by = c("avg_moonlight_s", "trmt_bin"),
    type = "response",
    re.form = NULL,
    allow.new.levels = TRUE
  ) %>%
    as_tibble() %>%
    mutate(
      Treatment = factor(
        trmt_bin,
        levels = c(-1, 1),
        labels = c("Dark", "Lit")
      ),
      Experiment = experiment_label,
      typical_effort_hours = typical_effort_hours,
      typical_background = typical_background
    ) %>%
    arrange(Treatment, avg_moonlight_s)
}

# -------------------------------------------------------------------------
# Generate predictions from current selected models
# -------------------------------------------------------------------------

pred_moon_rob <- moon_predictions(
  model = m10_bg,
  data = rob_db,
  experiment_label = "Moth decoys"
)

pred_moon_spkr <- moon_predictions(
  model = m10_s_bg,
  data = spkr_db,
  experiment_label = "Speaker lures"
)

# -------------------------------------------------------------------------
# Check prediction ranges
# -------------------------------------------------------------------------

bind_rows(pred_moon_rob, pred_moon_spkr) %>%
  group_by(Experiment) %>%
  summarise(
    minimum_estimate = min(estimate, na.rm = TRUE),
    minimum_lower_CI = min(conf.low, na.rm = TRUE),
    maximum_upper_CI = max(conf.high, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------------------------------------------------
# Plotting function
# -------------------------------------------------------------------------

plot_moon_effect <- function(pred_data, panel_title, y_title = NULL) {
  
  ggplot(
    pred_data,
    aes(
      x = avg_moonlight_s,
      y = estimate,
      color = Treatment,
      fill = Treatment,
      group = Treatment
    )
  ) +
    geom_ribbon(
      aes(
        ymin = conf.low,
        ymax = conf.high
      ),
      alpha = 0.20,
      color = NA
    ) +
    geom_line(linewidth = 1.15) +
    scale_color_manual(
      values = c(
        "Dark" = "grey20",
        "Lit"  = "grey60"
      )
    ) +
    scale_fill_manual(
      values = c(
        "Dark" = "grey20",
        "Lit"  = "grey60"
      )
    ) +
    scale_y_continuous(
      trans = "log1p"
    ) +
    labs(
      x = "Moonlight intensity (standardized)",
      y = y_title,
      title = panel_title,
      color = "Treatment",
      fill = "Treatment"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

# -------------------------------------------------------------------------
# Build individual panels
# -------------------------------------------------------------------------

p_moon_rob <- plot_moon_effect(
  pred_moon_rob,
  panel_title = "A",
  y_title = "Predicted number of faint bat passes"
)

p_moon_spkr <- plot_moon_effect(
  pred_moon_spkr,
  panel_title = "B",
  y_title = NULL
)

# -------------------------------------------------------------------------
# Combine panels
# -------------------------------------------------------------------------

p_moon_combined <-
  p_moon_rob + p_moon_spkr +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste(
      "bat calls at moth decoys (A) and speaker lures (B)",
      "by treatment and moonlight"
    ),
    subtitle = paste(
      "Lines show community-average model predictions;",
      "shaded areas show 95% confidence intervals;",
      "models include offsets for sampling effort and species-specific background activity"
    )
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10)
  )

p_moon_combined

# -------------------------------------------------------------------------
# Save figure
# -------------------------------------------------------------------------

dir.create(
  "figures/rob_spkr_model/v3",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = "figures/rob_spkr_model/v3/moonlight_marginal_effects_m10_models.png",
  plot = p_moon_combined,
  width = 9,
  height = 5,
  dpi = 300
)


ggsave(
  plot = p_moon_combined,
  filename = "figures/rob_spkr_model/v3/moon_effects_combined.tiff",
  width = 8,
  height = 4,
  units = "in",
  dpi = 600,
  compression = "lzw", # this compression just works with tiff files. 
  bg = "white"
)


# Old script  -------------------------------------------------------------
# 
# 
# 
# 

# # currently the best model seems to be m10. but we still have to add the offset. 
# 
# m12 <- glmmTMB(
#   n_calls ~ trmt_bin + jday_s + I(jday_s^2) +
#     nit_avg_wspm_s_s + avg_moonlight_s + elev_max_s +
#     yr_s + t_leps_s +
#     offset(log(effort_hours)) +
#     (1 | site) +
#     (1 + trmt_bin | sp) +
#     trmt_bin * jday_s +
#     trmt_bin * I(jday_s^2) +
#     trmt_bin * avg_moonlight_s,
#   family = nbinom2,
#   data = rob_db
# )
# 
# check_singularity(m12)
# m12$sdr$pdHess
# m12$fit$message
# check_zeroinflation(m12)
# calculate_c_hat(m12)   
# performance_mae(m12)
# performance::r2(m12)
# DHARMa::simulateResiduals(m12, n = 1000, plot = TRUE)
# performance::check_collinearity(m12)
# summary(m12)
# anova( m10, m11, m12)
# 



# # the section below is an old part of the script and might be removed in future iterations 
# 
# # 
# # model m11 --------------------------------------------------------------
# 
# 
# library(dplyr)
# library(ggplot2)
# library(marginaleffects)
# library(tidyr)
# 
# # -------------------------
# # 1. Species IDs to exclude
# # -------------------------
# drop_ids <- c("hilo", "hif", "hifrag", "lof", "mysp", "euma", "pahe")
# 
# # -------------------------
# # 2. Lookup table
# # rob_db already has sp_label
# # -------------------------
# sp_lookup <- rob_db %>%
#   select(sp, sp_label) %>%
#   distinct()
# 
# # keep only species you want
# keep_sp <- sp_lookup %>%
#   filter(!sp %in% drop_ids) %>%
#   pull(sp)
# 
# # -------------------------
# # 3. Typical covariate values
# # -------------------------
# typical <- list(
#   jday_s = 0,
#   nit_avg_wspm_s_s = 0,
#   nit_avg_temp_c_s = 0,
#   avg_moonlight_s = 0,
#   t_leps_s = 0,
#   elev_max_s = mean(rob_db$elev_max_s, na.rm = TRUE),
#   yr_s = mean(rob_db$yr_s, na.rm = TRUE)
# )
# 
# # -------------------------
# # 4. Community-level predictions
# # -------------------------
# pred_comm <- predictions(
#   m11,
#   newdata = datagrid(
#     trmt_bin = c(-1, 1),
#     jday_s = typical$jday_s,
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     t_leps_s = typical$t_leps_s,
#     elev_max_s = typical$elev_max_s,
#     yr_s = typical$yr_s
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit")),
#     panel = "Community"
#   )
# 
# # -------------------------
# # 5. Species-level predictions
# # predictions() returns sp, but not sp_label
# # so join sp_lookup AFTER prediction
# # -------------------------
# pred_sp <- predictions(
#   m11,
#   newdata = datagrid(
#     sp = sort(unique(rob_db$sp[rob_db$sp %in% keep_sp])),
#     trmt_bin = c(-1, 1),
#     jday_s = typical$jday_s,
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     t_leps_s = typical$t_leps_s,
#     elev_max_s = typical$elev_max_s,
#     yr_s = typical$yr_s
#   ),
#   type = "response",
#   re.form = NULL
# ) %>%
#   left_join(sp_lookup, by = "sp") %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit")),
#     panel = if_else(is.na(sp_label), sp, sp_label)
#   )
# 
# # -------------------------
# # 6. Order panels by treatment effect
# # -------------------------
# panel_order <- pred_sp %>%
#   select(panel, trmt, estimate) %>%
#   pivot_wider(names_from = trmt, values_from = estimate) %>%
#   mutate(light_effect = Lit - Dark) %>%
#   arrange(desc(light_effect)) %>%
#   pull(panel)
# 
# panel_levels <- c("Community", panel_order)
# 
# # -------------------------
# # 7. Combine predictions
# # -------------------------
# pred_all <- bind_rows(
#   pred_comm %>% select(panel, trmt, estimate, conf.low, conf.high),
#   pred_sp %>% select(panel, trmt, estimate, conf.low, conf.high)
# ) %>%
#   mutate(panel = factor(panel, levels = panel_levels))
# 
# # -------------------------
# # 8. Raw data for plotting
# # rob_db already has sp_label, so DO NOT join again
# # -------------------------
# raw_sp <- rob_db %>%
#   filter(sp %in% keep_sp) %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit")),
#     panel = if_else(is.na(sp_label), sp, sp_label)
#   )
# 
# raw_comm <- rob_db %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit")),
#     panel = "Community"
#   )
# 
# raw_all <- bind_rows(
#   raw_comm %>% select(panel, trmt, c_buzz),
#   raw_sp %>% select(panel, trmt, c_buzz)
# ) %>%
#   mutate(panel = factor(panel, levels = panel_levels))
# 
# # we don't need all the poits we can do a subset 
# 
# raw_all_sub <- raw_all %>% sample_frac(0.25)
# 
# # -------------------------
# # 9. Plot
# # -------------------------
# p_trmt_m11 <- ggplot() +
#   geom_jitter(
#     data = raw_all_sub,
#     aes(x = trmt, y = c_buzz),
#     width = 0.08, height = 0,
#     alpha = 0.08, size = 0.7
#   ) +
#   geom_line(
#     data = pred_all,
#     aes(x = trmt, y = estimate, group = 1),
#     linewidth = 0.5
#   ) +
#   geom_point(
#     data = pred_all,
#     aes(x = trmt, y = estimate),
#     size = 1.8
#   ) +
#   geom_errorbar(
#     data = pred_all,
#     aes(x = trmt, ymin = conf.low, ymax = conf.high),
#     width = 0.06,
#     linewidth = 0.4
#   ) +
#   facet_wrap(~ panel, scales = "free_y", ncol = 4) +
#   scale_y_continuous(trans = "log1p") +
#   labs(
#     x = "Treatment",
#     y = "Predicted feeding buzzes",
#     title = ""
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(
#     legend.position = "none",
#     axis.title = element_text(size = 12),
#     axis.text  = element_text(size = 12),
#     plot.title = element_text(size = 12),
#     strip.text = element_text(size = 11),
#     text = element_text(family = "serif")
#   )
# 
# p_trmt_m11
# 
# ggsave(
#   filename = "figures/rob_spkr_model/robomoth_light.tiff",
#   plot = p_trmt_m11,
#   width = 18,   # adjust width
#   height = 15,   # adjust height
#   dpi = 300,    # publication quality
#   units = "cm", # inches
#   device = "tiff",
#   compression = "lzw"  # smaller file size
# )
# 
# 
# # marginal effects for the jday
# 
# sdjday<- sd(rob_db$jday, na.rm = TRUE)          # raw SD of temp_c (ignoring NAs)
# meanjday<- mean(rob_db$jday, na.rm = TRUE)  
# 
# # sequence across observed scaled range
# j_seq <- seq(
#   min(rob_db$jday_s, na.rm = TRUE),
#   max(rob_db$jday_s, na.rm = TRUE),
#   length.out = 100
# )
# 
# # typical values for other predictors
# typical <- list(
#   nit_avg_wspm_s_s = 0,
#   nit_avg_temp_c_s = 0,
#   avg_moonlight_s = 0,
#   t_leps_s = 0,
#   elev_max_s = mean(rob_db$elev_max_s, na.rm = TRUE),
#   yr_s = mean(rob_db$yr_s, na.rm = TRUE)
# )
# 
# # marginal predictions
# pred_jday <- predictions(
#   m11,
#   newdata = datagrid(
#     jday_s = j_seq,
#     trmt_bin = c(-1, 1),
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     t_leps_s = typical$t_leps_s,
#     elev_max_s = typical$elev_max_s,
#     yr_s = typical$yr_s
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit")),
#     jday = jday_s * (2*sdjday) + meanjday  # back-transform to raw Julian day
#   )
# 
# # raw data for plotting behind the marginal curves
# raw_jday <- rob_db %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit"))
#   )
# # I dont't need all the raw points for the plot just a fraction 
# raw_jday_sub <- raw_jday %>% sample_frac(0.25)
# 
# # plot
# p_jday_raw <- ggplot() +
#   geom_point(
#     data = raw_jday_sub,
#     aes(x = jday, y = c_buzz ),
#     alpha = 0.20,
#     size = 0.9,
#     position = position_jitter(width = 0.3, height = 0)
#   ) +
#   geom_ribbon(
#     data = pred_jday,
#     aes(x = jday, y = estimate, ymin = conf.low, ymax = conf.high, fill = trmt),
#     alpha = 0.25,
#     color = NA
#   ) +
#   geom_line(
#     data = pred_jday,
#     aes(x = jday, y = estimate, color = trmt),
#     linewidth = 1.2
#   ) +
#   scale_color_manual(values = c("Dark" = "grey20", "Lit" = "grey70")) +
#   scale_fill_manual(values = c("Dark" = "grey20", "Lit" = "grey70"))+
#   scale_y_continuous(trans = "log1p") +
#   labs(
#     x = "Julian day",
#     y = "Predicted feeding buzzes",
#     color = "Treatment",
#     fill = "Treatment",
#     title = "Seasonal effect of artificial light on feeding buzzes"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# p_jday_raw
# 
# ggsave(
#   filename = "figures/rob_spkr_model/p_jday_raw.tiff",
#   plot = p_jday_raw,
#   width = 22,
#   height = 15,
#   dpi = 600,
#   compression = "lzw",
#   units = "cm"
# )
# 
# pred_diff <- pred_jday %>%
#   select(jday_s, trmt, estimate) %>%
#   tidyr::pivot_wider(names_from = trmt, values_from = estimate) %>%
#   mutate(diff_lit_dark = Lit - Dark)
# 
# ggplot(pred_diff, aes(x = jday_s, y = diff_lit_dark)) +
#   geom_line(linewidth = 1) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   labs(
#     x = "Season (scaled Julian day)",
#     y = "Difference in predicted feeding buzzes (Lit - Dark)",
#     title = "Seasonal change in the treatment effect of light"
#   ) +
#   theme_bw(base_size = 13)
# 
# 
# 
# 
# # effect of night average wind speed
# 
# 
# # back transform wind
# 
# w_seq <- seq(
#   min(rob_db$nit_avg_wspm_s_s, na.rm = TRUE),
#   max(rob_db$nit_avg_wspm_s_s, na.rm = TRUE),
#   length.out = 100
# )
# 
# sd_wind <- sd(rob_db$nit_avg_wspm_s, na.rm = TRUE)
# mean_wind <- mean(rob_db$nit_avg_wspm_s, na.rm = TRUE)
# 
# typical <- list(
#   jday_s = 0,
#   nit_avg_temp_c_s = 0,
#   avg_moonlight_s = 0,
#   t_leps_s = 0,
#   elev_max_s = mean(rob_db$elev_max_s, na.rm = TRUE),
#   yr_s = mean(rob_db$yr_s, na.rm = TRUE)
# )
# 
# pred_wind <- predictions(
#   m11,
#   newdata = datagrid(
#     nit_avg_wspm_s_s = w_seq,
#     trmt_bin = 1,   # or -1 (doesn't matter — no interaction)
#     jday_s = typical$jday_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     t_leps_s = typical$t_leps_s,
#     elev_max_s = typical$elev_max_s,
#     yr_s = typical$yr_s
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     wind = nit_avg_wspm_s_s * (2 * sd_wind) + mean_wind
#   )
# 
# raw_wind <- rob_db
# 
# raw_wind_sub <- rob_db %>%
#   group_by(trmt_bin) %>%
#   slice_sample(prop = 0.15) %>%
#   ungroup()
# 
# 
# p_wind <- ggplot() +
#   geom_point(
#     data = raw_wind_sub,
#     aes(x = nit_avg_wspm_s, y = c_buzz),
#     alpha = 0.05,
#     size = 0.5,
#     color = "grey40",
#     position = position_jitter(width = 0.1, height = 0)
#   ) +
#   geom_ribbon(
#     data = pred_wind,
#     aes(x = wind, ymin = conf.low, ymax = conf.high),
#     fill = "grey70",
#     alpha = 0.3
#   ) +
#   geom_line(
#     data = pred_wind,
#     aes(x = wind, y = estimate),
#     color = "black",
#     linewidth = 1.2
#   ) +
#   scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
#   labs(
#     x = "Nightly average wind speed (m/s)",
#     y = "Feeding buzzes",
#     title = ""
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# p_wind
# 
# ggsave(
#   "figures/rob_spkr_model/p_wind.tiff",
#   plot = p_wind,
#   width = 22,
#   height = 15,
#   dpi = 600,
#   compression = "lzw",
#   units = "cm"
# )
# 
# 
# # moon marginal effect 
# # 
# # for unsealing predictor 
# sd_moon <- sd(rob_db$avg_moonlight, na.rm = TRUE)
# mean_moon <- mean(rob_db$avg_moonlight, na.rm = TRUE)
# 
# moon_seq <- seq(
#   min(rob_db$avg_moonlight_s, na.rm = TRUE),
#   max(rob_db$avg_moonlight_s, na.rm = TRUE),
#   length.out = 100
# )
# 
# 
# 
# typical <- list(
#   jday_s = 0,
#   nit_avg_wspm_s_s = 0,
#   nit_avg_temp_c_s = 0,
#   t_leps_s = 0,
#   elev_max_s = mean(rob_db$elev_max_s, na.rm = TRUE),
#   yr_s = mean(rob_db$yr_s, na.rm = TRUE)
# )
# 
# pred_moon <- predictions(
#   m11,
#   newdata = datagrid(
#     avg_moonlight_s = moon_seq,
#     trmt_bin = c(-1, 1),
#     jday_s = typical$jday_s,
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     t_leps_s = typical$t_leps_s,
#     elev_max_s = typical$elev_max_s,
#     yr_s = typical$yr_s
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit")),
#     moon_raw = avg_moonlight_s * (2 * sd_moon) + mean_moon
#   )
# 
# 
# 
# raw_moon <- rob_db %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit"))
#   )
# raw_moon_sub <- rob_db %>%
#   group_by(trmt_bin) %>%
#   slice_sample(prop = 0.15) %>%
#   ungroup()
# 
# p_moon <- ggplot() +
#   geom_point(
#     data = raw_moon_sub,
#     aes(x = avg_moonlight, y = c_buzz),
#     alpha = 0.05,
#     size = 1,
#     color = "grey40"
#   ) +
#   geom_ribbon(
#     data = pred_moon,
#     aes(x = moon_raw, ymin = conf.low, ymax = conf.high, fill = trmt),
#     alpha = 0.2,
#     color = NA
#   ) +
#   geom_line(
#     data = pred_moon,
#     aes(x = moon_raw, y = estimate, color = trmt),
#     linewidth = 1.1
#   ) +
#   scale_color_manual(values = c("Dark" = "grey20", "Lit" = "grey70")) +
#   scale_fill_manual(values = c("Dark" = "grey20", "Lit" = "grey70")) +
#   scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
#   labs(
#     x = "Average moonlight",
#     y = "Feeding buzzes",
#     color = "Treatment",
#     fill = "Treatment",
#     title = "Marginal effect of moonlight on feeding buzzes"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# p_moon
# 
# ggsave(
#   "figures/rob_spkr_model/p_moon.tiff", 
#   width = 22,
#   height = 15, 
#   dpi = 600,
#   compression = "lzw",
#   units = "cm"
# )
# 
# 
# # Marginal effects plot for year
# 
# # typical values for the other predictors
# typical <- list(
#   jday_s = 0,
#   nit_avg_wspm_s_s = 0,
#   nit_avg_temp_c_s = 0,
#   avg_moonlight_s = 0,
#   t_leps_s = 0,
#   elev_max_s = mean(rob_db$elev_max_s, na.rm = TRUE)
# )
# 
# # marginal predictions for year by treatment
# p_year <- predictions(
#   m11,
#   newdata = datagrid(
#     yr_s = c(-1, 1),          # 2021 = -1, 2022 = 1
#     trmt_bin = c(-1, 1),      # Dark vs Lit
#     jday_s = typical$jday_s,
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     t_leps_s = typical$t_leps_s,
#     elev_max_s = typical$elev_max_s
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     trmt = factor(trmt_bin,
#                   levels = c(-1, 1),
#                   labels = c("Dark", "Lit")
#     ),
#     year = factor(yr_s,
#                   levels = c(-1, 1),
#                   labels = c("2021", "2022")
#     )
#   )
# 
# plot_year <- ggplot(
#   p_year,
#   aes(x = year, y = estimate, group = trmt, color = trmt)
# ) +
#   geom_errorbar(
#     aes(ymin = conf.low, ymax = conf.high),
#     width = 0.1,
#     position = pd,
#     alpha = 0.5,
#     linewidth = 0.7
#   ) +
#   geom_line(
#     position = position_dodge(width = 0.2),
#     linewidth = 1
#   ) +
#   geom_point(
#     size = 3,
#     position = pd
#   ) +
#   scale_color_manual(values = c("Dark" = "grey20", "Lit" = "grey50")) +
#   labs(
#     x = "Year",
#     y = "Predicted feeding buzzes",
#     title = "Interannual variation in feeding buzzes by treatment",
#     color = "Treatment"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     text = element_text(family = "serif", size = 14),
#     plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
#     axis.title = element_text(size = 18, face = "bold"),
#     axis.text = element_text(size = 13),
#     legend.position = "bottom"
#   )
# 
# plot_year
# 
# ggsave(
#   "figures/rob_spkr_model/p_year.tiff",
#   plot = plot_year,
#   width = 22,
#   height = 15,
#   dpi = 600,
#   compression = "lzw",
#   units = "cm"
# )
# 
# # marginal effects for elevation
# 
# sd_elev <- sd(rob_db$elev_max, na.rm = TRUE)
# mean_elev <- mean(rob_db$elev_max, na.rm = TRUE)
# 
# pred_elev <- predictions(
#   m11,
#   newdata = datagrid(
#     elev_max_s = seq(
#       min(rob_db$elev_max_s, na.rm = TRUE),
#       max(rob_db$elev_max_s, na.rm = TRUE),
#       length.out = 100
#     ),
#     trmt_bin = c(-1, 1),
#     jday_s = typical$jday_s,
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     t_leps_s = typical$t_leps_s,
#     yr_s = typical$yr_s
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit")),
#     elev_raw = elev_max_s * (2 * sd_elev) + mean_elev
#   )
# 
# 
# elevation<-ggplot(pred_elev, aes(x = elev_raw, y = estimate, color = trmt)) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = trmt), alpha = 0.2, color = NA) +
#   geom_line(linewidth = 1.2) +
#   scale_color_manual(values = c("Dark" = "grey20", "Lit" = "grey70")) +
#   scale_fill_manual(values = c("Dark" = "grey20", "Lit" = "grey70")) 
# 
# ggsave(
#   "figures/rob_spkr_model/p_elev.tiff",
#   plot = elevation,
#   width = 22,
#   height = 15,
#   dpi = 600,
#   compression = "lzw",
#   units = "cm"
# )
# 
# 
# #############################################################################################
# ##################Spkr analysis############################################################
#  
# # model for the speaker data 
# 
# glimpse(spkr_db)
# 
# # correct the treatment variable to be binary -1 to 1
# spkr_db <- spkr_db %>%
#   mutate(trmt_bin = if_else(treatmt == "dark", -1, 1))
# 
# m0_sp <- glmmTMB(
#   c_buzz ~ trmt_bin + (1 | site),
#   family = poisson(link = "log"),
#   data = spkr_db
# )
# 
# check_singularity(m0_sp)
# m0_sp$sdr$pdHess
# m0_sp$fit$message
# check_zeroinflation(m0_sp)
# calculate_c_hat(m0_sp)   
# performance_mae(m0_sp)
# range(rob_db$c_buzz)
# hist(rob_db$c_buzz, breaks = 100 )
# performance::r2(m0_sp)
# DHARMa::simulateResiduals(m0_sp, n = 1000, plot = TRUE)
# performance::check_collinearity(m0_sp)
# summary(m0_sp)
# 
# 
# 
# m1_sp <- glmmTMB(
#   c_buzz ~ trmt_bin + (1 | site),
#   family = nbinom2,
#   data = spkr_db
# )
# 
# 
# check_singularity(m1_sp)
# m1_sp$sdr$pdHess
# m1_sp$fit$message
# check_zeroinflation(m1_sp)
# calculate_c_hat(m1_sp)   
# performance_mae(m1_sp)
# range(rob_db$c_buzz)
# hist(rob_db$c_buzz, breaks = 100 )
# performance::r2(m1_sp)
# DHARMa::simulateResiduals(m1_sp, n = 1000, plot = TRUE)
# performance::check_collinearity(m1_sp)
# summary(m1_sp)
# 
# anova(m0_sp, m1_sp)
# 
# # so negative binomial is better than poisson, now we can add the seasonality
# 
# m2_sp <- glmmTMB(
#   c_buzz ~ trmt_bin +
#     jday_s + I(jday_s^2) +
#     (1 | site),
#   family = nbinom2,
#   data = spkr_db
# )
# 
# check_singularity(m1_sp)
# m1_sp$sdr$pdHess
# m1_sp$fit$message
# check_zeroinflation(m1_sp)
# calculate_c_hat(m1_sp)   
# performance_mae(m1_sp)
# range(rob_db$c_buzz)
# hist(rob_db$c_buzz, breaks = 100 )
# performance::r2(m1_sp)
# DHARMa::simulateResiduals(m1_sp, n = 1000, plot = TRUE)
# performance::check_collinearity(m1_sp)
# summary(m2_sp)
# anova(m1_sp, m2_sp)
# 
# # including seasonality improves the model, light still affects feeding buzzes at speaker lures. likelihood and AIC still indicate the model m2_sp is better. now we try covariates 
# 
# 
# m3_sp <- glmmTMB(
#   c_buzz ~ trmt_bin +
#     jday_s + I(jday_s^2) +
#     nit_avg_wspm_s_s +
#     nit_avg_temp_c_s +
#     avg_moonlight_s +
#     elev_max_s +
#     yr_s +
#     (1 | site),
#   family = nbinom2,
#   data = spkr_db
# )
# 
# m4_sp <- glmmTMB(
#   c_buzz ~ trmt_bin +
#     jday_s + I(jday_s^2) +
#     nit_avg_wspm_s_s +
#     nit_avg_temp_c_s +
#     avg_moonlight_s +
#     elev_max_s +
#     yr_s +
#     t_leps_s +
#     (1 | site),
#   family = nbinom2,
#   data = spkr_db
# )
# 
# m5_sp <- glmmTMB(
#   c_buzz ~ trmt_bin +
#     jday_s + I(jday_s^2) +
#     nit_avg_wspm_s_s +
#     nit_avg_temp_c_s +
#     avg_moonlight_s +
#     elev_max_s +
#     yr_s +
#     t_leps_s +
#     (1 | site) +
#     trmt_bin * t_leps_s,
#   family = nbinom2,
#   data = spkr_db
# )
# 
# m6_sp <- glmmTMB(
#   c_buzz ~ trmt_bin +
#     jday_s + I(jday_s^2) +
#     nit_avg_wspm_s_s +
#     nit_avg_temp_c_s +
#     avg_moonlight_s +
#     elev_max_s +
#     yr_s +
#     t_leps_s +
#     (1 | site) +
#     (1 | sp),
#   family = nbinom2,
#   data = spkr_db
# )
# 
# 
# # model with all predictors and no interactions but treatmentr as random slop s by species
# 
# m7_sp <- glmmTMB(
#   c_buzz ~ trmt_bin +
#     jday_s + I(jday_s^2) +
#     nit_avg_wspm_s_s +
#     nit_avg_temp_c_s +
#     avg_moonlight_s +
#     elev_max_s +
#     yr_s +
#     t_leps_s +
#     (1 | site) +
#     (1 + trmt_bin | sp),
#   family = nbinom2,
#   data = spkr_db
# )
# 
# check_singularity(m7_sp)
# m7_sp$sdr$pdHess
# m7_sp$fit$message
# check_zeroinflation(m7_sp)
# calculate_c_hat(m7_sp)
# performance_mae(m7_sp)
# range(spkr_db$c_buzz)
# hist(spkr_db$c_buzz, breaks = 100)
# performance::r2(m7_sp)
# DHARMa::simulateResiduals(m7_sp, n = 1000, plot = TRUE)
# performance::check_collinearity(m7_sp)
# summary(m7_sp)                            
# 
# anova(m1_sp, m2_sp, m3_sp, m4_sp, m5_sp, m6_sp, m7_sp)
# 
# 
# 
# # now we test interactions between treatment and seasonality
# 
# m8_sp <- glmmTMB(
#   c_buzz ~ trmt_bin +
#     jday_s + I(jday_s^2) +
#     nit_avg_wspm_s_s +
#     nit_avg_temp_c_s +
#     avg_moonlight_s +
#     elev_max_s +
#     yr_s +
#     t_leps_s +
#     (1 | site) +
#     (1 + trmt_bin | sp) +
#     trmt_bin * jday_s + trmt_bin * I(jday_s^2),
#   family = nbinom2,
#   data = spkr_db
# )
# 
# check_singularity(m8_sp)
# m8_sp$sdr$pdHess
# m8_sp$fit$message
# check_zeroinflation(m8_sp)
# calculate_c_hat(m8_sp)
# performance_mae(m8_sp)
# range(spkr_db$c_buzz)
# hist(spkr_db$c_buzz, breaks = 100)
# performance::r2(m8_sp)
# DHARMa::simulateResiduals(m8_sp, n = 1000, plot = TRUE)
# performance::check_collinearity(m8_sp)
# summary(m8_sp)
# anova(m1_sp, m2_sp, m3_sp, m4_sp, m5_sp, m6_sp, m7_sp, m8_sp)
# 
# 
# # adding interaction with seasonality improves the model. it does increases collinearity of temperature so we might have to remove
# 
# # now we check the interaction with moonlight
# 
# m9_sp <- glmmTMB(
#   c_buzz ~ trmt_bin +
#     jday_s + I(jday_s^2) +
#     nit_avg_wspm_s_s +
#     nit_avg_temp_c_s +
#     avg_moonlight_s +
#     elev_max_s +
#     yr_s +
#     t_leps_s +
#     (1 | site) +
#     (1 + trmt_bin | sp) +
#     trmt_bin * jday_s + trmt_bin * I(jday_s^2) + trmt_bin * avg_moonlight_s,
#   family = nbinom2,
#   data = spkr_db
# )
# 
# check_singularity(m9_sp)
# m9_sp$sdr$pdHess
# m9_sp$fit$message
# check_zeroinflation(m9_sp)
# calculate_c_hat(m9_sp)
# performance_mae(m9_sp)
# range(spkr_db$c_buzz)
# hist(spkr_db$c_buzz, breaks = 100)
# performance::r2(m9_sp)
# DHARMa::simulateResiduals(m9_sp, n = 1000, plot = TRUE)
# performance::check_collinearity(m9_sp)
# summary(m9_sp)
# anova(m1_sp, m2_sp, m3_sp, m4_sp, m5_sp, m6_sp, m7_sp, m8_sp, m9_sp)
# 
# 
# 
# 
# # marginal effects  -------------------------------------------------------
# 
# 
# # species IDs to exclude
# drop_ids <- c("hilo", "hif", "hifrag", "lof", "mysp", "euma", "pahe")
# 
# # lookup table
# sp_lookup <- spkr_db %>%
#   select(sp, sp_label) %>%
#   distinct()
# 
# # keep only clearly identified species with labels
# keep_sp <- sp_lookup %>%
#   filter(!sp %in% drop_ids, !is.na(sp_label)) %>%
#   pull(sp)
# 
# typical <- list(
#   jday_s = 0,
#   nit_avg_wspm_s_s = 0,
#   nit_avg_temp_c_s = 0,
#   avg_moonlight_s = 0,
#   t_leps_s = 0,
#   elev_max_s = mean(spkr_db$elev_max_s, na.rm = TRUE),
#   yr_s = mean(spkr_db$yr_s, na.rm = TRUE)
# )
# 
# pred_comm_spkr <- predictions(
#   m9_sp,
#   newdata = datagrid(
#     trmt_bin = c(-1, 1),
#     jday_s = typical$jday_s,
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     t_leps_s = typical$t_leps_s,
#     elev_max_s = typical$elev_max_s,
#     yr_s = typical$yr_s
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit")),
#     panel = "Community"
#   )
# 
# 
# pred_sp_spkr <- predictions(
#   m9_sp,
#   newdata = datagrid(
#     sp = sort(unique(spkr_db$sp[spkr_db$sp %in% keep_sp])),
#     trmt_bin = c(-1, 1),
#     jday_s = typical$jday_s,
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     t_leps_s = typical$t_leps_s,
#     elev_max_s = typical$elev_max_s,
#     yr_s = typical$yr_s
#   ),
#   type = "response",
#   re.form = NULL
# ) %>%
#   left_join(sp_lookup, by = "sp") %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit")),
#     panel = sp_label
#   )
# 
# panel_order <- pred_sp_spkr %>%
#   select(panel, trmt, estimate) %>%
#   pivot_wider(names_from = trmt, values_from = estimate) %>%
#   mutate(light_effect = Lit - Dark) %>%
#   arrange(desc(light_effect)) %>%
#   pull(panel)
# 
# panel_levels <- c("Community", panel_order)
# 
# 
# pred_all <- bind_rows(
#   pred_comm_spkr %>% select(panel, trmt, estimate, conf.low, conf.high),
#   pred_sp_spkr %>% select(panel, trmt, estimate, conf.low, conf.high)
# ) %>%
#   mutate(panel = factor(panel, levels = panel_levels))
# 
# 
# set.seed(123)
# 
# raw_sp_sub <- spkr_db %>%
#   filter(sp %in% keep_sp) %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit")),
#     panel = sp_label
#   ) %>%
#   group_by(panel, trmt) %>%
#   slice_sample(prop = 0.15) %>%
#   ungroup()
# 
# raw_comm_sub <- spkr_db %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit")),
#     panel = "Community"
#   ) %>%
#   group_by(trmt) %>%
#   slice_sample(prop = 0.10) %>%
#   ungroup()
# 
# raw_all <- bind_rows(
#   raw_comm_sub %>% select(panel, trmt, c_buzz),
#   raw_sp_sub %>% select(panel, trmt, c_buzz)
# ) %>%
#   mutate(panel = factor(panel, levels = panel_levels))
# 
# p_trmt_spkr <- ggplot() +
#   geom_jitter(
#     data = raw_all,
#     aes(x = trmt, y = c_buzz),
#     width = 0.08,
#     height = 0,
#     alpha = 0.08,
#     size = 0.5,
#     color = "grey45"
#   ) +
#   geom_line(
#     data = pred_all,
#     aes(x = trmt, y = estimate, group = 1, color = trmt),
#     linewidth = 0.7
#   ) +
#   geom_point(
#     data = pred_all,
#     aes(x = trmt, y = estimate, color = trmt),
#     size = 1.8
#   ) +
#   geom_errorbar(
#     data = pred_all,
#     aes(x = trmt, ymin = conf.low, ymax = conf.high, color = trmt),
#     width = 0.06,
#     linewidth = 0.4
#   ) +
#   facet_wrap(~ panel, scales = "free_y", ncol = 4) +
#   scale_color_manual(values = c("Dark" = "grey20", "Lit" = "grey70")) +
#   scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
#   labs(
#     x = "Treatment",
#     y = "Feeding buzzes",
#     color = "Treatment",
#     title = "Effect of artificial light on feeding buzzes in the acoustic lure experiment"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# p_trmt_spkr
# 
# 
# ggsave(
#   filename = "figures/rob_spkr_model/spkr_light.tiff",
#   plot = p_trmt_spkr,
#   width = 22,
#   height = 15,
#   dpi = 600,
#   compression = "lzw",
#   units = "cm"
# )
# 
# 
# # the following graph need both the speaker and the robolight data to be plotted together, so both models and predictions need to be ran to work it out.
# 
# pred_all_spkr <- bind_rows(
#   pred_comm_spkr %>% select(panel, trmt, estimate, conf.low, conf.high),
#   pred_sp_spkr %>% select(panel, trmt, estimate, conf.low, conf.high)
# ) %>%
#   mutate(
#     experiment = "Speaker"
#   )
# 
# pred_all_rob <- bind_rows(
#   pred_comm %>% select(panel, trmt, estimate, conf.low, conf.high),
#   pred_sp %>% select(panel, trmt, estimate, conf.low, conf.high)
# ) %>%
#   mutate(
#     experiment = "Robomoth"
#   )
# 
# shared_panels <- intersect(
#   unique(as.character(pred_all_spkr$panel)),
#   unique(as.character(pred_all_rob$panel))
# )
# 
# pred_compare <- bind_rows(pred_all_rob, pred_all_spkr) %>%
#   filter(panel %in% shared_panels)
# 
# 
# library(dplyr)
# library(tidyr)
# 
# panel_order <- pred_compare %>%
#   select(panel, experiment, trmt, estimate) %>%
#   pivot_wider(names_from = trmt, values_from = estimate) %>%
#   mutate(light_effect = Lit - Dark) %>%
#   group_by(panel) %>%
#   summarise(mean_effect = mean(light_effect, na.rm = TRUE), .groups = "drop") %>%
#   arrange(desc(mean_effect)) %>%
#   pull(panel)
# 
# panel_levels <- c("Community", setdiff(panel_order, "Community"))
# 
# pred_compare <- pred_compare %>%
#   mutate(
#     panel = factor(panel, levels = panel_levels),
#     experiment = factor(experiment, levels = c("Robomoth", "Speaker"))
#   )
# 
# library(ggplot2)
# 
# pd <- position_dodge(width = 0.18)
# 
# p_compare <- ggplot(
#   pred_compare,
#   aes(x = trmt, y = estimate, group = experiment, color = experiment)
# ) +
#   geom_errorbar(
#     aes(ymin = conf.low, ymax = conf.high),
#     width = 0.06,
#     linewidth = 0.45,
#     position = pd
#   ) +
#   geom_line(
#     linewidth = 0.7,
#     position = pd
#   ) +
#   geom_point(
#     size = 2,
#     position = pd
#   ) +
#   facet_wrap(~ panel, scales = "free_y", ncol = 4) +
#   scale_color_manual(values = c("Robomoth" = "grey20", "Speaker" = "grey60")) +
#   scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
#   labs(
#     x = "Treatment",
#     y = "Predicted feeding buzzes",
#     color = "Experiment",
#     title = "Treatment effects on feeding buzzes across experiments"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(
#     legend.position = "bottom",
#     axis.title = element_text(size = 12),
#     axis.text = element_text(size = 11),
#     strip.text = element_text(size = 10),
#     plot.title = element_text(hjust = 0.5),
#     text = element_text(family = "serif")
#   )
# 
# p_compare
# 
# ggsave(
#   filename = "figures/rob_spkr_model/p_compare.png",
#   plot = p_compare,
#   width = 22,
#   height = 15,
#   dpi = 600,
#   # compression = "lzw",
#   units = "cm"
# )
# 
# 
# 
# # julian day graph 
# 
# jday_seq <- seq(
#   min(spkr_db$jday_s, na.rm = TRUE),
#   max(spkr_db$jday_s, na.rm = TRUE),
#   length.out = 100
# )
# 
# 
# typical <- list(
#   nit_avg_wspm_s_s = 0,
#   nit_avg_temp_c_s = 0,
#   avg_moonlight_s = 0,
#   t_leps_s = 0,
#   elev_max_s = mean(spkr_db$elev_max_s, na.rm = TRUE),
#   yr_s = mean(spkr_db$yr_s, na.rm = TRUE)
# )
# 
# pred_jday <- predictions(
#   m9_sp,
#   newdata = datagrid(
#     jday_s = jday_seq,
#     trmt_bin = c(-1, 1),
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     t_leps_s = typical$t_leps_s,
#     elev_max_s = typical$elev_max_s,
#     yr_s = typical$yr_s
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit"))
#   )
# 
# mean_jday <- mean(spkr_db$jday, na.rm = TRUE)
# sd_jday   <- sd(spkr_db$jday, na.rm = TRUE)
# 
# pred_jday <- pred_jday %>%
#   mutate(
#     jday = jday_s * (2 * sd_jday) + mean_jday
#   )
# 
# p_jday <- ggplot(pred_jday,
#                  aes(x = jday, y = estimate, color = trmt, fill = trmt)
# ) +
#   geom_ribbon(
#     aes(ymin = conf.low, ymax = conf.high),
#     alpha = 0.2,
#     color = NA
#   ) +
#   geom_line(linewidth = 1) +
#   scale_color_manual(values = c("Dark" = "grey20", "Lit" = "grey70")) +
#   scale_fill_manual(values = c("Dark" = "grey20", "Lit" = "grey70")) +
#   labs(
#     x = "Julian day",
#     y = "Predicted feeding buzzes",
#     color = "Treatment",
#     fill = "Treatment",
#     title = "Seasonal effect of Julian day on feeding buzzes"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     text = element_text(family = "serif"),
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# p_jday
# 
# 
# set.seed(123)
# 
# raw_sub <- spkr_db %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit"))
#   ) %>%
#   sample_frac(0.05)
# 
# p_jday +
#   geom_point(
#     data = raw_sub,
#     aes(x = jday, y = c_buzz),
#     alpha = 1,
#     size = 0.5,
#     inherit.aes = FALSE
#   )
# 
# ggsave(
#   filename = "figures/rob_spkr_model/p_jday_spkr.tiff",
#   plot = p_jday,
#   width = 22,
#   height = 15,
#   dpi = 600,
#   compression = "lzw",
#   units = "cm"
# )
# 
# 
# 
# # moon_light marginal effects plot 
# 
# moon_seq <- seq(
#   min(spkr_db$avg_moonlight_s, na.rm = TRUE),
#   max(spkr_db$avg_moonlight_s, na.rm = TRUE),
#   length.out = 100
# )
# 
# typical <- list(
#   jday_s = 0,
#   `jday_s^2` = 0,   # optional but not used directly
#   nit_avg_wspm_s_s = 0,
#   nit_avg_temp_c_s = 0,
#   t_leps_s = 0,
#   elev_max_s = mean(spkr_db$elev_max_s, na.rm = TRUE),
#   yr_s = mean(spkr_db$yr_s, na.rm = TRUE)
# )
# 
# 
# # pred_moon <- predictions(
# #   m9_sp,
# #   newdata = datagrid(
# #     avg_moonlight_s = moon_seq,
# #     trmt_bin = c(-1, 1),
# #     
# #     # main variable held constant
# #     jday_s = typical$jday_s,
# #     `I(jday_s^2)` = typical$jday_s^2,   # 🔴 REQUIRED FIX
# #     
# #     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
# #     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
# #     t_leps_s = typical$t_leps_s,
# #     elev_max_s = typical$elev_max_s,
# #     yr_s = typical$yr_s,
# #     
# #     # grouping variables (safe to include)
# #     site = unique(spkr_db$site)[1],
# #     sp = unique(spkr_db$sp)[1]
# #   ),
# #   type = "response",
# #   re.form = NA
# # ) %>%
# #   mutate(
# #     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit"))
# #   )
# # 
# 
# 
# 
# pred_moon <- predictions(
#   m9_sp,
#   newdata = datagrid(
#     avg_moonlight_s = moon_seq,
#     trmt_bin = c(-1, 1),
#     jday_s = typical$jday_s,
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     t_leps_s = typical$t_leps_s,
#     elev_max_s = typical$elev_max_s,
#     yr_s = typical$yr_s,
#     site = unique(spkr_db$site)[1],
#     sp = unique(spkr_db$sp)[1]
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit"))
#   )
# 
# mean_moon <- mean(spkr_db$avg_moonlight, na.rm = TRUE)
# sd_moon   <- sd(spkr_db$avg_moonlight, na.rm = TRUE)
# 
# pred_moon <- pred_moon %>%
#   mutate(
#     avg_moonlight = avg_moonlight_s * (2 * sd_moon) + mean_moon
#   )
# 
# 
# 
# p_moon <- ggplot(
#   pred_moon,
#   aes(x = avg_moonlight, y = estimate, color = trmt, fill = trmt)
# ) +
#   geom_ribbon(
#     aes(ymin = conf.low, ymax = conf.high),
#     alpha = 0.2,
#     color = NA
#   ) +
#   geom_line(linewidth = 1) +
#   scale_color_manual(values = c("Dark" = "grey20", "Lit" = "grey70")) +
#   scale_fill_manual(values = c("Dark" = "grey20", "Lit" = "grey70")) +
#   scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
#   labs(
#     x = "Average moonlight",
#     y = "Predicted feeding buzzes",
#     color = "Treatment",
#     fill = "Treatment",
#     title = "Effect of moonlight on feeding buzzes"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     text = element_text(family = "serif"),
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# p_moon
# 
# 
# set.seed(123)
# 
# raw_moon_sub <- spkr_db %>%
#   mutate(
#     trmt = factor(trmt_bin, levels = c(-1, 1), labels = c("Dark", "Lit"))
#   ) %>%
#   sample_frac(0.05)
# 
# 
# p_moon +
#   geom_point(
#     data = raw_moon_sub,
#     aes(x = avg_moonlight, y = c_buzz),
#     inherit.aes = FALSE,
#     alpha = 0.05,
#     size = 0.5,
#     color = "grey40"
#   )
# 
# ggsave(
#   filename = "figures/rob_spkr_model/p_moon_spkr.tiff",
#   plot = p_moon,
#   width = 22,
#   height = 15,
#   dpi = 600,
#   # compression = "lzw",
#   units = "cm"
# )
# 
# 
# # now we do the marginal effects graph for wind speed. 
# 
# 
# # sequence across the observed scaled wind range
# wind_seq <- seq(
#   min(spkr_db$nit_avg_wspm_s_s, na.rm = TRUE),
#   max(spkr_db$nit_avg_wspm_s_s, na.rm = TRUE),
#   length.out = 100
# )
# 
# # raw wind column is nit_avg_wspm_s (despite the confusing name)
# mean_wind <- mean(spkr_db$nit_avg_wspm_s, na.rm = TRUE)
# sd_wind   <- sd(spkr_db$nit_avg_wspm_s, na.rm = TRUE)
# 
# # typical values for the other predictors
# typical <- list(
#   trmt_bin = 0,   # midpoint between Dark (-1) and Lit (1)
#   jday_s = 0,
#   nit_avg_temp_c_s = 0,
#   avg_moonlight_s = 0,
#   t_leps_s = 0,
#   elev_max_s = mean(spkr_db$elev_max_s, na.rm = TRUE),
#   yr_s = mean(spkr_db$yr_s, na.rm = TRUE)
# )
# 
# pred_wind <- predictions(
#   m9_sp,
#   newdata = datagrid(
#     nit_avg_wspm_s_s = wind_seq,
#     trmt_bin = typical$trmt_bin,
#     jday_s = typical$jday_s,
#     `I(jday_s^2)` = typical$jday_s^2,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     t_leps_s = typical$t_leps_s,
#     elev_max_s = typical$elev_max_s,
#     yr_s = typical$yr_s,
#     site = unique(spkr_db$site)[1],
#     sp = unique(spkr_db$sp)[1]
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     wind = nit_avg_wspm_s_s * (2 * sd_wind) + mean_wind
#   )
# 
# p_wind <- ggplot(pred_wind, aes(x = wind, y = estimate)) +
#   geom_ribbon(
#     aes(ymin = conf.low, ymax = conf.high),
#     fill = "grey70",
#     alpha = 0.25
#   ) +
#   geom_line(
#     color = "grey20",
#     linewidth = 1
#   ) +
#   # scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
#   labs(
#     x = "Nightly average wind speed (m/s)",
#     y = "Predicted feeding buzzes",
#     title = "Marginal effect of wind speed on feeding buzzes"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     text = element_text(family = "serif"),
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# p_wind
# 
# ggsave(
#   filename = "figures/rob_spkr_model/p_wind_spkdr.tiff",
#   plot = p_wind,
#   width = 22,
#   height = 15,
#   dpi = 600,
#   compression = "lzw",
#   units = "cm"
# )
# 
# # now the temperature marginal effects plot
# 
# 
# 
# # sequence across the observed scaled temperature range
# temp_seq <- seq(
#   min(spkr_db$nit_avg_temp_c_s, na.rm = TRUE),
#   max(spkr_db$nit_avg_temp_c_s, na.rm = TRUE),
#   length.out = 100
# )
# 
# # raw temperature mean and SD
# mean_temp <- mean(spkr_db$nit_avg_temp_c, na.rm = TRUE)
# sd_temp   <- sd(spkr_db$nit_avg_temp_c, na.rm = TRUE)
# 
# # typical values for the other predictors
# typical <- list(
#   trmt_bin = 0,   # midpoint between Dark (-1) and Lit (1)
#   jday_s = 0,
#   nit_avg_wspm_s_s = 0,
#   avg_moonlight_s = 0,
#   t_leps_s = 0,
#   elev_max_s = mean(spkr_db$elev_max_s, na.rm = TRUE),
#   yr_s = mean(spkr_db$yr_s, na.rm = TRUE)
# )
# 
# pred_temp <- predictions(
#   m9_sp,
#   newdata = datagrid(
#     nit_avg_temp_c_s = temp_seq,
#     trmt_bin = typical$trmt_bin,
#     jday_s = typical$jday_s,
#     `I(jday_s^2)` = typical$jday_s^2,
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     t_leps_s = typical$t_leps_s,
#     elev_max_s = typical$elev_max_s,
#     yr_s = typical$yr_s,
#     site = unique(spkr_db$site)[1],
#     sp = unique(spkr_db$sp)[1]
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     temp_c = nit_avg_temp_c_s * (2 * sd_temp) + mean_temp
#   )
# 
# 
# p_temp <- ggplot(pred_temp, aes(x = temp_c, y = estimate)) +
#   geom_ribbon(
#     aes(ymin = conf.low, ymax = conf.high),
#     fill = "grey70",
#     alpha = 0.25
#   ) +
#   geom_line(
#     color = "grey20",
#     linewidth = 1
#   ) +
#   # scale_y_continuous() +
#   labs(
#     x = "Nightly average temperature (°C)",
#     y = "Predicted feeding buzzes",
#     title = "Marginal effect of temperature on feeding buzzes"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     text = element_text(family = "serif"),
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# p_temp
# 
# ggsave(
#   filename = "figures/rob_spkr_model/p_temp_spkr.tiff",
#   plot = p_temp,
#   width = 22,
#   height = 15,
#   dpi = 600,
#   compression = "lzw",
#   units = "cm"
# )
# 
# # now the elevation marginal effects plot for the speaker data.
# 
# # elevation spkr ----------------------------------------------------------
# 
# 
# # sequence across the observed scaled elevation range
# elev_seq <- seq(
#   min(spkr_db$elev_max_s, na.rm = TRUE),
#   max(spkr_db$elev_max_s, na.rm = TRUE),
#   length.out = 100
# )
# 
# # raw elevation mean and SD
# mean_elev <- mean(spkr_db$elev_max, na.rm = TRUE)
# sd_elev   <- sd(spkr_db$elev_max, na.rm = TRUE)
# 
# # typical values for the other predictors
# typical <- list(
#   trmt_bin = 0,
#   jday_s = 0,
#   nit_avg_wspm_s_s = 0,
#   nit_avg_temp_c_s = 0,
#   avg_moonlight_s = 0,
#   t_leps_s = 0,
#   yr_s = mean(spkr_db$yr_s, na.rm = TRUE)
# )
# 
# pred_elev <- predictions(
#   m9_sp,
#   newdata = datagrid(
#     elev_max_s = elev_seq,
#     trmt_bin = typical$trmt_bin,
#     jday_s = typical$jday_s,
#     `I(jday_s^2)` = typical$jday_s^2,
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     t_leps_s = typical$t_leps_s,
#     yr_s = typical$yr_s,
#     site = unique(spkr_db$site)[1],
#     sp = unique(spkr_db$sp)[1]
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     elev_max = elev_max_s * (2 * sd_elev) + mean_elev
#   )
# 
# p_elev <- ggplot(pred_elev, aes(x = elev_max, y = estimate)) +
#   geom_ribbon(
#     aes(ymin = conf.low, ymax = conf.high),
#     fill = "grey70",
#     alpha = 0.25
#   ) +
#   geom_line(
#     color = "grey20",
#     linewidth = 1
#   ) +
#   # scale_y_continuous(trans = log1p()) + # this is not working
#   labs(
#     x = "Maximum elevation (m)",
#     y = "Predicted feeding buzzes",
#     title = "Marginal effect of elevation on feeding buzzes"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     text = element_text(family = "serif"),
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# p_elev
# 
# ggsave(
#   filename = "figures/rob_spkr_model/p_elev_spkr.tiff",
#   plot = p_elev,
#   width = 22,
#   height = 15,
#   dpi = 600,
#   compression = "lzw",
#   units = "cm"
# )
# 
# 
# 
# # year marginal spkr ------------------------------------------------------
# 
# 
# # typical values for the other predictors
# typical <- list(
#   trmt_bin = 0,   # midpoint between Dark (-1) and Lit (1)
#   jday_s = 0,
#   nit_avg_wspm_s_s = 0,
#   nit_avg_temp_c_s = 0,
#   avg_moonlight_s = 0,
#   t_leps_s = 0,
#   elev_max_s = mean(spkr_db$elev_max_s, na.rm = TRUE)
# )
# 
# pred_year <- predictions(
#   m9_sp,
#   newdata = datagrid(
#     yr_s = c(-1, 1),
#     trmt_bin = typical$trmt_bin,
#     jday_s = typical$jday_s,
#     `I(jday_s^2)` = typical$jday_s^2,
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     t_leps_s = typical$t_leps_s,
#     elev_max_s = typical$elev_max_s,
#     site = unique(spkr_db$site)[1],
#     sp = unique(spkr_db$sp)[1]
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     year = factor(yr_s, levels = c(-1, 1), labels = c("2022", "2023"))
#   )
# 
# p_year <- ggplot(pred_year, aes(x = year, y = estimate, group = 1)) +
#   geom_errorbar(
#     aes(ymin = conf.low, ymax = conf.high),
#     width = 0.08,
#     linewidth = 0.5,
#     color = "grey30"
#   ) +
#   geom_line(
#     linewidth = 0.8,
#     color = "grey30"
#   ) +
#   geom_point(
#     size = 2.5,
#     color = "grey20"
#   ) +
#   scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
#   labs(
#     x = "Year",
#     y = "Predicted feeding buzzes",
#     title = "Marginal effect of year on feeding buzzes"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     text = element_text(family = "serif"),
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# p_year
# 
# ggsave(
#   filename = "figures/rob_spkr_model/p_year_spkr.tiff",
#   plot = p_year,
#   width = 22,
#   height = 15,
#   dpi = 600,
#   compression = "lzw",
#   units = "cm"
# )
# 
# 
# # marginal effects leps spkr ----------------------------------------------
# 
# 
# # marginal effects for the number of leps in the area.
# 
# # sequence across the observed scaled Lepidoptera range
# leps_seq <- seq(
#   min(spkr_db$t_leps_s, na.rm = TRUE),
#   max(spkr_db$t_leps_s, na.rm = TRUE),
#   length.out = 100
# )
# 
# # raw Lepidoptera mean and SD
# mean_leps <- mean(spkr_db$t_leps, na.rm = TRUE)
# sd_leps   <- sd(spkr_db$t_leps, na.rm = TRUE)
# 
# # typical values for the other predictors
# typical <- list(
#   trmt_bin = 0,   # midpoint between Dark (-1) and Lit (1)
#   jday_s = 0,
#   nit_avg_wspm_s_s = 0,
#   nit_avg_temp_c_s = 0,
#   avg_moonlight_s = 0,
#   elev_max_s = mean(spkr_db$elev_max_s, na.rm = TRUE),
#   yr_s = mean(spkr_db$yr_s, na.rm = TRUE)
# )
# 
# pred_leps <- predictions(
#   m9_sp,
#   newdata = datagrid(
#     t_leps_s = leps_seq,
#     trmt_bin = typical$trmt_bin,
#     jday_s = typical$jday_s,
#     `I(jday_s^2)` = typical$jday_s^2,
#     nit_avg_wspm_s_s = typical$nit_avg_wspm_s_s,
#     nit_avg_temp_c_s = typical$nit_avg_temp_c_s,
#     avg_moonlight_s = typical$avg_moonlight_s,
#     elev_max_s = typical$elev_max_s,
#     yr_s = typical$yr_s,
#     site = unique(spkr_db$site)[1],
#     sp = unique(spkr_db$sp)[1]
#   ),
#   type = "response",
#   re.form = NA
# ) %>%
#   mutate(
#     t_leps = t_leps_s * (2 * sd_leps) + mean_leps
#   )
# 
# p_leps <- ggplot(pred_leps, aes(x = t_leps, y = estimate)) +
#   geom_ribbon(
#     aes(ymin = conf.low, ymax = conf.high),
#     fill = "grey70",
#     alpha = 0.25
#   ) +
#   geom_line(
#     color = "grey20",
#     linewidth = 1
#   ) +
#   scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
#   labs(
#     x = "Lepidoptera abundance",
#     y = "Predicted feeding buzzes",
#     title = "Marginal effect of Lepidoptera abundance on feeding buzzes"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     text = element_text(family = "serif"),
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# p_leps
# 
# ggsave(
#   filename = "figures/rob_spkr_model/p_leps_spkr.tiff",
#   plot = p_leps,
#   width = 22,
#   height = 15,
#   dpi = 600,
#   compression = "lzw",
#   units = "cm"
# )
# 
# 
