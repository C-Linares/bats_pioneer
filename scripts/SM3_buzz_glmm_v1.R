#=============================================================================
# Script Name:    SM3_buzz_glmm_v1.R
# Purpose:        Fit GLMMs for bat paper using weather, and moon data, and light spectra.
#
# Version:        v1 - first go at the models

# #
# Author:         Carlos Linares
# Collaborator:   
# Created:        2025-02-11
# Contact:        carlosgarcialina@u.boisestate.edu
#
# Notes:
#   - Data prepped using prep_for_glmm_v2.R
#   - Outputs include marginal effect and random effects plots stored in the figures folder. 
#
# Inputs:
#   - data_for_analysis/prep_for_glmm_v2/sm3_buzz_v2.csv         : Bat call counts, predictors, standardized 
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
  "data.table", "janitor", "patchwork", "flextable", 
  "gtsummary", "lubridate", "beepr", "emmeans", "ggpubr", "ggrepel"
)

# Load all packages in one go
invisible(lapply(packages, library, character.only = TRUE))


# funtions ----------------------------------------------------------------


# Define a function to calculate and print c-hat
calculate_c_hat <- function(model) {
  residual_deviance <- deviance(model)
  residual_df <- df.residual(model)
  c_hat_deviance <- residual_deviance / residual_df
  print(c_hat_deviance)
}


# load data  --------------------------------------------------------------
sm3_buzz_v2 <- read.csv("data_for_analysis/prep_for_glmm_v2/sm3_buzz_v2.csv")


# confirm that trmt and treatment agree
sm3_buzz_v2 %>%
  count(treatmt, trmt_bin)


ggplot(sm3_buzz_v2, aes(t_buzzes)) +
  geom_histogram(binwidth = 1) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(
    x = "Feeding buzzes per species-night",
    y = "Number of observations"
  ) +
  theme_classic()

sm3_buzz_v2 %>%
  arrange(desc(t_buzzes)) %>%  # we have two nights with feeding buzzes above 2000 we might have to double check. 
  select(
    site, noche, sp_clean, treatmt,
    t_buzzes, identified_files, eff.hrs
  ) %>%
  slice_head(n = 30)

# modeling ----------------------------------------------------------------


# Note: the model we are using as the global model is m1.8 for this data





# model with treatment as binomial 

# first try poisson and just the data summary. 

m0<-glmmTMB(t_buzzes ~ trmt_bin + (1 | site),
            data = sm3_buzz_v2, family = poisson(link = "log"))

# check model
check_singularity(m0)
check_zeroinflation(m0)
calculate_c_hat(m0)  # indicates the is overdispersion c_hat= 43.42771
# performance_mae(m0)
# range(bm2$n)
summary(m0)



# model with a poisson distribution is overdispersed c-hat=43.42771 so we try negative binomial.It also shows zero inflation but this should be an effect of the over-dispersion. 

m1<- glmmTMB(
  #fixed effects
  t_buzzes ~ trmt_bin + (1 | site),
  data = sm3_buzz_v2,
  family = nbinom2(link = "log"))

# m1.1 <- glmmTMB(
#   #fixed effects
#   t_buzzes ~ trmt_bin + (1 | site),
#   data = sm3_buzz_v2,
#   family = nbinom1(link = "log"))
# 


# check model
check_singularity(m1)
check_zeroinflation(m1)
calculate_c_hat(m1)  # now c_hat= 0.62 improved fromthe last model. 
performance::r2(m1)
summary(m1)

AIC(m1, m1.1)
anova(m0, m1, m1.1) # the negattive binomial nb2 model is better thatn the nb1 and the poisson so we will chose that one. 

# model with a negative binomial has a better c_hat and AIC than the poisson model.
# now we try adding seasonality as fixed effects. R conditional is around 0.14 and marginal 0.095

m1.1<- glmmTMB(
  #fixed effects
  t_buzzes ~ trmt_bin +
    #random effects
    (1 | site) + (1  | sp_clean),  
  data = sm3_buzz_v2,
  family = nbinom2(link = "log")
)

check_singularity(m1.1)
check_zeroinflation(m1.1)
calculate_c_hat(m1.1)  # 
performance::r2(m1.1)
DHARMa::simulateResiduals(m1.1, n = 1000, plot = TRUE)
summary(m1.1)

anova(m0, m1, m1.1) # the model 
# now we add seasonality to the model 

m1.2 <- glmmTMB(
  #fixed effects
  t_buzzes ~ trmt_bin + jday_s + I(jday_s^2) +
    #random effects
    (1 | site) + (1  | sp_clean),  
  #interactions
  data = sm3_buzz_v2,
  family = nbinom2(link = "log")
)


check_singularity(m1.2)
m1.2$sdr$pdHess
m1.2$fit$message
check_zeroinflation(m1.2)
calculate_c_hat(m1.2)  # 
performance::r2(m1.2)
DHARMa::simulateResiduals(m1.2, n = 1000, plot = TRUE)
anova(m1,m1.1, m1.2)  
summary(m1.2)
res_m1.2 <- simulateResiduals(m1.2, n = 1000)

testUniformity(res_m1.2)
testDispersion(res_m1.2)
testZeroInflation(res_m1.2)
testOutliers(res_m1.2, type = "bootstrap")
# it seems like the model with seasonality is better thatn the model withouth according to AIC m1.1 142153 vs m1.2  141542.
# chisq is also smaller for model m1.2. 
# now we try the environmental predictors.



m1.3<-glmmTMB(
  #fixed effects
  t_buzzes ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1  | sp_clean),
  data = sm3_buzz_v2,
  family = nbinom2(link = "log")
)

check_singularity(m1.3)
m1.3$sdr$pdHess
m1.3$fit$message
check_zeroinflation(m1.3)
calculate_c_hat(m1.3)  
# performance_mae(m1.3)
# range(bm2$n)
# hist(bm2$n, breaks = 500)
performance::r2(m1.3)
res_m1.3<-DHARMa::simulateResiduals(m1.3, n = 1000, plot = TRUE)
testUniformity(res_m1.3)
testDispersion(res_m1.3)
testZeroInflation(res_m1.3)
testOutliers(res_m1.3, type = "bootstrap")
AIC(m1.2, m1.3) # 
anova(m1.2, m1.3) 
summary(m1.3)


# now lets try the model treatment inside random slope 

m1.4<- glmmTMB(
  #fixed effects
  t_buzzes ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin | sp_clean),  
  data = sm3_buzz_v2,
  family = nbinom2(link = "log")
)

check_singularity(m1.4)
m1.4$sdr$pdHess
m1.4$fit$message
check_zeroinflation(m1.4)
calculate_c_hat(m1.4) 
# performance_mae(m1.4)
# range(bm2$n)
# hist(bm2$n, breaks = 500)
performance::r2(m1.4)
res_m1.4<-DHARMa::simulateResiduals(m1.4, n = 1000, plot = TRUE)
testUniformity(res_m1.4)
testDispersion(res_m1.4)
testZeroInflation(res_m1.4)
testOutliers(res_m1.4, type = "bootstrap")
anova( m1.2, m1.3, m1.4)
summary(m1.4)
beep(sound = 4 )

# for model m1.4 aic shows that including treatment inside the random slopes improves the model fit. 

#now let's try seasonality inside the random slopes

m1.5<-glmmTMB(
  #fixed effects
  t_buzzes ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp_clean),  
  data = sm3_buzz_v2,
  family = nbinom2(link = "log")
)

check_singularity(m1.5)
m1.5$sdr$pdHess
m1.5$fit$message
check_zeroinflation(m1.5)
calculate_c_hat(m1.5) 
# performance_mae(m1.5)
# range(bm2$n)
# hist(bm2$n, breaks = 500)
performance::r2(m1.5)
res_m1.5<-DHARMa::simulateResiduals(m1.5, n = 1000, plot = TRUE)
testUniformity(res_m1.5)
testDispersion(res_m1.5)
testZeroInflation(res_m1.5)
testOutliers(res_m1.5, type = "bootstrap")
anova( m1.2, m1.3,m1.4,m1.5)
summary(m1.5)
beep(sound = 4 )


# allowing the experimental-light effect to vary among species substantialy improved model fit. Although the average fixed effect was not significant.
# we mild overdispersion and, a small deficit of zeros, modest uniformity, and some ouliers. 

# now let's try for interactions of treatment with seasonalty 
# we ran this model with both the square of jday interaction and the linear but we see the square is not significant so we will keep the linear interaction.

m1.6<-glmmTMB(
  #fixed effects
  t_buzzes ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp_clean) +
    #interactions
    jday_s * trmt_bin ,
  data = sm3_buzz_v2,
  family = nbinom2(link = "log")
)


check_singularity(m1.6)
m1.6$sdr$pdHess
m1.6$fit$message
check_zeroinflation(m1.6)
calculate_c_hat(m1.6)
# performance_mae(m1.6)
# range(bm2$n)
# hist(bm2$n, breaks = 500)
performance::r2(m1.6)
res_m1.6<-DHARMa::simulateResiduals(m1.6, n = 1000, plot = TRUE)
testUniformity(res_m1.6)
testDispersion(res_m1.6)
testZeroInflation(res_m1.6)
testOutliers(res_m1.6, type = "bootstrap")
performance::check_collinearity(m1.6)
anova( m1.2, m1.3,m1.4,m1.5, m1.6)
summary(m1.6)
beep(sound=4)


# adding the treatment and season  interaction we improve the fit slightly

# now we try to add the moon and treatment interaction but we will remove the I(jday_s^2) * trmt_bin interaction because it is not significant and it is not improving the model fit. we will rerun the m1.6 so it also has the trmt_bin*avg_moonlight_s interaction.

m1.7<-glmmTMB(
  #fixed effects
  t_buzzes ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp_clean) +
    #interactions
    jday_s * trmt_bin + trmt_bin*avg_moonlight_s,  
  data = sm3_buzz_v2,
  family = nbinom2(link = "log")
)


check_singularity(m1.7)
m1.7$sdr$pdHess
m1.7$fit$message
check_zeroinflation(m1.7)
calculate_c_hat(m1.7) 
# performance_mae(m1.7)
# range(bm2$n)
# hist(bm2$n, breaks = 500)
performance::r2(m1.7)
res_m1.7<-DHARMa::simulateResiduals(m1.7, n = 1000, plot = TRUE)
testUniformity(res_m1.7)
testDispersion(res_m1.7)
testZeroInflation(res_m1.7)
testOutliers(res_m1.7, type = "bootstrap")
performance::check_collinearity(m1.7)
anova( m1.2, m1.3,m1.4,m1.5,m1.6,m1.7)
summary(m1.7)
beep(sound = 4)

# gtsummary::tbl_regression(c(m1.7), exponentiate = TRUE, 
#                           label = list(
#                             trmt_bin = "Treatment (Lit vs Dark)",
#                             jday_s = "Julian Day (scaled)",
#                             # I(jday_s^2)= "Julian Day Squared (scaled)",
#                             avg_moonlight_s = "Average Moonlight (scaled)",
#                             nit_avg_temp_c_s = "Average Temperature (scaled)",
#                             nit_avg_wspm_s_s = "Average Wind Speed (scaled)",
#                             yr_s = "Year (scaled)",
#                             t_lepidoptera_s = "Lepidoptera Abundance (scaled)"
#                           )) 

#now we try the interaction with year. Although there's no effect of light as fixed effect, species specific responses to light are still strong. 
#the light response becomes more poisitive as the season progresses and less positive as light increases. insects also have a positive effect on buzz activity. 

# now we try the interaction with year. 

m1.8<-glmmTMB(
  #fixed effects
  t_buzzes ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp_clean) +
    #interactions
    jday_s * trmt_bin + trmt_bin*avg_moonlight_s + trmt_bin*yr_s,  
  data = sm3_buzz_v2,
  family = nbinom2(link = "log")
)

check_singularity(m1.8)
m1.8$sdr$pdHess
m1.8$fit$message
check_zeroinflation(m1.8)
calculate_c_hat(m1.8) 
# performance_mae(m1.8)
# range(bm2$n)
# hist(bm2$n, breaks = 500)
performance::r2(m1.8)
res_m1.8<-DHARMa::simulateResiduals(m1.8, n = 1000, plot = TRUE)
testUniformity(res_m1.8)
testDispersion(res_m1.8)
testZeroInflation(res_m1.8)
testOutliers(res_m1.8, type = "bootstrap")
performance::check_collinearity(m1.8)
anova( m1.2, m1.3,m1.4,m1.5,m1.6,m1.7, m1.8)
summary(m1.8)
beep(sound = 4)

t1<-gtsummary::tbl_regression(m1.8, exponentiate = TRUE, 
                          label = list(
                            trmt_bin = "Treatment",
                            jday_s = "Julian Day",
                             `I(jday_s^2)`= "Julian Day Squared",
                            avg_moonlight_s = "Moonlight",
                            nit_avg_temp_c_s = "Temperature",
                            nit_avg_wspm_s_s = "Wind Speed",
                            yr_s = "Year",
                            t_lepidoptera_s = "Lepidoptera Abundance"
                          )) 
as_flex_table(t1) %>%
  autofit()

# this model seems strong and it seems might be a good candidate for the final model. after runnin m1.9 it seems like this is the final model. 


# lets try interaction with lepidptera

m1.9<-glmmTMB(
  #fixed effects
  t_buzzes ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp_clean) +
    #interactions
    jday_s * trmt_bin + trmt_bin*avg_moonlight_s  + trmt_bin*t_lepidoptera_s,  
  data = sm3_buzz_v2,
  family = nbinom2(link = "log")
)

check_singularity(m1.9)
m1.9$sdr$pdHess
m1.9$fit$message
check_zeroinflation(m1.9)
calculate_c_hat(m1.9) 
# performance_mae(m1.9)
# range(bm2$n)
# hist(bm2$n, breaks = 500)
performance::r2(m1.9)
res_m1.9<-DHARMa::simulateResiduals(m1.9, n = 1000, plot = TRUE)
testUniformity(res_m1.9)
testDispersion(res_m1.9)
testZeroInflation(res_m1.9)
testOutliers(res_m1.9, type = "bootstrap")
performance::check_collinearity(m1.9)
anova( m1.2, m1.3,m1.4,m1.5,m1.6, m1.8, m1.9)
summary(m1.9)

# including the lepidoptera does not improve the model so we keep m1.7



# marginal effects plots --------------------------------------------------

# Change sp_label below if your scientific-name column has another name
species_key <- sm3_buzz_v2 %>%
  distinct(sp_clean, sp_label) %>%
  mutate(sp_clean = as.character(sp_clean))


# 2. Extract species random effects --------------------------------------

re_sp <- ranef(m1.8, condVar = TRUE) %>%
  as.data.frame() %>%
  filter(
    component == "cond",
    grpvar == "sp_clean"
  )

random_intercepts <- re_sp %>%
  filter(term == "(Intercept)") %>%
  transmute(
    sp_clean = as.character(grp),
    random_intercept = condval
  )

random_slopes <- re_sp %>%
  filter(term == "trmt_bin") %>%
  transmute(
    sp_clean = as.character(grp),
    random_slope = condval,
    random_slope_se = condsd
  )

# 3. Fixed treatment effect ----------------------------------------------

fixed_effects <- fixef(m1.8)$cond

fixed_intercept <- unname(fixed_effects["(Intercept)"])
fixed_treatment <- unname(fixed_effects["trmt_bin"])

fixed_treatment_se <- sqrt(
  vcov(m1.8)$cond["trmt_bin", "trmt_bin"]
)


# 4. Predictions and uncertainty table ----------------------------------

species_effects <- random_intercepts %>%
  left_join(random_slopes, by = "sp_clean") %>%
  left_join(species_key, by = "sp_clean") %>%
  mutate(
    # Species-specific conditional coefficients
    intercept_species = fixed_intercept + random_intercept,
    treatment_beta = fixed_treatment + random_slope,
    
    # Approximate uncertainty of species treatment coefficient
    treatment_beta_se = sqrt(
      fixed_treatment_se^2 + random_slope_se^2
    ),
    
    # trmt_bin changes from -1 to +1: contrast equals 2 × beta
    log_response_ratio = 2 * treatment_beta,
    log_ratio_se = 2 * treatment_beta_se,
    
    log_ratio_low = log_response_ratio - 1.96 * log_ratio_se,
    log_ratio_high = log_response_ratio + 1.96 * log_ratio_se,
    
    # Convert from log scale
    response_ratio = exp(log_response_ratio),
    ratio_low = exp(log_ratio_low),
    ratio_high = exp(log_ratio_high),
    
    percent_change = 100 * (response_ratio - 1),
    percent_low = 100 * (ratio_low - 1),
    percent_high = 100 * (ratio_high - 1),
    
    # Predictions at average standardized covariate values
    predicted_dark = exp(intercept_species - treatment_beta),
    predicted_lit = exp(intercept_species + treatment_beta),
    
    evidence = case_when(
      percent_low > 0 ~ "Supported increase",
      percent_high < 0 ~ "Supported decrease",
      TRUE ~ "Uncertain"
    )
  ) %>%
  arrange(percent_change)

uncertainty_table <- species_effects %>%
  transmute(
    species_code = sp_clean,
    species = sp_label,
    treatment_beta,
    treatment_beta_se,
    predicted_dark,
    predicted_lit,
    response_ratio,
    ratio_low,
    ratio_high,
    percent_change,
    percent_low,
    percent_high,
    evidence
  )

uncertainty_table


uncertainty_table_display <- uncertainty_table %>%
  mutate(
    across(
      c(
        treatment_beta,
        treatment_beta_se,
        predicted_dark,
        predicted_lit,
        response_ratio,
        ratio_low,
        ratio_high,
        percent_change,
        percent_low,
        percent_high
      ),
      ~ round(.x, 2)
    )
  )

uncertainty_table_display

flextable(uncertainty_table_display)

# uncertainty table m1.8 for pub ------------------------------------------

# 1. Format the data frame first
table_data_formatted <- uncertainty_table %>%
  # Combine columns into interval strings matching your target layout
  mutate(
    `Lit:dark ratio (95% CI)` = sprintf("%.2f (%.2f, %.2f)", 
                                        response_ratio, ratio_low, ratio_high),
    `Percent change (95% CI)` = sprintf("%+.1f%% (%+.1f%%, %+.1f%%)", 
                                        percent_change, percent_low, percent_high),
    `Absolute change`         = predicted_lit - predicted_dark
  ) %>%
  # Select and rename columns exactly as shown in the target table
  select(
    Species = species,
    `Predicted calls: dark` = predicted_dark,
    `Predicted calls: lit`  = predicted_lit,
    `Absolute change`,
    `Lit:dark ratio (95% CI)`,
    `Percent change (95% CI)`,
    Evidence = evidence
  )

# 2. Build and style the flextable
ft <- flextable(table_data_formatted) %>%
  
  # --- Number Formatting ---
  # Round numerical columns to 2 decimal places (adjust digits if needed)
  colformat_double(j = 2:4, digits = 2) %>% 
  
  # --- Text Formatting ---
  # Italicize species names
  italic(j = "Species", part = "body") %>% 
  
  # --- Alignment ---
  align(j = 1, align = "left", part = "all") %>%
  align(j = 2:4, align = "right", part = "all") %>%
  align(j = 5:7, align = "left", part = "all") %>%
  
  # --- Manuscript Styling (APA / Booktabs Style) ---
  theme_booktabs() %>%                    # Horizontal rules at top, bottom, and under header
  # font(fontname = "Arial", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  padding(padding.top = 4, padding.bottom = 4, part = "all") %>%
  autofit()

# View table in RStudio Viewer
ft


save_as_docx(
  "table.s2" = ft,
  path = "figures/sm3_glmm_v1/table/ts3.docx")

save_as_image(ft, path = "figures/sm3_glmm_v1/table/ts3.png")


######

species_order <- species_effects$sp_clean

pred_trmt <- species_effects %>%
  select(
    sp_clean,
    sp_label,
    predicted_dark,
    predicted_lit,
    percent_change,
    evidence
  ) %>%
  pivot_longer(
    cols = c(predicted_dark, predicted_lit),
    names_to = "treatment",
    values_to = "estimate"
  ) %>%
  mutate(
    trmt_bin = if_else(treatment == "predicted_dark", -1, 1),
    sp_clean = factor(sp_clean, levels = species_order)
  )

raw_plot <- sm3_buzz_v2 %>%
  mutate(
    sp_clean = factor(sp_clean, levels = species_order)
  )

facet_labels <- species_effects %>%
  mutate(
    facet_text = sprintf(
      "%s (%+.1f%%)",
      sp_label,
      percent_change
    )
  ) %>%
  select(sp_clean, facet_text) %>%
  deframe()


p_trmt <- ggplot() +
  
  # Raw feeding-buzz observations
  geom_jitter(
    data = raw_plot,
    aes(
      x = trmt_bin,
      y = t_buzzes
    ),
    width = 0.15,
    alpha = 0.15,
    size = 0.6,
    color = "grey60"
  ) +
  
  # Model predictions
  geom_line(
    data = pred_trmt,
    aes(
      x = trmt_bin,
      y = estimate,
      group = sp_clean,
      color = evidence
    ),
    linewidth = 1
  ) +
  
  geom_point(
    data = pred_trmt,
    aes(
      x = trmt_bin,
      y = estimate,
      color = evidence
    ),
    size = 2.4
  ) +
  
  facet_wrap(
    ~ sp_clean,
    scales = "free_y",
    labeller = as_labeller(facet_labels)
  ) +
  
  scale_color_manual(
    values = c(
      "Supported increase" = "black",
      "Supported decrease" = "grey35",
      "Uncertain" = "grey75"
    ),
    drop = FALSE
  ) +
  
  scale_x_continuous(
    breaks = c(-1, 1),
    labels = c("Dark", "Lit")
  ) +
  
  scale_y_continuous(
    trans = "log1p"
  ) +
  
  labs(
    x = "Treatment",
    y = "Feeding buzzes per night",
    color = "Estimated response"
  ) +
  
  theme_minimal() +
  
  theme(
    strip.text = element_text(
      size = 13,
      face = "italic"
    ),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    legend.position = "bottom"
  )

p_trmt



ggsave(
  filename = "figures/sm3_glmm_v1/buzz_species_treatment_m1.8.tiff",
  plot = p_trmt,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)





# Uncertainty graph


buzz_forest_data <- uncertainty_table %>%
  mutate(
    evidence = case_when(
      percent_low > 0  ~ "Increase supported",
      percent_high < 0 ~ "Decrease supported",
      TRUE             ~ "Uncertain"
    ),
    
    # Species with the greatest increase appear at the top
    species = factor(
      species,
      levels = species[order(percent_change)]
    )
  )

# Forest plot --------------------------------------------------------------

p_buzz_forest <- ggplot(
  buzz_forest_data,
  aes(
    x = percent_change,
    y = species,
    color = evidence
  )
) +
  
  # Reference value indicating no treatment effect
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey45"
  ) +
  
  # Approximate 95% confidence intervals
  geom_errorbar(
    aes(
      xmin = percent_low,
      xmax = percent_high
    ),
    orientation = "y",
    width = 0.14,
    linewidth = 0.7
  ) +
  
  # Percentage-change estimates
  geom_point(size = 2.7) +
  
  scale_color_manual(
    values = c(
      "Decrease supported" = "#D62728",
      "Increase supported" = "#228B22",
      "Uncertain" = "grey55"
    ),
    breaks = c(
      "Decrease supported",
      "Increase supported",
      "Uncertain"
    )
  ) +
  
  # Display species names in italics
  scale_y_discrete(
    labels = function(x) {
      parse(text = paste0("italic('", x, "')"))
    }
  ) +
  
  scale_x_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.03, 0.04))
  ) +
  
  labs(
    title = "Species-specific feeding-buzz responses to artificial light",
    subtitle = paste(
      "Points show conditional estimates;",
      "bars show approximate 95% confidence intervals"
    ),
    x = "Predicted change from dark to lit (%)",
    y = NULL,
    color = "Evidence"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 10,
      margin = margin(b = 12)
    ),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(
      size = 10,
      margin = margin(t = 8)
    ),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.width = unit(0.7, "cm"),
    plot.margin = margin(10, 15, 10, 10)
  )

p_buzz_forest

# create a summary table of total buzzes by species 

summary_table <- sm3_buzz_v2 %>%
  group_by(sp_clean) %>%
  summarise(
    total_buzzes = sum(t_buzzes)
  )



# margina Jday  -----------------------------------------------------------
# Create a lookup between species codes and scientific names
species_key <- sm3_buzz_v2 %>%
  distinct(sp_clean, sp_label) %>%
  mutate(sp = as.character(sp_clean))

# Convert the lookup table into a named vector
species_labels <- setNames(
  species_key$sp_label,
  species_key$sp
)

species_labels
# add the percent change to each sp 
species_labels_effect <- species_effects %>%
  mutate(
    sp = as.character(sp_clean),
    effect_text = sprintf(
      "%s (%+.1f%%)",
      species_labels[sp_clean],
      percent_change
    )
  ) %>%
  select(sp, effect_text) %>%
  deframe()



jday_seq <- seq(
  min(sm3_buzz_v2$jday_s, na.rm = TRUE),
  max(sm3_buzz_v2$jday_s, na.rm = TRUE),
  length.out = 100
)

pred_jday <- predictions(
  m1.8,
  newdata = datagrid(
    jday_s = jday_seq,
    trmt_bin = c(-1, 1),
    sp_clean = unique(sm3_buzz_v2$sp_clean),
    avg_moonlight_s = 0,
    nit_avg_temp_c_s = 0,
    nit_avg_wspm_s_s = 0,
    yr_s = 0,
    t_lepidoptera_s = 0
  ),
  type = "response"
)

jday_mean <- mean(sm3_buzz_v2$jday, na.rm = TRUE)
jday_sd   <- sd(sm3_buzz_v2$jday, na.rm = TRUE)

pred_jday <- pred_jday %>%
  mutate(
    jday = jday_s * (2 * jday_sd) + jday_mean,
    trmt = ifelse(trmt_bin == -1, "Dark", "Lit"),
    sp_label = str_replace(sp_clean, "([A-Z])[A-Z]+([A-Z]+)", "\\1. \\2")
  )

bm_jday_plot <- sm3_buzz_v2 %>%
  mutate(
    trmt = ifelse(trmt_bin == -1, "Dark", "Lit"),
    sp_label = str_replace(sp_clean, "([A-Z])[A-Z]+([A-Z]+)", "\\1. \\2")
  )


p_jday <- ggplot() +
  geom_point(
    data = bm_jday_plot,
    aes(x = jday, y = t_buzzes, color = trmt),
    alpha = 0.10,
    size = 0.7
  ) +
  geom_line(
    data = pred_jday,
    aes(x = jday, y = estimate, color = trmt),
    linewidth = 1.2
  ) +
  facet_wrap(~ sp_clean, scales = "free_y") +
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
  filename = "figures/glmm_v5/jday_v1.png",
  plot = p_jday,
  width = 10,
  height = 8,
  dpi = 300
)

# here we modify the graph to follow the asthetics from the treatment graph. 

pred_jday <- pred_jday %>%
  mutate(
    jday = jday_s * (2 * jday_sd) + jday_mean,
    
    trmt = factor(
      trmt_bin,
      levels = c(-1, 1),
      labels = c("Dark", "Lit")
    ),
    
    # Same species order as p_trmt
    sp = factor(
      as.character(sp_clean),
      levels = as.character(species_order)
    )
  ) %>%
  arrange(sp, trmt, jday)

bm_jday_plot <- sm3_buzz_v2 %>%
  mutate(
    trmt = factor(
      trmt_bin,
      levels = c(-1, 1),
      labels = c("Dark", "Lit")
    ),
    
    # Same species order as p_trmt
    sp = factor(
      as.character(sp_clean),
      levels = as.character(species_order)
    )
  )

p_jday <- ggplot() +
  geom_point(
    data = bm_jday_plot,
    aes(
      x = jday,
      y = t_buzzes,
      color = trmt
    ),
    alpha = 0.10,
    size = 0.7
  ) +
  
  geom_line(
    data = pred_jday,
    aes(
      x = jday,
      y = estimate,
      color = trmt,
      group = interaction(sp, trmt)
    ),
    linewidth = 1.2
  ) +
  
  facet_wrap(
    ~ sp,
    scales = "free_y",
    labeller = labeller(sp = species_labels)
  ) +
  
  scale_color_manual(
    values = c(
      "Dark" = "grey20",
      "Lit"  = "grey70"
    )
  ) +
  
  scale_y_continuous(
    trans = "log1p"
  ) +
  
  labs(
    x = "Julian day",
    y = "Predicted bat calls",
    color = "Treatment",
    title = ""
  ) +
  
  theme_minimal() +
  
  theme(
    strip.text = element_text(
      size = 12,
      face = "italic"
    ),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    legend.position = "bottom"
  )

p_jday


ggsave(
  filename = "figures/sm3_glmm_v1/jday_v1_fig5.tiff",
  plot = p_jday,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# rest of the marginal effects alternative. we might change it latter.

# Colors used across figures
trmt_colors <- c(
  "Dark" = "#3B3B3B",
  "Lit"  = "grey70"
)

# Use the central 95% of each predictor to avoid plotting sparse extremes
prediction_range <- function(x, n = 100) {
  limits <- quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  seq(limits[1], limits[2], length.out = n)
}


make_predictions <- function(variable,
                             values = NULL,
                             by_treatment = FALSE) {
  
  if (is.null(values)) {
    values <- prediction_range(sm3_buzz_v2[[variable]])
  }
  
  treatment_values <- if (by_treatment) c(-1, 1) else 0
  
  newdata <- expand_grid(
    predictor_value = values,
    trmt_bin = treatment_values
  ) %>%
    mutate(
      jday_s            = 0,
      avg_moonlight_s   = 0,
      nit_avg_temp_c_s  = 0,
      nit_avg_wspm_s_s  = 0,
      yr_s              = 0,
      t_lepidoptera_s   = 0,
      
      # Existing grouping levels are required even though REs are excluded
      site = sm3_buzz_v2$site[which(!is.na(sm3_buzz_v2$site))[1]],
      sp_clean = sm3_buzz_v2$sp_clean[
        which(!is.na(sm3_buzz_v2$sp_clean))[1]
      ]
    )
  
  # Replace the focal predictor's reference value
  newdata[[variable]] <- newdata$predictor_value
  newdata$predictor_value <- NULL
  
  predictions(
    m1.8,
    newdata = newdata,
    type = "response",
    re.form = NA
  ) %>%
    as.data.frame() %>%
    mutate(
      treatment = case_when(
        trmt_bin == -1 ~ "Dark",
        trmt_bin == 1  ~ "Lit",
        TRUE           ~ "Treatment average"
      )
    )
}


# Interactions with treatment
pred_date <- make_predictions(
  "jday_s",
  by_treatment = TRUE
)

pred_moon <- make_predictions(
  "avg_moonlight_s",
  by_treatment = TRUE
)

pred_year <- make_predictions(
  "yr_s",
  values = sort(unique(na.omit(sm3_buzz_v2$yr_s))),
  by_treatment = TRUE
)

# Significant main effects without treatment interactions
pred_temp <- make_predictions("nit_avg_temp_c_s")

pred_wind <- make_predictions("nit_avg_wspm_s_s")

pred_moths <- make_predictions("t_lepidoptera_s")


plot_interaction <- function(data, variable, x_label) {
  
  ggplot(
    data,
    aes(
      x = .data[[variable]],
      y = estimate,
      color = treatment,
      fill = treatment
    )
  ) +
    geom_ribbon(
      aes(
        ymin = conf.low,
        ymax = conf.high,
        group = treatment
      ),
      alpha = 0.18,
      color = NA
    ) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = trmt_colors) +
    scale_fill_manual(values = trmt_colors) +
    labs(
      x = x_label,
      y = "Predicted feeding buzzes",
      color = "Treatment",
      fill = "Treatment"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

# plot function 
plot_interaction <- function(data, variable, x_label) {
  
  ggplot(
    data,
    aes(
      x = .data[[variable]],
      y = estimate,
      color = treatment,
      fill = treatment
    )
  ) +
    geom_ribbon(
      aes(
        ymin = conf.low,
        ymax = conf.high,
        group = treatment
      ),
      alpha = 0.18,
      color = NA
    ) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = trmt_colors) +
    scale_fill_manual(values = trmt_colors) +
    labs(
      x = x_label,
      y = "Predicted feeding buzzes",
      color = "Treatment",
      fill = "Treatment"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
}


plot_main_effect <- function(data, variable, x_label) {
  
  ggplot(
    data,
    aes(
      x = .data[[variable]],
      y = estimate
    )
  ) +
    geom_ribbon(
      aes(
        ymin = conf.low,
        ymax = conf.high
      ),
      fill = "#377EB8",
      alpha = 0.20
    ) +
    geom_line(
      color = "#377EB8",
      linewidth = 1
    ) +
    labs(
      x = x_label,
      y = "Predicted feeding buzzes"
    ) +
    theme_classic(base_size = 11)
}

# individual panels
p_date <- plot_interaction(
  pred_date,
  "jday_s",
  "Julian date (standardized)"
)

p_moon <- plot_interaction(
  pred_moon,
  "avg_moonlight_s",
  "Mean moonlight intensity (standardized)"
)

p_temp <- plot_main_effect(
  pred_temp,
  "nit_avg_temp_c_s",
  "Nighttime temperature (standardized)"
)

p_wind <- plot_main_effect(
  pred_wind,
  "nit_avg_wspm_s_s",
  "Nighttime wind speed (standardized)"
)

p_moths <- plot_main_effect(
  pred_moths,
  "t_lepidoptera_s",
  "Moth abundance (standardized)"
)

if ("yr" %in% names(sm3_buzz_v2)) {
  
  year_key <- sm3_buzz_v2 %>%
    distinct(yr_s, yr) %>%
    drop_na() %>%
    arrange(yr_s)
  
  p_year <- plot_interaction(
    pred_year,
    "yr_s",
    "Year"
  ) +
    geom_point(size = 2.3) +
    scale_x_continuous(
      breaks = year_key$yr_s,
      labels = year_key$yr
    )
  
} else {
  
  p_year <- plot_interaction(
    pred_year,
    "yr_s",
    "Year (standardized)"
  ) +
    geom_point(size = 2.3)
}

p_all_covariates <-
  (p_date | p_moon) /
  (p_temp | p_wind) /
  (p_year | p_moths) +
  plot_annotation(
    title = "Environmental predictors of bat feeding-buzz activity",
    subtitle = paste(
      "Lines and points show population-level predictions;",
      "shaded areas show 95% confidence intervals"
    ),
    tag_levels = "A"
  ) &
  theme(
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold")
  )

p_all_covariates

ggsave("figures/sm3_glmm_v1/predictors_feddingbuzz.tiff", plot = p_all_covariates, width = 10, height = 8, dpi = 300)
