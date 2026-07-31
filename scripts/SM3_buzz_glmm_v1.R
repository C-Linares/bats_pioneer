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

gtsummary::tbl_regression(m1.8, exponentiate = TRUE, 
                          label = list(
                            trmt_bin = "Treatment (Lit vs Dark)",
                            jday_s = "Julian Day (scaled)",
                             I(jday_s^2)= "Julian Day Squared (scaled)",
                            avg_moonlight_s = "Average Moonlight (scaled)",
                            nit_avg_temp_c_s = "Average Temperature (scaled)",
                            nit_avg_wspm_s_s = "Average Wind Speed (scaled)",
                            yr_s = "Year (scaled)",
                            t_lepidoptera_s = "Lepidoptera Abundance (scaled)"
                          )) 

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

# lets run a model without the correlation 
m1.9.1<-glmmTMB(
  #fixed effects
  n ~ trmt_bin + jday_s + I(jday_s^2) + avg_moonlight_s + nit_avg_temp_c_s +
    nit_avg_wspm_s_s + yr_s  + t_lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) || sp) +
    #interactions
    jday_s * trmt_bin + I(jday_s^2) * trmt_bin + trmt_bin*avg_moonlight_s  + trmt_bin*t_lepidoptera_s,  
  data = bm2,
  family = nbinom2(link = "log")
)

check_singularity(m1.9.1)
m1.9.1$sdr$pdHess
m1.9.1$fit$message
check_zeroinflation(m1.9.1)
calculate_c_hat(m1.9.1) 
# performance_mae(m1.9.1)
# range(bm2$n)
# hist(bm2$n, breaks = 500)
performance::r2(m1.9.1)
DHARMa::simulateResiduals(m1.9.1, n = 1000, plot = TRUE)
performance::check_collinearity(m1.9.1)
anova( m1.2, m1.3,m1.4,m1.5,m1.6, m1.8,m1.9, m1.9.1)
summary(m1.9.1)




