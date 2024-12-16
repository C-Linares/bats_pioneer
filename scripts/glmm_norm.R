# ---------------------------
##
## Script name:  glmm_norm
##
## Purpose of script: run models on the normalized data
##
## Author: Carlos Linares, \
## Date Created: 12/13/2024
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes: 
##   
##
## # inputs ------------------------------------------------------------------
#   data_for_analysis/prep_for_glmm/normalized_bm.csv  #product of the prep_for_glmm.R script. summary daily of bat call counts


# outputs ----------------------

# figures and model data


# libraries  --------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(lme4)
library(sjPlot)
library(ggeffects)
library(car)
library(glmmTMB)
library(corrplot)
library(effects)
library(reshape2)


# data --------------------------------------------------------------------


normalized_bm<-read_csv('data_for_analysis/prep_for_glmm/normalized_bm.csv')
head(normalized_bm)
normalized_bm$date<- ymd(normalized_bm$noche)


moon<-read.csv('data_for_analysis/moon_pred.csv')
moon$date<- as_date(moon$date)
moon.adj<-moon %>% mutate(
  phase = ifelse(above_horizon==FALSE,0,phase),
  fraction= ifelse(above_horizon==FALSE,0,fraction),
  l.illum= ifelse(above_horizon==FALSE,0,l.illum)
)
head(moon.adj)
moon.adj$date<-ymd(moon.adj$date)

# dayly average of l.illum

l.illum <- moon.adj %>% group_by(date) %>% summarise(l.illum = mean(l.illum))
l.illum$date<-ymd(l.illum$date)
head(l.illum)

weather<-read.csv('data_for_analysis/weather/nigh_averages.csv', header = T) #load nightly averages
weather$date<-as_date(weather$date)
head(weather)
weather$date<-ymd(weather$date)


#merge data by date
bm<-left_join(normalized_bm, l.illum, by='date')
bm<-left_join(bm, weather, by='date')
head(bm)



#correlation matrix just numeric variables

numeric_cols<- sapply(bm, is.numeric) # separate all the num col
cor1<-bm[,numeric_cols] #keeps just the numeric
c1<- cor(cor1,use="pairwise.complete.obs")
corrplot(c1, order= 'AOE')

# List of variable names to be scaled
variables_to_scale <- c(
  "jday",
  "l.illum",
  "avg_wind_speed",
  "avg_temperature"
)

# Loop over each variable, scale it, and assign it back to the data frame with a new name
for (var in variables_to_scale) {
  bm[[paste0(var, "_s")]] <- scale(bm[[var]], center = TRUE, scale = TRUE)
}

# make year between -1:1
bm <- bm %>%
  mutate(yr_s = case_when(
    year == 2021 ~ -1,
    year == 2022 ~ 0,
    TRUE ~ 1
  ))


# explore data ------------------------------------------------------------

ggplot(normalized_bm, aes(x = normalized_activity)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black") +
  facet_wrap(~ year) +  # Facet by year
  labs(
    title = "Histogram of Normalized Activity by Year",
    x = "Normalized Activity",
    y = "Frequency"
  ) +
  theme_minimal()

# models ------------------------------------------------------------------


# Fit a Mixed-Effects Model (GLM)

m1.1 <- lmer(
  normalized_activity ~ jday_s + I(jday_s ^ 2) + l.illum_s +
    avg_wind_speed_s + avg_temperature_s + yr_s + (1 |pair_group),
  data = bm
)

# Summary of the model
summary(m1.1)


plot_model(m1.1)

plot(residuals(m1.1))
ranef(m1.1)


m1.2 <- lmer(
  normalized_activity ~ jday_s + I(jday_s ^ 2) + l.illum_s + yr_s + (1 |pair_group)+(1 | sp),
  data = bm
)

summary(m1.2)
plot(residuals(m1.2))

AIC(m1.1,m1.2)




# glmm  -------------------------------------------------------------------


# Fit a zero-inflated beta GLMM
m2.1_proportional_zi <- glmmTMB(
  normalized_activity ~ jday_s + I(jday_s ^ 2) + l.illum_s +
    avg_wind_speed_s + avg_temperature_s + yr_s + (1 | pair_group),
  family = beta_family(link = "logit"),  # Beta distribution with logit link
  zi = ~1,  # Zero-inflation model
  data = bm
)

# Check the model summary
summary(m2.1_proportional_zi)


plot(residuals(m2.1_proportional_zi))
hist(residuals(m2.1_proportional_zi), breaks = 30, main = "Histogram of Residuals", xlab = "Residuals")
qqnorm(residuals(m2.1_proportional_zi))
qqline(residuals(m2.1_proportional_zi))

plot_models(m2.1_proportional_zi)
