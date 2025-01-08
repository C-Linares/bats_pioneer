# ---------------------------
##
## Script name:  glmm_norm
##
## Purpose of script: run models on the normalized data. Data that has been normalized by dividing the calls in lit sites by the calls in dark sites. 
##
## Author: Carlos Linares, 
## Date Created: 12/13/2024
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes: 
## 
## 
##   
##
## # inputs ------------------------------------------------------------------

# -data_for_analysis/prep_for_glmm/normalized_bm.csv  #product of the prep_for_glmm.R script. summary daily of bat call counts


# outputs ----------------------

# figures and model data


# libraries  --------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")

# Use pacman to load libraries
pacman::p_load(tidyverse, magrittr, lme4, sjPlot, ggeffects, 
               car, glmmTMB, corrplot, effects, reshape2)

# save past images

load('working_env/glmm_norm.RData')


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

ggplot(bm, aes(x = normalized_activity)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black", alpha=.5) +
  facet_wrap(~ year) +  # Facet by year
  labs(
    title = "Histogram of Normalized Activity by Year",
    x = "Normalized Activity",
    y = "Frequency"
  ) +
  theme_minimal()

ggplot(bm, aes(x = normalized_activity, colour = sp)) +
  # geom_histogram(aes(y = ..density..), binwidth = 0.2, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_density(color = "blue", size = 1) +
  facet_wrap(~ year) +
  labs(
    title = "Histogram and Density Curve of Normalized Activity by Year",
    x = "Normalized Activity",
    y = "Density"
  ) +
  theme_minimal()

ggplot(bm, aes(x = normalized_activity, color = sp)) +
  # geom_histogram(aes(y = ..density..), binwidth = 0.2, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_density(size = 1) +
  facet_wrap(~ year) +
  labs(
    title = "Histogram and Density Curve of Normalized Activity by Year and Species",
    x = "Normalized Activity",
    y = "Density",
    color = "Species"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggplot(bm, aes(x = bin_act)) +
  # geom_histogram(aes(y = ..density..), binwidth = 0.2, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_density(color = "blue", size = 1) +
  facet_wrap(~ year) +
  labs(
    title = "Histogram and Density Curve of if exp>contrl=1 by Year",
    x = "Normalized Activity",
    y = "Density"
  ) +
  theme_minimal()

# Summarize the data to count the number of 1's and 0's for each year
summary_data <- bm %>%
  group_by(year, bin_act, sp) %>%
  summarize(count = n(), .groups = 'drop')

# View the summarized data
print(summary_data)

# Create the box plot
ggplot(summary_data, aes(x = factor(year), y = count)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  labs(title = "Count of bin_act by Year",
       x = "Year",
       y = "Count",
       fill = "bin_act") +
  scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()



# models ------------------------------------------------------------------



# logistic regression -----------------------------------------------------


library(broom)
library(ResourceSelection)

m3.1 <- glm(
  bin_act ~ jday_s + I(jday_s^2) + yr_s + avg_wind_speed_s + avg_temperature_s + l.illum_s ,
    data = bm,
  family = binomial
)
summary(m3.1)

# Calculate c-hat for GLM
residual_deviance_glm <- deviance(m3.1)
residual_df_glm <- df.residual(m3.1)
c_hat_glm <- residual_deviance_glm / residual_df_glm
print(c_hat_glm)
AIC(m3.1)

pseudoR2 <- (m3.1$null.deviance - m3.1$deviance) / m3.1$null.deviance  # I got this formula from http://r.qcbs.ca/workshop06/book-en/binomial-glm.html

pseudo_r2 # model explains nothing...

# Plot residuals
par(mfrow = c(2, 2))
plot(m3.1)

# Extract coefficients and confidence intervals
coef_df <- tidy(m3.1, conf.int = TRUE)

# Plot coefficients with ggplot2
ggplot(coef_df, aes(x = term, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(title = "Model Coefficients with Confidence Intervals", x = "Predictor", y = "Estimate") +
  theme_minimal() 

  # coord_flip()

plot_model(m3.1)



#glmm

# glmm binomial -----------------------------------------------------------

# Fit the GLMM model
m3.1_glmm <- glmer(
  bin_act ~ jday_s + I(jday_s^2) + yr_s + avg_wind_speed_s + avg_temperature_s + l.illum_s +
    (1 | pair_group) + (1+yr_s+jday_s+I(jday_s^2) | sp),  # Specify random effects
  data = bm,
  family = binomial
)


# Summary of the GLMM model
summary(m3.1_glmm)


# Calculate c-hat for glmm
residual_deviance_glm <- deviance(m3.1_glmm)
residual_df_glm <- df.residual(m3.1_glmm)
c_hat_glm <- residual_deviance_glm / residual_df_glm
print(c_hat_glm)

AIC(m3.1,m3.1_glmm)

install.packages("MuMIn")
install.packages("performance")
library(lme4)
library(MuMIn)
library(performance)

# Calculate R-squared values
r_squared <- r.squaredGLMM(m3.1_glmm) # seems like the model does not explain much of the variance in the data. 
print(r_squared)

plot_model(m3.1_glmm)
plot_model(m3.1_glmm, type = "re")



trmt<- c("lit", "dark") 
trmt_bin_s <- c(-1,1)
trmt_bin_s

# random effects 
ran.efs <- ranef( m3.1_glmm )$sp # get the random effects 
randint<- ran.efs[,1]
randslope<- ran.efs[,2]

# confidence intervals fix eff
cint<-confint(m3.1_glmm)#[1:2,] # get fixed effects from the model
fixint<-cint[1,3]
fixslope<-cint[2,3]

# jday random effects plot. 

n=100
sl=100
ones = rep(1,100)

# obs. values
ord.day <- seq(min(bm[, "jday"]), max(bm[, "jday"]), length.out = n)
# extraer el jday sqr 
ord.day.sqr<-ord.day^2
#std pred
jday.s <- scale(ord.day)
jday.sqr<- scale(ord.day.sqr)


#extract fixed coef jday
#pull out random effects at the sp level #
ran.efs <- ranef( m3.1_glmm )$sp
ranefs<- ran.efs[, c(1,3,4)]

#pull out fixed effects
fix.efs <- fixef( m3.1_glmm )
#view
fix.efs



rss <- ranefs
rss[, 1] <- rss[, 1] + fix.efs[1]
rss[, 2] <- rss[, 2] + fix.efs[2]# keep doing this for each of the random effects. 
rss[, 3] <- rss[, 3] + fix.efs[3]

rss # we added te fix to the random effe

a<-t(rss)
a

b<-t( cbind( ones, jday.s, jday.sqr))
b

indpred<- exp( as.matrix(a) %*% as.matrix(b) )

abunddf <- data.frame(t(indpred), jday.s, ord.day)



# save the environment ---------------------------------------------------


#save
save.image('working_env/glmm_norm.RData')
