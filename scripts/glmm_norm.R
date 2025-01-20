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

# -data_for_analysis/prep_for_glmm/normalized_bm.csv  #product of the prep_for_glmm.R script. activity normalized by control site 


# outputs ----------------------

# figures and model data from normalized activity


# libraries  --------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")

# Use pacman to load libraries
pacman::p_load(tidyverse, magrittr, lme4, sjPlot, ggeffects, 
               car, glmmTMB, corrplot, effects, reshape2, DHARMa)

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




# correlation -------------------------------------------------------------

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


for (var in variables_to_scale) {
    # Check if the column exists and is numeric
    if (var %in% names(bm) && is.numeric(bm2[[var]])) {
      mean_val <- mean(bm[[var]], na.rm = TRUE)
      sd_val <- sd(bm[[var]], na.rm = TRUE)
      bm[[paste0(var, "_s")]] <- (bm[[var]] - mean_val) / (2 * sd_val)
    } else {
      warning(paste("Column", var, "is not numeric or does not exist in the data frame."))
    }
  }
  

# make year between -1:1
bm <- bm %>%
  mutate(yr_s = case_when(
    year == 2021 ~ -1,
    year == 2022 ~ 0,
    TRUE ~ 1
  ))


# Check for missing values
sum(is.na(bm$j_diff))       # Number of missing values in j_diff
sum(is.na(bm$jday_s))       # Number of missing values in jday_s
sum(is.na(bm$yr_s))         # Number of missing values in yr_s
sum(is.na(bm$l.illum_s))    # Number of missing values in l.illum_s
sum(is.na(bm$pair_group))   # Number of missing values in pair_group



# explore data ------------------------------------------------------------

# Visualize the distribution of j_diff
hist(bm$j_diff, main="Distribution of j_diff", xlab="j_diff")

hist(bm$jday_s, main="Distribution of jday_s", xlab="jday_s")

hist(bm$l.illum_s, main="Distribution of l.illum_s", xlab="l.illum_s")

hist(bm$normalized_activity, main="Distribution of normalized_activity", xlab="normalized_activity")

hist(bm$bin_act, main="Distribution of bin_act", xlab="bin_act")


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


plot(m3.1_glmm)

# Simulate residuals using DHARMa for GLMM
sim_residuals <- simulateResiduals(m3.1_glmm)
plotSimulatedResiduals(sim_residuals)

t1 <- testUniformity(sim_residuals) # model is not uniform...
t2<-testZeroInflation(sim_residuals)
t3<-testDispersion(sim_residuals)
# model jdiff -------------------------------------------------------------

m2<- glm(
  j_diff ~ jday_s + I(jday_s^2) + yr_s + l.illum_s,  # Specify random effects
  data = bm
)

summary(m2)

m3.1 <- lmer(
  j_diff ~ jday_s + I(jday_s^2) + yr_s + l.illum_s  + 
    (1 | pair_group) + (1|sp),  # Specify random effects
  data = bm
)

update(m3.1, ~ . + (1 +yr_s| sp)) # update model to include year in the random slope

summary(m3.1)

fixef(m3.1)
ranef(m3.1)

plot(m3.1)


AIC(m3.1,m2,m3.1_glmm) # compare models

plot_model(m3.1)





# model normalize activity ------------------------------------------------


# Adjust values at the boundaries
bm$normalized_activity[bm$normalized_activity == 0] <- 1e-6 
bm$normalized_activity[bm$normalized_activity == 1] <- 1 - 1e-6 

# Fit the Beta GLMM
m4_beta <- glmmTMB(
  normalized_activity ~ jday_s + I(jday_s^2) + yr_s + l.illum_s + 
    (1 | pair_group) + (1 | sp),
  data = bm,
  family = beta_family(link = "logit")
)

m4_beta

summary(m4_beta)

m4.1 <- glmmTMB(
  normalized_activity ~ jday_s + I(jday_s^2) + yr_s + l.illum_s + 
    (1 | pair_group) + (1 + yr_s | sp),
  data = bm,
  family = betabinomial(link = "logit")
)
m4.1
summary(m4.1)


# the previous model indicate that the fit is singular. 

# Calculate c-hat for m4_beta
residual_deviance_glm <- deviance(m4_beta) # there's no deviance for a beta family model. 
# residual_df_glm <- df.residual(m3.1_glmm)
# c_hat_glm <- residual_deviance_glm / residual_df_glm
# print(c_hat_glm)

plot(fitted(m4_beta), residuals(m4_beta)) # from this plot is seems like the model does not fit well.

plot(density(residuals(m4_beta))) # check the distribution of the residuals. bimodal

#beta binomial model 
m4.2<- glmmTMB(
  normalized_activity ~ jday_s + I(jday_s^2) + yr_s + l.illum_s + 
    (1 | pair_group) + (1 | sp),
  data = bm,
  family = beta_binomial(link = "logit")
)

# above model failed...
summary(m4.2)

# marginal effect plots ---------------------------------------------------




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
