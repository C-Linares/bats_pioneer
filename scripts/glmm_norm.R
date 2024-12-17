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
if (!require("pacman")) install.packages("pacman")

# Use pacman to load libraries
pacman::p_load(tidyverse, magrittr, lme4, sjPlot, ggeffects, 
               car, glmmTMB, corrplot, effects, reshape2)

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

ggplot(bm, aes(x = normalized_activity)) +
  # geom_histogram(aes(y = ..density..), binwidth = 0.2, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_density(color = "blue", size = 1) +
  facet_wrap(~ year) +
  labs(
    title = "Histogram and Density Curve of Normalized Activity by Year",
    x = "Normalized Activity",
    y = "Density"
  ) +
  theme_minimal()

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

pseudo_r2 <- 1 - (m3.1$deviance / m3.1$null.deviance)
pseudo_r2

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

plot_model(m3.1_glmm)



library(effects)

# Plot marginal effects
all_effects <- allEffects(m3.1_glmm)
plot(all_effects)


library(sjPlot)

# Plot marginal effects
sjp.glmer(m3.1_glmm, type = "eff")


install.packages("lattice")
library(lattice)

# Plot random effects
dotplot(ranef(m3.1_glmm, condVar = TRUE), scales = list(relation = "free"))
ranef(m3.1_glmm)


trmt<- c("lit", "dark") 
trmt_bin_s <- c(-1,1)
trmt_bin_s

# correcting random eff
ran.efs <- ranef( m3.1_glmm )$sp # get the random effects 
randint<- ran.efs[,1]
randslope<- ran.efs[,2]

cint<-confint(m3.1_glmm)#[1:2,] # get fixed effects from the model
fixint<-cint[1,3]
fixslope<-cint[2,3]


#y = int + random.int[sp] + beta[1]treatment + random.slope[sp] * treatment

pred<- c(exp((randint+fixint)+(randslope+fixslope)*trmt_bin_s[1]),
         exp((randint+fixint)+(randslope+fixslope)*trmt_bin_s[2]))


abunddf <- data.frame(pred, trmt_bin, trmt_bin_s)

sp_names <- rownames(ran.efs)
sp_doubled <- rep(sp_names, each = 2)

# Check the length to match the number of rows in the dataframe
if (length(sp_doubled) == nrow(abunddf)) {
  # Add the new column to the dataframe
  abunddf$sp <- sp_doubled
} else {
  stop("The length of sp_doubled does not match the number of rows in abunddf.")
}

# View the updated dataframe
head(abunddf)

p7<-ggplot( abunddf, aes( x = trmt_bin, y = pred, colour = sp, group = sp, shape = sp) ) +  
  theme_classic( base_size = 12) +
  # scale_color_viridis_d(direction = -1, option = "A")+
  scale_shape_manual(values = 1:14) +  # Customize shapes for each 'sp', adjust as needed
  ylab( "bat calls" ) +
  xlab( "year" ) +
  geom_point(size=3)+
  geom_line()+
  # labels
  ylab("Bat calls") +
  xlab("") 
p7

#black and greys graph
p7.1 <- ggplot(abunddf, aes(x = trmt_bin, y = pred, colour = sp, group = sp, shape = sp)) +  
  theme_classic(base_size = 12) +
  scale_color_grey(start = 0, end = 1) +  # Use grayscale colors
  scale_shape_manual(values = 1:14) +  # Customize shapes for each 'sp'
  ylab("Bat calls") +
  xlab("Year") +
  geom_point(size = 3) +
  geom_line() +
  theme(
    plot.background = element_rect(fill = "black"),  # Set plot background to black
    panel.background = element_rect(fill = "black"),  # Set panel background to black
    axis.text = element_text(color = "white"),  # Set axis text to white
    axis.title = element_text(color = "white"),  # Set axis titles to white
    legend.background = element_rect(fill = "black"),  # Set legend background to black
    legend.text = element_text(color = "white"),  # Set legend text to white
    legend.title = element_text(color = "white")  # Set legend title to white
  )

p7.1