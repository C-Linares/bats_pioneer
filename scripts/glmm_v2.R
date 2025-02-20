# ---------------------------
##
## Script name:  glmm_v2.R
##
## Purpose of script: Running models for the bat paper. 
## v2: updated data for weather from craters of the moon, bat traits from the ratio ear/arm and new moon data calculated with the moonlit package. 
##
## Author: Carlos Linares, Jen Cruz (collaborator)
##
## Date Created: 2/11/2025
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes: 
##   
##
## ---------------------------
## # inputs ------------------------------------------------------------------
#   data_for_analysis/prep_for_glmm/bm.csv product of the prep_for_glmm.R script. summary daily of bat call counts

# bm.ai<- read_csv('data_for_analysis/') # activity index for the bat species.


# outputs ----------------------

# Marginal effect plots and random effects plots.



# libraries  --------------------------------------------------------------


if (!require("pacman")) install.packages("pacman")

# Use pacman to load libraries
pacman::p_load(
  "tidyverse",
  "magrittr",
  "lme4",
  "sjPlot",
  "ggeffects",
  "car",
  "glmmTMB",#model
  "corrplot",
  "effects",
  "reshape2",
  "DHARMa",#check assumptions 
  "marginaleffects", #plot models
  "MuMIn", #models
  "performance",
  "viridis",
  "data.table"
)

#load environment 
load(file = "working_env/glmm_v2.RData") 

#load data ---------------------------------------------------------------

bm <- read_csv('data_for_analysis/prep_for_glmm/bm.csv')

filtered_bm <- bm[!(bm$AUTO.ID. %in% c("Noise", "NoID")), ] # Noise needs to be filter out before analysis other wise there are more records in dark areas.
# rename col
colnames(filtered_bm)[2] <- "sp" # change from AUTO.ID to sp


#load activity index

bm.ai<- read_csv('data_for_analysis/prep_for_glmm/bm.miller.day.csv') # activity index for the bat species.
filtered_bm.ai<-bm.ai[!(bm.ai$sp %in% c("Noise", "NoID")), ] # Noise needs to be filter out before analysis
colnames(filtered_bm.ai)
filtered_bm.ai$...1<-NULL # remove unnecessary variable. 
# predictors --------------------------------------------------------------


# bat traits arm/ear ratio

btrait<-read.csv('data_for_analysis/bat_traits/bat_eararm_v2.csv') # load bat traits')

# Creaters weather (night).

crmo.wet.night<-read_csv("data_for_analysis/weather/craters_weater/craters_night.csv") # load creaters night weather. 


# Moon

moon.int.night<-fread("data_for_analysis/moon_pred/moon.int.night.csv")
moon.int.night$date<-as_date(moon.int.night$date)

# insects 

c_bugs<-read_csv("data_for_analysis/insect_wranglin/c_bugs.csv") # load insect data.
colnames(c_bugs)[3]<-"yr" # rename site

# merge -------------------------------------------------------------------

# activity index with bat data. 

bm2 <- left_join(filtered_bm, filtered_bm.ai[, c("sp", "site", "noche", "activity_min")], by = c("sp", "site", "noche"))

# merge with traits
bm2<- left_join(bm2, select(btrait, six_sp, ear.arm), by = c("sp" = "six_sp"))

summary(bm2)

# merge with weather

bm2<- left_join(bm2, crmo.wet.night, by=c("noche"="date"))


# merge with moon
# summarize moon.int.night by date.
m.moon.int<- moon.int.night %>%
  group_by(date) %>%
  summarise(
    moonlight = mean(moonlightModel, na.rm = TRUE),
    mphase = mean(moonPhase, na.rm = TRUE),
    twilight = mean(twilightModel, na.rm = TRUE),
    tillum = mean(illumination, na.rm = TRUE)
  )


#
bm2<- left_join(bm2, m.moon.int, by=c("noche"="date")) # merge with moon data



# Identify rows with NA values
rows_with_na <- bm2[rowSums(is.na(bm2)) > 0, ] # some wind and temp have NAs because we are missing august 2021 weather data.

# make the rows with NAs to 0 
bm2<- bm2 %>%
  mutate(
    moonlight = ifelse(is.na(moonlight), 0, moonlight),
    mphase = ifelse(is.na(mphase), 0, mphase),
    twilight = ifelse(is.na(twilight), 0, twilight),
    tillum = ifelse(is.na(tillum), 0, tillum),
  )

summary(bm2)

# merge with insects

# calculate week for bm2 
bm2$wk<-week(bm2$noche) # week of the year

bm2<- left_join(bm2, c_bugs, by = c("site", "wk", "yr")) # merge with insects data

# check for NAs
summary(bm2)  

#make NAs the 0 for t.insect, and t. lepidoptera given we don't have that data.

bm2<- bm2 %>%
  mutate(
    t.insect = ifelse(t.insect==0, NA, t.insect),
    t.lepidoptera = ifelse(t.lepidoptera==0, NA, t.lepidoptera)
  )

# we have several weeks in 2023 where there is no data 5659 lines.


# correlation  ------------------------------------------------------------
# before modelling we have to check for correlation between the predictors as VIF is not adequate for negative binomial models.
# check for correlation 

numeric_cols<- sapply(bm2, is.numeric) # separate all the num col
cor1<-bm2[,numeric_cols] #keeps just the numeric

c1<- cor(cor1,use="pairwise.complete.obs")
corrplot(c1, order= 'AOE')

# from the plot it seems like the tmp and wind speed, twilight, total illumination, moon phase and moonlight are correlated not to be included in the model at the same time. Insect variables like total lepidoptera and and total insects are also correlated with moon phase and light but less than 0.4

# standardize predictors ---------------------------------------------
# calculate jday
bm2$jday<-yday(bm2$noche)

# List of variable names to be scaled
variables_to_scale <- c(
  "ear.arm",
  "nit_avg_tempC",
  "nit_avg_wspm.s",
  "moonlight",
  "mphase",
  "twilight",
  "tillum",
  "jday",
  "t.insect",
  "t.lepidoptera"
)

 
# # Loop over each variable, scale it by dividing by two standard deviations, and assign it back to the data frame with a new name
for (var in variables_to_scale) {
  # Check if the column exists and is numeric
  if (var %in% names(bm2) && is.numeric(bm2[[var]])) {
    mean_val <- mean(bm2[[var]], na.rm = TRUE)
    sd_val <- sd(bm2[[var]], na.rm = TRUE)
    bm2[[paste0(var, "_s")]] <- (bm2[[var]] - mean_val) / (2 * sd_val)
  } else {
    warning(paste("Column", var, "is not numeric or does not exist in the data frame."))
  }
}

summary(bm2) # visual check for new standardize variables.
# make treatment -1 to 1

bm2$trmt_bin <- ifelse(bm2$trmt_bin == 1, 1, -1)

# make year between -1:1
bm2 <- bm2 %>%
  mutate(yr_s = case_when(
    yr == 2021 ~ -1,
    yr == 2022 ~ 0,
    TRUE ~ 1
  ))

cattle <- c( # sites with cattle
  "long01",
  "long02",
  "long03" ,
  "long04",
  "long05",
  "vizc01",
  "vizc03" ,
  "vizc04",
  "vizc02"
)

bm2$moo<- ifelse(bm2$site %in% cattle, 1, -1)

# explore data ------------------------------------------------------------

summary(bm2)

# Plot the distribution of the count data
ggplot(bm2, aes(x = n)) +
  geom_histogram(binwidth = 40, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Bat Calls", x = "Number of Calls", y = "Frequency")

ggplot(bm2, aes(x=jday, y=n, col=treatmt))+
  geom_point()+
  facet_wrap(~site)+
  theme_blank()+
  scale_color_manual(values = c("#0033A0", "#D64309"))+
  labs(title = "Bat acoustic activity 2021-2023",
       x = "Julian Day",
       y = "n calls",
       color = "Treatment")
  
ggplot(bm2, aes(x=jday, y=n, col=treatmt))+
  geom_point()+
  facet_wrap(~sp, scales = "free_y")+
  theme_blank()+
  scale_color_manual(values = c("#0033A0", "#D64309"))+
  labs(title = "Bat acoustic activity 2021-2023",
       x = "Julian Day",
       y = "n calls",
       color = "Treatment")


# models ------------------------------------------------------------------

# Negative Binomial -------------------------------------------------------

m1.1nb <- glmmTMB(
  n ~ trmt_bin + jday_s + I(jday_s^2) + moonlight_s +
    nit_avg_wspm.s + yr_s +
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp) +
    jday_s * trmt_bin + I(jday_s^2) * trmt_bin + yr_s * trmt_bin,
  data = bm2,
  family = nbinom2(link = "log")
)

summary(m1.1nb)
plot_model(m1.1nb)
residual_deviance <- deviance(m1.1nb)
residual_df <- df.residual(m1.1nb)

# Calculate c-hat using residual deviance
c_hat_deviance <- residual_deviance / residual_df
print(c_hat_deviance)


# using a poisson model did not converged...
m1.1poisson <- glmmTMB(
  n ~ trmt_bin + jday_s + I(jday_s^2) + moonlight_s +
    nit_avg_wspm.s + yr_s +
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp) +
    jday_s * trmt_bin + I(jday_s^2) * trmt_bin + yr_s * trmt_bin,
  data = bm2,
  family = poisson(link = "log")
)



# Simulate residuals using DHARMa for GLMM
sim_residuals <- simulateResiduals(m1.1nb, plot = T,quantreg=T)

t1 <- testUniformity(sim_residuals) # model deviates from the expected distribution. What can we do about it? Jen said may not matter much. 
t2 <- testZeroInflation(sim_residuals)
t3 <- testDispersion(sim_residuals)
t4 <- testOutliers(sim_residuals) # model has some outliers... We have to find out where the outliers are coing from methodology or actual data or some other shit. 

# play with the comparison function in the marginal effects plot. 

# Calculate R-squared values

r_squared <- performance::r2(m1.1nb) 
print(r_squared)


# m1.2nb Model ear/arm ratio and insects 

m1.2nb <- glmmTMB(
  #fixed effects
  n ~ trmt_bin + jday_s + I(jday_s^2) + moonlight_s +
    nit_avg_wspm.s_s + yr_s + ear.arm_s + t.lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp) +
    #interactions
    jday_s * trmt_bin + I(jday_s^2) * trmt_bin + yr_s * trmt_bin,
  data = bm2,
  family = nbinom2(link = "log")
)

# Output model summary
summary(m1.2nb)

# Plot model
plot_model(m1.2nb)

# Calculate and print c-hat using residual deviance
c_hat_deviance <- deviance(m1.2nb) / df.residual(m1.2nb)
print(c_hat_deviance)

# Simulate residuals using DHARMa for GLMM
sim_residuals <- simulateResiduals(m1.2nb, plot = TRUE, quantreg = TRUE)

# Perform diagnostic tests on residuals
testUniformity(sim_residuals)
testZeroInflation(sim_residuals)
testDispersion(sim_residuals)
testOutliers(sim_residuals)

# Calculate and print R-squared values
r_squared <- performance::r2(m1.2nb)
print(r_squared)


# According to Allison I could just run a model with the light and that should be it. # tried that below and also tried adding moon as an interacting term. it did not improve the fit. 
m1.3nb <- glmmTMB(
  #fixed effects
  n ~ trmt_bin + jday_s + I(jday_s^2) + yr_s +
    (1 | site) + (1 + trmt_bin | sp) ,
  data = bm2,
  family = nbinom2(link = "log")
)
summary(m1.3nb)
plot_model(m1.3nb)
AIC(m1.3nb, m1.2nb,m1.1nb) # AIC is lower for the model with the insects and ear/arm ratio.

# Simulate residuals using DHARMa for GLMM
sim_residuals <- simulateResiduals(m1.3nb, plot = TRUE, quantreg = TRUE)

# Perform diagnostic tests on residuals
testUniformity(sim_residuals)
testZeroInflation(sim_residuals)
testDispersion(sim_residuals)
testOutliers(sim_residuals)

# Calculate and print R-squared values
r_squared <- performance::r2(m1.3nb)
print(r_squared)


# model with the light: insect interaction

m1.4nb <- glmmTMB(
  #fixed effects
  n ~ trmt_bin + jday_s + I(jday_s^2) + moonlight_s +
    nit_avg_wspm.s_s + yr_s + ear.arm_s + t.lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp) +
    #interactions
    jday_s * trmt_bin + I(jday_s^2) * trmt_bin + yr_s * trmt_bin + trmt_bin*t.lepidoptera_s + trmt_bin*moonlight_s,
  data = bm2,
  family = nbinom2(link = "log")
)
summary(m1.4nb)


# Marginal effect plots ---------------------------------------------------


# marginal with marginal effects package. 

p1.1<-plot_predictions(m1.1nb,                # The model we fit
                 type = "response",     # We'd like the predictions on the scale of the response
                 conf_level = 0.95,     # With a 95% confidence interval
                 condition = c("trmt_bin", "yr_s"), # Plot predictions for "l.illum_s" while holding all others at their mean
                 vcov = TRUE) +         # Compute and display the variance-covariance matrix
  xlab("") +            # Labels for axes
  ylab("bat calls") +
  theme_bw()                             # White background

p1.1 <- plot_predictions(m1.1nb,                # The model we fit
                         type = "response",     # We'd like the predictions on the scale of the response
                         conf_level = 0.95,     # With a 95% confidence interval
                         condition = c("trmt_bin", "yr_s"), # Plot predictions for "l.illum_s" while holding all others at their mean
                         vcov = TRUE) +         # Compute and display the variance-covariance matrix
  xlab("") +            # Labels for axes
  ylab("bat calls") +
  theme_bw() +          # White background
  scale_color_manual(values = c("-1" = "blue", "0" = "green", "1" = "red"), 
                     labels = c("-1" = "2021", "0" = "2022", "1" = "2023"), 
                     name = "Year") +
  scale_fill_manual(values = c("-1" = "blue", "0" = "green", "1" = "red"), 
                    labels = c("-1" = "2021", "0" = "2022", "1" = "2023"), 
                    name = "Year")

p1.2<-plot_predictions(m1.1nb,                # The model we fit
                 type = "response",     # We'd like the predictions on the scale of the response
                 conf_level = 0.95,     # With a 95% confidence interval
                 condition = c("jday_s","trmt_bin"), 
                 vcov = TRUE) +         # Compute and display the variance-covariance matrix
  xlab("") +            # Labels for axes
  ylab("bat calls") +
  theme_bw()        


# marginal effect for jday by sp
pred2 <- predictions(m1.1nb,
                     newdata = datagrid(sp = bm2$sp,
                                        jday_s = seq(min(bm2$jday_s), max(bm2$jday_s), length.out = 100)))
print(pred2)
p2 <- ggplot(pred2, aes(jday_s, estimate, color = sp)) +
  geom_line() +
  labs(y = "bat calls", x = "jday", title = "Quadratic growth model")+
  facet_wrap(~sp, scales = "free_y")
p2

p2 <- ggplot(pred2, aes(x = jday_s, y = estimate, color = sp, linetype = sp)) +
  geom_line(linewidth = 1) +  # Adjust line size for better visibility
  scale_color_viridis(discrete = TRUE, option = "D") +  # Use a colorblind-friendly palette
  labs(y = "Bat calls", x = "Julian day", title = "Marginal effect plot jday by sp") +
  theme_minimal() +  # Optional: adding a theme for better aesthetics
  theme(legend.title = element_blank())  # Remove legend title for cleaner appearance
p2

plot_marginal_effects <- function(model, x_var, group_var = NULL, title = NULL) {
  pred <- predictions(model, newdata = datagrid(!!x_var := unique(bm2[[x_var]]),
                                                sp = unique(bm2$sp)))
  
  p <- ggplot(pred, aes_string(x = x_var, y = "estimate", color = "sp")) +
    geom_line() +
    labs(y = "Bat calls", x = x_var, title = title) +
    scale_color_viridis(discrete = TRUE, option = "D") +
    theme_minimal()
  
  if (!is.null(group_var)) {
    p <- p + facet_wrap(as.formula(paste("~", group_var)), scales = "free_y")
  }
  return(p)
}


p2 <- plot_marginal_effects(m1.9nb, "jday_s", "sp", "Marginal effect of jday_s by species")

# marginal effect for year by sp

pred3 <- predictions(m1.1nb,
                     newdata = datagrid(sp = bm2$sp,
                                        yr_s = unique(bm2$yr_s)))

p3<- ggplot(pred3, aes(x =yr_s, y= estimate))+
              geom_point()+
  geom_errorbar( aes( ymin = conf.low, ymax = conf.high ) )+ 
  facet_wrap(~sp, scales="free")
  
            
p3

pred4<- predictions(m1.1nb,
                     newdata = datagrid(sp = bm2$sp,
                                        trmt_bin = unique(bm2$trmt_bin)))
p4 <- ggplot(pred4, aes(x = trmt_bin, y = estimate, col=sp)) +
  geom_point() +
  geom_line()
p4

p4+ facet_wrap(~sp)

preds<- predictions(m1.1nb)
preds$tmt<-ifelse(preds$trmt_bin==1, "light", "dark")
# now add the years -1=2021, 0=2022, 1=2023
preds$yr_s<-ifelse(preds$yr_s==-1, "2021", ifelse(preds$yr_s==0, "2022", "2023"))

p5.1<- ggplot(preds, aes(x = jday_s, y = estimate, col=tmt))+ # not sure what is happening here. 
  geom_line()+
  facet_wrap(~sp, scales = "free_y")
p5.1

p5<-ggplot(preds, aes(tmt, estimate, col=yr_s)) +
  geom_violin() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = .75), size=0.4, alpha=0.9) +  #The jitter.width argument controls the amount of jitter, and the dodge.width argument controls the spacing between the categories (treatments and years).
  facet_wrap(~sp, scales = "free_y")+
  labs(y="estimate bat calls", x="", title="Boxplot of Estimated bat calls by Treatment, Year, and Species")+
  scale_color_viridis(discrete = TRUE, option = "H")   # Use the viridis color palette
  # scale_colour_brewer(palette = "Spectral")+
  # theme_minimal()   # Optional: adding a theme for better aesthetics
  # theme(legend.position = "none")  # Optional: remove legend if not needed
p5


print(predictions(model = m1.9nb, 
                  newdata = datagrid( yr_s= c(-1, 0, 1),
                                      sp= unique(bm2$sp)),
                  conf_level = 0.95))

p6<- ggplot(preds, aes(x = jday_s, y = estimate, col=sp))+ # not sure what is happening here. 
  geom_point()+
  ylab("predicted bat calls")+
  scale_color_viridis(discrete = TRUE, option = "H") +  # Use the viridis color palette
  facet_grid(tmt~yr_s)
p6

# Create a list of plots and their corresponding filenames
figures_dir <- "figures/glmm_v1/marginleffects_out"


plots <- list(p1.1, p1.2, p2, p3, p4, p5.1, p5, p6)
filenames <- c("p1.1.png", "p1.2.png", "p2.png", "p3.png", "p4.png", "p5.1.png", "p5.png", "p6.png")

# Save all plots
lapply(seq_along(plots), function(i) {
  ggsave(filename = file.path(figures_dir, filenames[i]), plot = plots[[i]], device = "tiff", width = 15, height = 10, units = "cm")
})

library(sjPlot)
a1<-tab_model(m1.9nb, show.icc = TRUE, show.aic = TRUE, show.re.var = TRUE)
a1
plot_predictions(m1.9nb, condition = c("jday_s", "trmt_bin"), type = "response", vcov = T) +
  theme_minimal() +
  labs(y = "Predicted Bat Calls", x = "Julian Day", color = "Treatment")









# Marginal effect plots m1.2nb --------------------------------------------

# Below is marginal effect plot for moon illumination. 

p2.1<-plot_predictions(m1.2nb,                # The model we fit
                       type = "response",     # We'd like the predictions on the scale of the response
                       conf_level = 0.95,     # With a 95% confidence interval
                       condition = c("moonlight_s"), # Plot predictions 
                       vcov = TRUE) +         # Compute and display the variance-covariance matrix
  xlab("") +            # Labels for axes
  ylab("bat calls") +
  theme_bw()                             # White background
p2.1

# Create a data grid for predictions to replicate plot above.
pred2.1 <- predictions(m1.2nb, 
                       newdata = datagrid(moonlight_s = seq(min(bm2$moonlight_s), max(bm2$moonlight_s), 
                                                            length.out = 100)))
# Calculate non-standardized moonlight for graphics moonlight_s x (2xsd)+mean
sdmoon<-sd(bm2$moonlight)
mmoon<-mean(bm2$moonlight)
pred2.1$moonlight<-pred2.1$moonlight_s*(2*sdmoon)+mmoon


p2.1<-ggplot(pred2.1, aes(x = moonlight, y = estimate)) +
     geom_line() +
     geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
     labs(y = "Bat Calls", x = "Moonlight", title = "Marginal Effect of Moonlight") +
     theme_minimal() + 
     scale_color_viridis(discrete = TRUE, option = "D") +  # Use a colorblind-friendly palette
     theme(legend.title = element_blank())  # Remove legend title for cleaner appearance
p2.1

# Marginal effect plot for wind speed.
pred2.2<- predictions(m1.2nb, 
                       newdata = datagrid(nit_avg_wspm.s_s = seq(min(bm2$nit_avg_wspm.s_s), max(bm2$nit_avg_wspm.s_s), 
                                                            length.out = 100)),
                      re.form = NA) # re.form = NA to get the marginal effect of wind speed without random effects.)

# Calculate non-standardized wind speed for graphics
pred2.2$wspm.s<-pred2.2$nit_avg_wspm.s_s*(2*sd(bm2$nit_avg_wspm.s))+mean(bm2$nit_avg_wspm.s) # recalculate the non-standardized wind speed for graphics.

p2.2 <- ggplot(pred2.2, aes(x = wspm.s, y = estimate, ymin= conf.low, ymax=conf.high)) +
        geom_line()+
        geom_ribbon( alpha = 0.2)+
        labs(y = "Bat Calls", x = "Wind Speed", title = "Marginal Effect of Wind Speed") +
        theme_minimal() +
        scale_color_viridis(discrete = TRUE, option = "D") +  # Use a colorblind-friendly palette
        theme(legend.title = element_blank())+  # Remove legend title for cleaner appearance
p2.2


# Marginal effect plot for year. 

t<-plot_predictions(m1.2nb,                # The model we fit
                       type = "response",     # We'd like the predictions on the scale of the response
                       conf_level = 0.95,     # With a 95% confidence interval
                       condition = c("trmt_bin", "yr_s"), # Plot predictions for "l.illum_s" while holding all others at their mean
                       vcov = TRUE) +         # Compute and display the variance-covariance matrix
  xlab("") +            # Labels for axes
  ylab("bat calls") +
  theme_bw()       

t
pred2.3 <- predictions(m1.2nb, 
                       newdata = datagrid(yr_s = unique(bm2$yr_s),
                                          trmt_bin = unique(bm2$trmt_bin),
                                          sp=NA),
                       # re.form = NA,
                       conf_level = 0.95) 

pred2.3

# calculate non-standardized year for graphics
pred2.3 <- pred2.3 %>%
  mutate(yr = case_when(
    yr_s == -1 ~ 2021,
    yr_s == 0 ~ 2022,
    TRUE ~ 2023
  ))
# add dark and light 

pred2.3$treatment<-ifelse(pred2.3$trmt_bin==1, "light", "dark")

# Create the plot using ggplot2
p2.3 <- ggplot(pred2.3, aes(x = factor(yr), y = estimate, colour = factor(treatment))) +
  geom_point() +
  geom_line() +
  geom_smooth(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(x = "Year", y = "Predicted bat calls", color = "Treatment") +
  theme_bw()

p2.3

# marginal effect ear/arm ratio. 
pred2.4<- predictions(m1.2nb, 
                      by= "ear.arm_s",
                       newdata = datagrid(ear.arm_s = seq(min(bm2$ear.arm_s), max(bm2$ear.arm_s),
                                                            length.out = 100)),
                      re.form = NA) # re.form = NA to get the marginal effect of wind speed without random effects.)
#calculate non-standardized ear/arm ratio for graphics
pred2.4$ear.arm<-pred2.4$ear.arm_s*(2*sd(bm2$ear.arm))+mean(bm2$ear.arm) # recalculate the non-standardized ear/arm ratio for graphics.


p2.4<- ggplot(pred2.4, aes(x = ear.arm, y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(y = "Bat Calls", x = "Ear/Arm Ratio", title = "Marginal Effect of Ear/Arm Ratio") +
  theme_minimal() +
  scale_color_viridis(discrete = TRUE, option = "D") +  # Use a colorblind-friendly palette
  theme(legend.title = element_blank())  # Remove legend title for cleaner appearance

p2.4

# marginal effect for total lepidoptera.

# Remove NA values from the t.lepidoptera_s variable
cleaned_t_lepidoptera_s <- na.omit(bm2$t.lepidoptera_s)

# Create the sequence without NA values
t_lepidoptera_s_seq <- seq(min(cleaned_t_lepidoptera_s), max(cleaned_t_lepidoptera_s), length.out = 100)

# Generate predictions using the cleaned sequence
pred2.5 <- predictions(m1.2nb, 
                       by = "t.lepidoptera_s",
                       newdata = datagrid(t.lepidoptera_s = t_lepidoptera_s_seq),
                       re.form = NA) # re.form = NA to get the marginal effect without random effects

# Calculate non-standardized total lepidoptera for graphics

pred2.5$t.lepidoptera<-pred2.5$t.lepidoptera_s*(2*sd(bm2$t.lepidoptera,na.rm = T))+mean(bm2$t.lepidoptera, na.rm=T) # recalculate the non-standardized total lepidoptera for graphics.


p2.5 <- ggplot(pred2.5, aes(x = t.lepidoptera, y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(y = "Bat Calls", x = "Total Lepidoptera", title = "Marginal Effect of Total Lepidoptera") +
  theme_minimal() +
  scale_color_viridis(discrete = TRUE, option = "D") +  # Use a colorblind-friendly palette
  theme(legend.title = element_blank())  # Remove legend title for cleaner appearance
p2.5

# marginal effect for interaction treatment and yday

pred2.6 <- predictions(m1.2nb,
                       re.form = NA)
#adding treatment labels
pred2.6$tmt<-ifelse(pred2.6$trmt_bin==1, "light", "dark")
#calculating non-standardized jday for graphics
pred2.6$jday<-pred2.6$jday_s*(2*sd(bm2$jday))+mean(bm2$jday) # recalculate the non-standardized jday for graphics.


p2.6<-ggplot(pred2.6, aes(x = jday, y = estimate, colour =tmt )) +
  geom_point() +
  geom_smooth(method = 'auto', colour="black") +
  facet_wrap(~tmt) +
  labs(y = "Bat Calls", x = "Julian Day", title = "Marginal Effect of Julian") +
  theme_minimal() +
  scale_color_viridis(discrete = T, option = "H") +  # Use a colorblind-friendly palette
  theme(legend.title = element_blank())  # Remove legend title for cleaner appearance

p2.6
# now the random effects. 

# treatment by sp.

pred2.7 <- predictions(m1.2nb) # re.form = NULL to get the marginal effect without random effects

pred2.7 <- predictions(m1.2nb,
                     newdata = datagrid(sp = bm2$sp,
                                        trmt_bin = unique(bm2$trmt_bin)))

p2.7 <- ggplot(pred2.7, aes(x = trmt_bin, y = estimate, color = sp)) +
  geom_point() +
  geom_line() +
  labs(y = "Bat Calls", x = "Treatment", title = "Marginal Effect of Treatment by Species") +
  theme_minimal() +
  scale_color_viridis(discrete = TRUE, option = "D") +  # Use a colorblind-friendly palette
  theme(legend.title = element_blank())  # Remove legend title for cleaner appearance

p2.7
p2.7i<-p2.7+ facet_wrap(~sp, scales = "free_y")
p2.7i



# Create a list of plots and their corresponding filenames
figures_dir <- "figures/glmm_v1/marginleffects_out"


plots <- list(p1.1, p1.2, p2, p3, p4, p5.1, p5, p6)
filenames <- c("p1.1.png", "p1.2.png", "p2.png", "p3.png", "p4.png", "p5.1.png", "p5.png", "p6.png")

# Save all plots
lapply(seq_along(plots), function(i) {
  ggsave(filename = file.path(figures_dir, filenames[i]), plot = plots[[i]], device = "tiff", width = 15, height = 10, units = "cm")
})

library(sjPlot)
a1<-tab_model(m1.9nb, show.icc = TRUE, show.aic = TRUE, show.re.var = TRUE)
a1
plot_predictions(m1.9nb, condition = c("jday_s", "trmt_bin"), type = "response", vcov = T) +
  theme_minimal() +
  labs(y = "Predicted Bat Calls", x = "Julian Day", color = "Treatment")


# plots -------------------------------------------------------------------



# partial predictor manual jday -------------------------------------------



# jday --------------------------------------------------------------------


# values to use
n <- 100
int <- rep(1, n)

# obs. values
ord.day <- seq(min(bm2[, "jday"]), max(bm2[, "jday"]), length.out = n)
# extraer el jday sqr 
ord.day.sqr<-ord.day^2
#std pred
jday.s <- scale(ord.day)
jday.sqr<- scale(ord.day.sqr)
#extract fixed coef jday

#predicted ab
#
fday<-confint(m1.5nb)


predabund <- exp(t(fday[c(1,3,4),]) %*% t(cbind(int, jday.s, jday.sqr)))


#Data

abunddf <- data.frame(t(predabund), jday.s, ord.day) # make a df with all the above

colnames(abunddf)[1:3] <- c( "lowCI", "highCI", "Mean")

ggplot(abunddf, aes(x = ord.day, y = Mean)) +
  theme_classic(base_size = 17) +
  ylab("bat calls") +
  xlab("jday") +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.3, aes(ymin = lowCI, ymax = highCI))





# trmt_bin partial predictor ----------------------------------------------


# obs. values
trmt_bin<-c("dark", "light") # how do you plot the obs. values when you have
trmt_bin_s<- c(-1,1) # this should be -1 and 1  and remember to rerun the model for that 



# summary_m1.5nb <- summary(m1.5nb) # save the summary of the model
# feff<-summary_m1.5nb$coefficients$cond # get fixed eff
# c0<-feff["(Intercept)","Estimate"] # get intercept
# b2<-feff["trmt_bin", "Estimate"]   # get slope for treatment
cint<- confint(m1.5nb)[c(1,2 ), ]
cint

# mean.pred <- c(exp( trmt_bin_s[1]*f.trmt[1]), exp(trmt_bin_s[2]*(f.trmt[1]+f.trmt[2]) )) example 
# mean.pred <- c(exp( trmt_bin_s[1]*c0), exp(trmt_bin_s[2]*(c0+b2) ))

mean.pred <- c(exp(cint[1,3]+trmt_bin_s[1]*cint[2,3]),
               exp(cint[1,3]+trmt_bin_s[2]*cint[2,3]))


lowCI <-  c(exp( cint[1,1]+trmt_bin_s[1]*cint[2,1]), exp(cint[1,1]+trmt_bin_s[2]*cint[2,1] ))

highCI <- c(exp( cint[1,2]+trmt_bin_s[1]*cint[2,2]), exp(cint[1,2]+trmt_bin_s[2]*cint[2,2] ))
abunddf <- data.frame(mean.pred, lowCI, highCI, trmt_bin)

colnames(abunddf )[1:3] <- c( "Mean", "lowCI", "highCI" )

ggplot( abunddf, aes( x = trmt_bin, y = Mean) ) +  
  theme_classic( base_size = 17) +
  ylab( "bat calls" ) +
  xlab( "treatment" ) +
  geom_point()+
  geom_errorbar( aes( ymin = lowCI, ymax = highCI ) )




# lunar illumination partial predictor ------------------------------------

# obs. values
n<-100
int <- rep(1, n) #generate an interval of ones


l.illum<-seq(min(bm2$l.illum), max(bm2$l.illum),length.out = n) # get the observed values
# scale  l.illum
l.illum_s<- scale(l.illum)

# extract fixed coefficients


cint<-confint(m1.5nb)
cint[c(1,6),]

predcalls <- exp(t(cint[c(1,6),]) %*% t(cbind(int, l.illum_s)))



#Data

l.illumdf <- data.frame(t(predcalls), l.illum_s, l.illum) # make a df with all the above

colnames(l.illumdf)[1:3] <- c( "lowCI", "highCI", "Mean")

ggplot(l.illumdf, aes(x = l.illum, y = Mean)) +
  theme_classic(base_size = 12) +
  ylab("bat calls") +
  xlab("lunar illumination") +
  geom_line(size = .4) +
  geom_ribbon(alpha = 0.3, aes(ymin = lowCI, ymax = highCI))

# Improved plot
p1<-ggplot(l.illumdf, aes(x = l.illum, y = Mean)) +
  geom_line(size = .75, color = "black") + # Adjust line size and color
  geom_ribbon(aes(ymin = lowCI, ymax = highCI), fill = "grey", alpha = 0.5) + # Customize ribbon fill and transparency
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, ), # Center and bold title
    axis.title = element_text(), # Bold axis titles
    axis.text = element_text(color = "black"), # Ensure axis text is visible
    axis.line = element_line(color = "black"), # Darken axis lines for better visibility
    panel.grid = element_blank(), # Remove grid lines for a cleaner look
    plot.margin = margin(10, 10, 10, 10) # Add margin around the plot
  ) +
  labs(
    title = "Effect of Lunar Illumination on Bat Calls", # Add a title
    y = "Bat Calls", # More descriptive y-axis label
    x = "Lunar Illumination " # More descriptive x-axis label
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # Adjust y-axis to avoid clipping
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) # Adjust x-axis for balance

p1
ggsave(filename = "l.illum_feff_v1.png",plot = p1,device = "png", path = 'figures/glmm_v1/' )





# temperature -------------------------------------------------------------

n=100
int # interval of 1s
# observed values 

tmp<-seq(min(bm2$avg_temperature, na.rm = T ), max(bm2$avg_temperature, na.rm = T), length.out= n)
tmp_s<- scale(tmp)

#fixed effect

cint<- confint(m1.5nb)
cint[c(1,7),]

predcalls <- exp(t(cint[c(1,7),]) %*% t(cbind(int, tmp_s)))

tmpdf <- data.frame(t(predcalls), tmp, tmp_s) # make a df with all the above
colnames(tmpdf)[1:3] <- c( "lowCI", "highCI", "Mean")  #This is to label the data frame appropiately. 

p2<-ggplot(tmpdf, aes(x = tmp, y = Mean)) +
  geom_line(size = .75, color = "black") + # Adjust line size and color
  geom_ribbon(aes(ymin = lowCI, ymax = highCI), fill = "grey", alpha = 0.5) + # Customize ribbon fill and transparency
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, ), # Center and bold title
    axis.title = element_text(), # Bold axis titles
    axis.text = element_text(color = "black"), # Ensure axis text is visible
    axis.line = element_line(color = "black"), # Darken axis lines for better visibility
    panel.grid = element_blank(), # Remove grid lines for a cleaner look
    plot.margin = margin(10, 10, 10, 10) # Add margin around the plot
  ) +
  labs(
    title = "Effect of average night temperature on Bat Calls", # Add a title
    y = "Bat Calls", # More descriptive y-axis label
    x = "Temperature °C" # More descriptive x-axis label
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # Adjust y-axis to avoid clipping
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) # Adjust x-axis for balance
p2
ggsave(filename = "temp_feff_v1.png",plot = p2,device = "png", path = 'figures/glmm_v1/' )


# wind --------------------------------------------------------------------


n=100
int # interval of 1s

# observed values 

wnd<-seq(min(bm2$avg_wind_speed, na.rm = T ), max(bm2$avg_wind_speed, na.rm = T), length.out= n)
wnd_s<- scale(wnd)

#fixed effect

cint<- confint(m1.5nb)
cint[c(1,8),]

predcalls <- exp(t(cint[c(1,8),]) %*% t(cbind(int, wnd_s)))

wnddf <- data.frame(t(predcalls), wnd, wnd_s) # make a df with all the above
colnames(wnddf)[1:3] <- c( "lowCI", "highCI", "Mean")  #This is to label the data frame appropiately. 

p3<-ggplot(wnddf, aes(x = tmp, y = Mean)) +
  geom_line(size = .75, color = "black") + # Adjust line size and color
  geom_ribbon(aes(ymin = lowCI, ymax = highCI), fill = "grey", alpha = 0.5) + # Customize ribbon fill and transparency
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, ), # Center and bold title
    axis.title = element_text(), # Bold axis titles
    axis.text = element_text(color = "black"), # Ensure axis text is visible
    axis.line = element_line(color = "black"), # Darken axis lines for better visibility
    panel.grid = element_blank(), # Remove grid lines for a cleaner look
    plot.margin = margin(10, 10, 10, 10) # Add margin around the plot
  ) +
  labs(
    title = "Effect of average night wind on Bat Calls", # Add a title
    y = "Bat Calls", # More descriptive y-axis label
    x = "wind m/s" # More descriptive x-axis label
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # Adjust y-axis to avoid clipping
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) # Adjust x-axis for balance
p3
ggsave(filename = "wnd_feff_v1.png",plot = p3,device = "png", path = 'figures/glmm_v1/' )


# elevation ---------------------------------------------------------------
# obs. values
n<-100
int <- rep(1, n) #generate an interval of ones


elv<-seq(min(bm2$elev_mean), max(bm2$elev_mean),length.out = n) # get the observed values
# scale  
elv_s<- scale(elv)

# extract fixed coefficients
cint<-confint(m1.5nb)
cint[c(1,10),]

predcalls <- exp(t(cint[c(1,10),]) %*% t(cbind(int, elv_s)))

#Data

elvdf <- data.frame(t(predcalls), elv, elv_s) # make a df with all the above
head(elvdf)
colnames(elvdf)[1:3] <- c( "lowCI", "highCI", "Mean")

ggplot(elvdf, aes(x = elv, y = Mean)) +
  theme_classic(base_size = 12) +
  ylab("bat calls") +
  xlab("elevation") +
  geom_line(size = .4) +
  geom_ribbon(alpha = 0.3, aes(ymin = lowCI, ymax = highCI))

# Improved plot
p4<-ggplot(elvdf, aes(x = elv, y = Mean)) +
  geom_line(size = .75, color = "black") + # Adjust line size and color
  geom_ribbon(aes(ymin = lowCI, ymax = highCI), fill = "grey", alpha = 0.5) + # Customize ribbon fill and transparency
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, ), # Center and bold title
    axis.title = element_text(), # Bold axis titles
    axis.text = element_text(color = "black"), # Ensure axis text is visible
    axis.line = element_line(color = "black"), # Darken axis lines for better visibility
    panel.grid = element_blank(), # Remove grid lines for a cleaner look
    plot.margin = margin(10, 10, 10, 10) # Add margin around the plot
  ) +
  labs(
    title = "Effect of elevation on Bat Calls", # Add a title
    y = "Bat Calls", # More descriptive y-axis label
    x = "elevation m " # More descriptive x-axis label
  ) 
  # scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # Adjust y-axis to avoid clipping
  # scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) # Adjust x-axis for balance

p4
ggsave(filename = "elv_feff_v1.png",plot = p4,device = "png", path = 'figures/glmm_v1/' )


# year marginal effect ----------------------------------------------------


# obs. values
yr<- unique(bm2$yr)
yr
yr_s<- unique(bm2$yr_s)
yr_s

cint <- confint(m1.5nb)[c(1, 9), ]
cint 

mcalls <- c(exp(cint[1, 3] + yr_s[1] * cint[2, 3]),
            exp(cint[1, 3] + yr_s[2] * cint[2, 3]),
            exp(cint[1, 3] + yr_s[3] * cint[2, 3]))# after Jen's correction

lowCI <-  c(exp(cint[1, 1] + yr_s[1] * cint[2, 1]),
            exp(cint[1, 1] + yr_s[2] * cint[2, 1]),
            exp(cint[1, 1] + yr_s[3] * cint[2, 1]))

highCI <- c(exp(cint[1, 2] + yr_s[1] * cint[2, 2]),
            exp(cint[1, 2] + yr_s[2] * cint[2, 2]),
            exp(cint[1, 2] + yr_s[3] * cint[2, 2]))


abunddf <- data.frame(mcalls, lowCI, highCI, yr)

colnames(abunddf )[1:3] <- c( "Mean", "lowCI", "highCI" )

p5<-ggplot( abunddf, aes( x = yr, y = Mean) ) +  
  theme_classic( base_size = 12) +
  ylab( "bat calls" ) +
  xlab( "year" ) +
  geom_point()+
  geom_errorbar( aes( ymin = lowCI, ymax = highCI ) )+
# labels
  ylab("Bat Calls") +
  xlab("Year") +
  
  # Customize axis text and title
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  
  # Optionally, add a title and subtitle
  ggtitle("Distribution of Bat Calls by Year")

p5
ggsave(filename = "yr_feff_v1.png",plot = p5,device = "png", path = 'figures/glmm_v1/' )

# random effects plots ----------------------------------------------------


sl=100
ones = rep(1,100)

# obs. values
ord.day <- seq(min(bm2[, "jday"]), max(bm2[, "jday"]), length.out = n)
# extraer el jday sqr 
ord.day.sqr<-ord.day^2
#std pred
jday.s <- scale(ord.day)
jday.sqr<- scale(ord.day.sqr)

#extract fixed coef jday
#pull out random effects at the sp level #
ran.efs <- ranef( m1.5nb )$cond$sp

#pull out fixed effects
fix.efs <- fixef( m1.5nb )$cond
#view
fix.efs

cint<-confint(m1.5nb)

rss <- ran.efs

# we correct the random effects adding the fixes effects. 

rss[, 1] <- rss[, 1] + fix.efs[1]  #adding fixed effects to each of the random effects
rss[, 2] <- rss[, 2] + fix.efs[2]
rss[, 3] <- rss[, 3] + fix.efs[3]
rss[, 4] <- rss[, 4] + fix.efs[4]

#view
rss
a<-rss[,c(1,3,4)]
t(a)                      # why we have to transpose the tables?

b<-t( cbind( ones, jday.s, jday.sqr))
b
indpred<- exp( as.matrix(a) %*% as.matrix(b) )

abunddf <- data.frame(t(indpred), jday.s, ord.day)

ggplot(abunddf, aes(x = ord.day, y = MYOLUC    )) +
  theme_classic(base_size = 17) +
  ylab("bat calls Myluc") +
  xlab("jday") +
  geom_line(linewidth = 1.5) 


# Create the melted data for plotting# Create the melted data for plotting# Create the melted data for plotting all columns
abunddf_melted <- melt(abunddf, id.vars = "ord.day",variable.name = "sp", value.name ="predicted calls" , measure.vars = 1:14)
unique(abunddf_melted$sp)
t<-left_join(abunddf_melted, batnames, )

abunddf_long<- pivot_longer(abunddf, cols= 1:14, names_to = "sp", values_to = "predicted calls")
unique(abunddf_long$sp)
head(abunddf_long)

# Create the ggplot object
p6<-ggplot(abunddf_melted, aes(x = ord.day, y = `predicted calls`, color = sp)) +
  scale_color_viridis_d()+
  # Add geom_point to plot points for each variable
  scale_shape_manual(values = 1:14)+ 
  geom_point(size = 2) +
    # Use different shapes
  # Set labels and title
  labs(title = "Abundance by Day", x = "Day", y = "bat calls") +
  # Set theme for better visuals (optional)
  theme_classic()
p6

p6 <- ggplot(abunddf_melted, aes(x = ord.day, y = `predicted calls`, color = sp, shape = sp)) +
  scale_color_manual(values = rep("white", length(unique(abunddf_melted$sp)))) +  # Set all points to white
  scale_shape_manual(values = 1:14) +  # Custom shapes for each bat species
  geom_point(size = 2) +
  # Set labels and title
  labs(title = "", x = "Day", y = "Bat Calls") +
  # Set theme for better visuals
  theme_classic(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "black"),  # Set plot background to black
    panel.background = element_rect(fill = "black"),  # Set panel background to black
    text = element_text(color = "white"),  # Set all text to white
    axis.text = element_text(color = "white"),  # Set axis text to white
    axis.title = element_text(color = "white"),  # Set axis titles to white
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.background = element_rect(fill = "black"),  # Set legend background to black
    legend.key = element_rect(fill = "black"),  # Set legend keys to black
    legend.text = element_text(color = "white"),  # Set legend text to white
    legend.title = element_text(color = "white")  # Set legend title to white
  )


print(p6)

p6.1<-ggplot(abunddf_long, aes(x = ord.day, y = `predicted calls`, color = sp)) +
  # Add geom_point to plot points for each variable
  geom_point(size = 1) +
  facet_wrap( ~ sp, scales = "free_y") +
  # Set labels and title
  labs(title = "Bat calls by julian day", x = "julian day", y = "bat calls") +
  # Set theme for better visuals (optional)
  theme_classic()
rm(abunddf)
p6.1




p6
ggsave(filename = "jday_raneff_v2.png",plot = p6,device = "png", path = 'figures/glmm_v1/' )
ggsave(filename = "jday_raneff_v2.png",plot = p6.1,device = "png", path = 'figures/glmm_v1/' )


# treatment  random effects plot

trmt<- c("lit", "dark") 
trmt_bin_s <- c(-1,1)
trmt_bin_s

# correcting random eff
ran.efs <- ranef( m1.5nb )$cond$sp # get the random effects 
randint<- ran.efs[,1]
randslope<- ran.efs[,2]

cint<-confint(m1.5nb)[1:2,] # get fixed effects from the model
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

ggsave(filename = "trmt_raneff_v1.png",plot = p7,device = "png", path = 'figures/glmm_v1/' )
ggsave(filename = "trmt_raneff_v2.png",plot = p7.1,device = "png", path = 'figures/glmm_v1/' )



# save models -------------------------------------------------------------

#save image 
save.image(file = "working_env/glmm_v2.RData")

save(m1.2nb, file = "models/m1.2nb.RData") # best model up to 2/19/2025

load("models/my_models.RData")




# trash -------------------------------------------------------------------

# the plot below takes long to run and does not display the data well.
# # Adjusted code to create a violin plot
# ggplot(bm2, aes(x = factor(jday), y = n, fill = treatmt)) +
#   geom_violin(trim = FALSE) +  # Violin plot to show data distribution
#   geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +  # Adding jitter for better visibility of individual points
#   facet_wrap(~ sp, scales = "free_y") +  # Facet by species with free y scales
#   theme_minimal() +  # Use a minimal theme for better aesthetics
#   scale_fill_manual(values = c("#0033A0", "#D64309")) +  # Custom color for treatments
#   labs(title = "Bat Acoustic Activity 2021-2023",
#        x = "Julian Day",
#        y = "Number of Calls",
#        fill = "Treatment") +  # Label adjustments
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability





# # marginal effect for trmt_bin by sp 
# plot_predictions(m1.9nb,
#                  type = "response",
#                  conf_level = 0.95,
#                  condition = c("sp", "trmt_bin"),
#                  vcov = TRUE)


# m1.7
# we run this one with the model we neede 

m1.7nb <- glmmTMB(
  n ~ Sum_Distance_s + jday_s + I(jday_s ^ 2)  + l.illum_s +
    avg_wind_speed_s + avg_temperature_s + yr_s +
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s ^ 2)+ | sp),
  data = bm2,
  nbinom2(link = "log")
)

summary(m1.7nb)

# Calculate residual deviance and residual degrees of freedom
residual_deviance <- deviance(m1.7nb)
residual_df <- df.residual(m1.7nb)

# Calculate c-hat using residual deviance
c_hat_deviance <- residual_deviance / residual_df
print(c_hat_deviance)

AIC(m1.5nb,m1.7nb)


# model m1.5nb plots

emmeans::emmeans(m1.5nb, "percent_s",type="response")

plot(ggeffects::ggpredict(m1.5nb, "percent_s [all]"))
plot(ggeffects::ggpredict(m1.5nb, "l.illum_s [all]"))
plot(ggeffects::ggpredict(m1.5nb, "avg_wind_speed_s [all]"))
plot(ggeffects::ggpredict(m1.5nb, "avg_temperature_s [all]"))


# rstan_models ------------------------------------------------------------


# Define the formula
formula <- n ~ trmt_bin + jday_s + I(jday_s^2)  + percent_s + l.illum_s +
  avg_wind_speed_s + avg_temperature_s + (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp)

# Fit the negative binomial model using rstanarm
m1.4_off <- stan_glmer(
  formula,
  family = neg_binomial_2(),  # Negative binomial family with log-link
  data = bm2,
  chains = 4,                 # Number of Markov chains
  cores = 4                    # Number of cores to use (adjust as needed)
)


loo_m1.4_off <- loo(m1.4_off)
print(loo_m1.4_off)

plot_model(m1.4_off)


# plot model m1.5nb --------------------------------------------------------

plot_model(m1.5nb, type = "est", show.values = TRUE, value.offset = 0.3)

plot_model(m1.5nb, type = "est", show.values = TRUE, value.offset = 0.3,
           ci.lvl = 0.95, dot.size = 3, line.size = 1) +
  theme_minimal() +
  labs(title = "Bat Vocal Activity 2021-2023",
       x = "Coefficient Estimate",
       y = "") +
  scale_color_manual(values = c("orange", "purple")) +
  theme(axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        legend.position = "bottom")

ff_jday <- Effect(c("jday_s", "I(jday_s^2)"), m1.5nb) # this does not work!


jday_coefficients <- summary_m1.5nb$coefficients$cond[, c("Estimate", "Std. Error", "Pr(>|z|)")]

# Plot the effects
plot(eff_jday, rug = TRUE, main = "Partial Predictor Plot for jday_s and I(jday_s^2)")



# plotting partial predictors manually



# values to use
n <- 100
int <- rep(1, n)

# obs. values
jday_s_values <- seq(from = min(bm2$jday), to = max(bm2$jday), length.out = 100)



# extraer el jday sqr 
#std pred
jday.s <- scale(jday_s_values)
jday.sqr<- scale(jday_s_values^2)

#extract fixed coef jday

t<-coef(m1.5nb) # ?????????????????????????????????????? what do we multiply in the step below????
t$cond$sp

c_inf<-confint(m1.5nb)

fix_eff<- c_inf[c(1,3,4), ]
# # intm<- summary_m1.5nb$coefficients$cond[1,1]
# # jday_coef <- summary_m1.5nb$coefficients$cond[3,1 ]
# # jday.sqr<- summary_m1.5nb$coefficients$cond[4,1]
# 
# fix_eff<-c(intm, jday_coef, jday.sqr)
# #predicted ab

predabund <- exp(t(fix_eff) %*% t(cbind( int, jday.s, jday.sqr)))
t(predabund)



#Data

abunddf <- data.frame( t(predabund), jday_s_values) # make a df with all the above

colnames(abunddf)[1:3] <- c( "lowCI", "highCI", "Mean")

ggplot(abunddf, aes(x = jday_s_values, y = Mean)) +
  theme_classic(base_size = 17) +
  ylab("bat calls") +
  xlab("jday") +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.3, aes(ymin = lowCI, ymax = highCI))







# scale data old way

# bm2$elevm_s<- scale(bm2$elev_mean, center = T, scale = T) Use the loop instead. 
# bm2$ndvim_s<- scale(bm2$ndvi_mean)
# bm2$percent_s<-scale(bm2$percent)
# bm2$jday_s<-scale(bm2$jday)
# bm2$PeakFreq_s<-scale(bm2$PeakFreq)
# bm2$l.illum_s<-scale(bm2$l.illum)
# bm2$wind_s<scale(bm2$avg_wind_speed)


# # Combine the plots using patchwork
# plots <- Filter(Negate(is.null), plots)# removes all the null or NAs from the plot. We have somefor missing data. in avg temp, wind and peak freq
# combined_plot <- wrap_plots(plots)
# # Print the combined plot
# print(combined_plot)

# variables_to_plot <- c("elev_mean", "ndvi_mean", "percent", "jday", 
#                        "PeakFreq", "l.illum", "avg_wind_speed","avg_temperature","phase","fraction" )
# # Generate the plots
# plots <- lapply(variables_to_plot, function(var) {
#   if (is.numeric(bm2[[var]])) {
#     ggplot(bm2, aes_string(x = var)) +
#       geom_histogram(fill = "blue", color = "black") +
#       theme_minimal() +
#       ggtitle(paste("Histogram of", var))
#   } else {
#     message(paste("Skipping non-numeric variable:", var))
#     NULL
#   }
# })
# 
# 
# # m1 <- glmer(n ~ (1|site),
# data = bm2, 
# family = poisson)
# exp(m1@beta)# base line for bats at all sites. 
# 
# summary(m1)
# 
# 
# 
# m1.1<- glmer(n ~ trmt_bin * jday_s + percent_s +elev_mean_s+  (1|site),
#              data = bm2, 
#              family = poisson)
# 
# plot_model(m1.1)
# 
# exp(m1.1@beta)
# summary(m1.1)
# 
# plot(fitted(m1.1), residuals(m1.1), main = "Residuals vs Fitted", 
#      xlab = "Fitted values", ylab = "Residuals")
# abline(h = 0, col = "red")
# 
# qqnorm(residuals(m1.1), main = "Normal Q-Q Plot")
# qqline(residuals(m1.1), col = "red")
# 
# # Histogram of residuals
# hist(residuals(m1.1), breaks = 30, main = "Histogram of Residuals", 
#      xlab = "Residuals")
# 
# # Scale-Location Plot
# plot(fitted(m1.1), sqrt(abs(residuals(m1.2))), main = "Scale-Location Plot",
#      xlab = "Fitted values", ylab = "Square Root of |Residuals|")
# abline(h = 0, col = "red")
# 
# # Calculate VIF
# vif(m1.1)
# 
#trying to make a quasi poisson.
# 
# m1.2 <- glmer(
#   n ~ trmt_bin + jday_s + I(jday_s^2) + ndvi_mean_s + percent_s + PeakFreq_s + l.illum_s +
#     avg_wind_speed_s + avg_temperature_s + (1 | site) + (1 + trmt_bin + ndvi_mean_s | sp ),
#   data = bm2,
#   family = 'poisson',
# )
# 
# summary(m1.2)
# 
# exp(coef(m1.2))
# 
# plot_model(m1.2,)
# 
# model_convergence <- m1.2@optinfo$conv$opt
# print(model_convergence)
# 
# # Check if the Hessian matrix is positive definite
# is_positive_definite <- all(eigen(m1.2@optinfo$derivs$Hessian)$values > 0)
# print(is_positive_definite)
# 
# 
# quasi_table <- function(m1.2,ctab=coef(summary(m1.2))) {
#   phi <- sum(residuals(m1.2, type="pearson")^2)/df.residual(m1.2)
#   qctab <- within(as.data.frame(ctab),
#                   {   `Std. Error` <- `Std. Error`*sqrt(phi)
#                   `z value` <- Estimate/`Std. Error`
#                   `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
#                   })
#   print(phi)
#   return(qctab)
# }
# 
# printCoefmat(quasi_table(m1.2),digits=2)
# 
# 
# 
# 
# 
# 
# 
# m1.3 <- glmer(
#   n ~ trmt_bin + jday_s + I(jday_s ^ 2) + ndvi_mean_s + percent_s + PeakFreq_s + l.illum_s +
#     avg_wind_speed_s + avg_temperature_s + (1 |site) + (1 + trmt_bin + ndvi_mean_s + jday_s + I(jday_s ^ 2) |sp),
#   data = bm2,
#   family = 'poisson',
#   #  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
# )
# 
# summary(m1.3)
# 
# printCoefmat(quasi_table(m1.3),digits=2)
# 
# model_convergence <- m1.3@optinfo$conv$opt
# print(model_convergence)
# 
# 
# # calculate c-hat
# # Calculate residual deviance and residual degrees of freedom
# residual_deviance <- deviance(m1.3)
# residual_df <- df.residual(m1.3)
# 
# # Calculate c-hat using residual deviance
# c_hat_deviance <- residual_deviance / residual_df
# print(c_hat_deviance)



# m1.4_quasi <- glm(
#   n ~ trmt_bin + jday_s + I(jday_s^2) + ndvi_mean_s + percent_s + PeakFreq_s + l.illum_s +
#     avg_wind_speed_s + avg_temperature_s,
#   data = bm2,
#   family = quasipoisson()
# )
# summary(m1.4_quasi)
# 
# 
# # Calculate residual deviance and residual degrees of freedom
# residual_deviance <- deviance(m1.4_quasi)
# residual_df <- df.residual(m1.4_quasi)
# 
# # Calculate c-hat using residual deviance
# c_hat_deviance <- residual_deviance / residual_df
# print(c_hat_deviance)
# 
# model with offset effort  -------------------------------------------------------

# na_species_count <- bm2 %>%
#   group_by(sp, Species, Four.letter.species.code) %>%
#   filter(is.na(Ear_Forearm_Ratio)) %>%
#   summarise(na_count = n()) %>%
#   arrange(desc(na_count))


# m1.4_off <- glmer.nb(
#   n ~ trmt_bin + jday_s + I(jday_s ^ 2) + ndvi_mean_s + percent_s  + l.illum_s +
#     avg_wind_speed_s + avg_temperature_s + (1 | site) + (1 + trmt_bin + ndvi_mean_s + jday_s + I(jday_s^2) | sp),
#   offset = eff.hrs_s,
#   data = bm2,
# )
# 
# summary(m1.4_off)
# 
# # Calculate residual deviance and residual degrees of freedom
# residual_deviance <- deviance(m1.4_off)
# residual_df <- df.residual(m1.4_off)
# 
# # Calculate c-hat using residual deviance
# c_hat_deviance <- residual_deviance / residual_df
# print(c_hat_deviance)
# 
# model_convergence <- m1.4_nb@optinfo$conv$opt
# print(model_convergence)
# 
# 
# 
# #rstan
# library(rstanarm)
# 
# 
# 
# 
# 
# 
# # plots with SJplot to see what I need to get manually
# # model m1.4_nb
# sjPlot::plot_model(m1.5nb, type = "re")
# 
# plot_jday <- effect("jday_s", m1.5nb, xlevels = list(jday_s = seq(
#   min(bm2$jday_s), max(bm2$jday_s), length.out = 100
# )))
# 
# plot_trmt_bin <- effect("trmt_bin", m1.4_nb, xlevels = list(trmt_bin = seq(
#   min(bm2$trmt_bin), max(bm2$trmt_bin), length.out = 100
# )))
# 
# plot_ndvi <- effect("ndvi_mean_s", m1.4_nb, xlevels = list(ndvi_mean_s = seq(
#   min(bm2$ndvi_mean_s), max(bm2$ndvi_mean_s), length.out = 100
# )))
# 
# plot_percent_s <- effect("percent_s", m1.4_nb, xlevels = list(percent_s = seq(
#   min(bm2$percent_s), max(bm2$percent_s), length.out = 100
# )))
# 
# plot_PeakFreq_s  <- effect("PeakFreq_s ", m1.4_nb, xlevels = list(PeakFreq_s  = seq(
#   min(bm2$PeakFreq_s,na.rm = TRUE), max(bm2$PeakFreq_s,na.rm = TRUE ), length.out = 100
# )))
# 
# plot_l.illum_s  <- effect("l.illum_s ", m1.4_nb, xlevels = list(l.illum_s  = seq(
#   min(bm2$l.illum_s,na.rm = TRUE), max(bm2$l.illum_s,na.rm = TRUE ), length.out = 100
# )))
# 
# plot_avg_wind_speed_s   <- effect("avg_wind_speed_s  ", m1.4_nb, xlevels = list(avg_wind_speed_s   = seq(
#   min(bm2$avg_wind_speed_s ,na.rm = TRUE), max(bm2$avg_wind_speed_s ,na.rm = TRUE ), length.out = 100
# )))
# 
# plot_avg_temperature_s     <- effect("avg_temperature_s    ", m1.4_nb, xlevels = list(avg_temperature_s     = seq(
#   min(bm2$avg_temperature_s   ,na.rm = TRUE), max(bm2$avg_temperature_s   ,na.rm = TRUE ), length.out = 100
# )))
# 
# plot(plot_jday, multiline = TRUE)
# plot(plot_trmt_bin, multiline = TRUE)
# plot(plot_ndvi, multiline = TRUE)
# plot(plot_percent_s, multiline = T)
# plot(plot_PeakFreq_s, multiline = T)
# plot(plot_l.illum_s, multiline = T)
# plot(plot_avg_wind_speed_s, multiline = T)
# plot(plot_avg_temperature_s, multiline = T)
# 

