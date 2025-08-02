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
  "data.table",
  "janitor" # to clean names
)

#load environment 
load(file = "working_env/glmm_v2.RData") 

#load data ---------------------------------------------------------------

bm <- read_csv('data_for_analysis/prep_for_glmm/bm.csv') %>% 
  clean_names() # load bat data.

filtered_bm <- bm %>%
  filter(!auto_id %in% c("Noise", "NoID")) %>% # remove calls mark as Noise and NoID
  rename(sp = auto_id) # safe rename

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

# light

light<-read_csv("data_for_analysis/lights/lightspectra_pioneer.csv") %>% 
  clean_names()# load light data.

# filter col V/H and keep all horizontal 

light<-light %>%
  filter(vert_horiz == "Horizontal") # keep only horizontal light

# calculate the mean light for m1-m3 
light$mwatts <- rowMeans(light[, c("watts_m1", "watts_m2", "watts_m3")], na.rm = TRUE) # calculate the mean of the three columns mwat_1, mwat_2, mwat_3

# keep necessary columns
light <- light %>%
  select(c("site", "lux", "yr", "mwatts")) # keep only the columns we need

# merge -------------------------------------------------------------------

# activity index with bat data. 
bm2 <- left_join(filtered_bm, filtered_bm.ai[, c("sp", "site", "noche", "activity_min")], by = c("sp", "site", "noche"))
nrow(bm2)

# merge with traits
bm2<- left_join(bm2, select(btrait, six_sp, ear.arm), by = c("sp" = "six_sp"))

summary(bm2)
nrow(bm2) # 17772 rows
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

#make 0 the NAs for t.insect, and t. lepidoptera given we don't have that data.

bm2<- bm2 %>%
  mutate(
    t.insect = ifelse(t.insect==0, NA, t.insect),
    t.lepidoptera = ifelse(t.lepidoptera==0, NA, t.lepidoptera)
  )

# we have several weeks in 2023 where there is no data 5659 lines.

# merge light 

bm2<- bm2 %>% 
  left_join(light, by = c("site", "yr")) # merge with light data

summary(bm2) #this has the predictors and the y's 
nrow(bm2) # 17772 rows
# correlation  ------------------------------------------------------------
# before modelling we have to check for correlation between the predictors as VIF is not adequate for negative binomial models.
# 1. Select numeric columns, optionally drop unique id/group columns
numeric_data <- bm2 %>% 
  select(where(is.numeric)) %>% 
  select(-any_of(c("yr", "trmt_bin"))) # add/remove columns as needed

# 2. Remove zero-variance columns
numeric_data <- numeric_data %>% select(where(~sd(., na.rm = TRUE) > 0))

# 3. Compute correlation matrix
cor_mat <- cor(numeric_data, use = "pairwise.complete.obs")

# 4. Visualize correlation matrix
corrplot(cor_mat, order = 'AOE')
# from the plot it seems like the tmp and wind speed, twilight, total illumination, moon phase and moonlight are correlated not to be included in the model at the same time. Insect variables like total lepidoptera and and total insects are also correlated with moon phase and light but less than 0.4. Spectral measurments are corre

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
  "t.lepidoptera",
  "lux",
  "mwatts"
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
    nit_avg_wspm.s_s + yr_s  + t.lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp) +
    #interactions
    jday_s * trmt_bin + I(jday_s^2) * trmt_bin + yr_s * trmt_bin + trmt_bin*t.lepidoptera_s + trmt_bin*moonlight_s,
  data = bm2,
  family = nbinom2(link = "log")
)
summary(m1.4nb)
check_zeroinflation(m1.4nb)
check_overdispersion(m1.4nb)
r2(m1.4nb)

anova(m1.4nb) 

# after meeting with Jesse we discuss running a model without the ear/arm ratio and try to run a model with the insects and the try to use all insects but with just the two years. to see what is the pattern observed. 

#select just years 2021-2022 from the bm2 data.
bm2.y21.22<-bm2[bm2$yr != 2023,]

m1.5nb<- glmmTMB(
  #fixed effects
  n ~ trmt_bin + jday_s + I(jday_s^2) + moonlight_s +
    nit_avg_wspm.s_s + yr_s + t.insect_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) + moonlight_s | sp) +
    #interactions
    jday_s * trmt_bin + I(jday_s^2) * trmt_bin + yr_s * trmt_bin,
  data = bm2.y21.22,
  family = nbinom2(link = "log")
)

summary(m1.5nb)
performance(m1.5nb)
r2(m1.5nb)
check_model(m1.5nb)
check_overdispersion(m1.5nb)
check_collinearity(m1.5nb)

per1<-compare_performance(m1.5nb, m1.4nb, m1.3nb, m1.2nb, m1.1nb)

plot_model(m1.5nb)


# model with light spectra instead of binomial 

m1.6nb <- glmmTMB(
  #fixed effects
  n ~ lux_s + jday_s + I(jday_s^2) + moonlight_s +
    nit_avg_wspm.s_s + yr_s  + t.lepidoptera_s +
    #random effects
    (1 | site) + (1 + lux_s + jday_s + I(jday_s^2) | sp) +
    #interactions
    jday_s * lux_s + I(jday_s^2) * lux_s + yr_s * lux_s + lux_s*t.lepidoptera_s + lux_s*moonlight_s,
  data = bm2,
  family = nbinom2(link = "log")
)
summary(m1.6nb)
check_zeroinflation(m1.6nb)
check_overdispersion(m1.6nb)
check_singularity(m1.6nb)
r2(m1.6nb)
DHARMa::simulateResiduals(m1.6nb, plot = TRUE, quantreg = TRUE)
anova(m1.6nb) # the model with light spectra is not better than the model with the insects and ear/arm ratio.) 


# acoustic activity index models ------------------------------------

# now instead of using the bat call counts we are going to use the activity index from the bm2 dataset (activity_min).

m2.1<-glmmTMB(
  #fixed effects
  activity_min ~ trmt_bin + jday_s + I(jday_s^2) + moonlight_s +
    nit_avg_wspm.s_s + yr_s + t.lepidoptera_s +
    #random effects
    (1 | site) + (1 + trmt_bin + jday_s + I(jday_s^2) | sp) +
    #interactions
    jday_s * trmt_bin + I(jday_s^2) * trmt_bin + yr_s * trmt_bin,
  data = bm2,
  family = nbinom2(link = "log")
)

summary(m2.1)
plot_model(m2.1)
r2(m2.1)

# Marginal effect plots ---------------------------------------------------


# marginal with marginal effects package. 

p1.1<-plot_predictions(m1.4nb,                # The model we fit
                 type = "response",     # We'd like the predictions on the scale of the response
                 conf_level = 0.95,     # With a 95% confidence interval
                 condition = c("trmt_bin", "yr_s"), # Plot predictions for "l.illum_s" while holding all others at their mean
                 vcov = TRUE) +         # Compute and display the variance-covariance matrix
  xlab("") +            # Labels for axes
  ylab("bat calls") +
  theme_bw()                             # White background

p1.1 <- plot_predictions(m1.4nb,                # The model we fit
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

p1.2<-plot_predictions(m1.4nb,                # The model we fit
                 type = "response",     # We'd like the predictions on the scale of the response
                 conf_level = 0.95,     # With a 95% confidence interval
                 condition = c("jday_s","trmt_bin"), 
                 vcov = TRUE) +         # Compute and display the variance-covariance matrix
  xlab("") +            # Labels for axes
  ylab("bat calls") +
  theme_bw()        


# marginal effect for jday by sp
pred2 <- predictions(m1.4nb,
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


# marginal effect for year by sp

pred3 <- predictions(m1.4nb,
                     newdata = datagrid(sp = bm2$sp,
                                        yr_s = unique(bm2$yr_s)))

p3<- ggplot(pred3, aes(x =yr_s, y= estimate))+
              geom_point()+
  geom_errorbar( aes( ymin = conf.low, ymax = conf.high ) )+ 
  facet_wrap(~sp, scales="free")
  
            
p3

pred4<- predictions(m1.4nb,
                     newdata = datagrid(sp = bm2$sp,
                                        trmt_bin = unique(bm2$trmt_bin)))
# add species names with actual species. 
pred4<-left_join(pred4, btrait, by = c("sp" = "six_sp"))
# add light and dark treatment
pred4$tmt<-ifelse(pred4$trmt_bin==1, "light", "dark")


p4 <- ggplot(pred4, aes(x = trmt_bin, y = estimate, col=sname)) +
  geom_point() +
  geom_line(aes(group = sname)) 
p4

p4 + facet_wrap( ~ sp, )

preds<- predictions(m1.4nb)
preds$tmt<-ifelse(preds$trmt_bin==1, "light", "dark")
# now add the years -1=2021, 0=2022, 1=2023
preds$yr_s<-ifelse(preds$yr_s==-1, "2021", ifelse(preds$yr_s==0, "2022", "2023"))
preds<-left_join(preds, btrait, by = c("sp" = "six_sp"))


p5.1<- ggplot(preds, aes(x = jday_s, y = estimate, col=tmt))+ # not sure what is happening here. 
  geom_point()+
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


p5 <- ggplot(preds, aes(x = tmt, y = estimate, color = yr_s)) +
  geom_violin(width = 0.8, alpha = 0.4, position = position_dodge(width = 0.75)) +
  geom_point(position = position_dodge(width = 0.75), size = 1.5, alpha = 0.9) +
  facet_wrap(~sp, scales = "free_y") +
  labs(
    y = "Estimated Bat Calls",
    x = "",
    title = "Estimated Bat Calls by Treatment, Year, and Species",
    color = "Year"
  ) +
  scale_color_brewer(palette = "Dark2") +  # Choose a harmonious, colorblind-friendly palette
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major.x = element_blank()
  )
p5


# print(predictions(model = m1.9nb, 
#                   newdata = datagrid( yr_s= c(-1, 0, 1),
#                                       sp= unique(bm2$sp)),
#                   conf_level = 0.95))

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



# write the preds object to be able to map it on R. 
write.csv(preds, "data_for_analysis/glmm_v2/preds.csv")






# Marginal effect plots m1.4nb --------------------------------------------

# Below is marginal effect plot for moon illumination. 

p2.1<-plot_predictions(m1.4nb,                # The model we fit
                       type = "response",     # We'd like the predictions on the scale of the response
                       conf_level = 0.95,     # With a 95% confidence interval
                       condition = c("moonlight_s"), # Plot predictions 
                       vcov = TRUE) +         # Compute and display the variance-covariance matrix
  xlab("") +            # Labels for axes
  ylab("bat calls") +
  theme_bw()                             # White background
p2.1

# Create a data grid for predictions to replicate plot above.
pred2.1 <- predictions(m1.4nb, 
                       newdata = datagrid(moonlight_s = seq(min(bm2$moonlight_s), max(bm2$moonlight_s), 
                                                            length.out = 100)))
# Calculate non-standardized moonlight for graphics moonlight_s x (2xsd)+mean
sdmoon<-sd(bm2$moonlight)
mmoon<-mean(bm2$moonlight)
pred2.1$moonlight<-pred2.1$moonlight_s*(2*sdmoon)+mmoon


p2.1<-ggplot(pred2.1, aes(x = moonlight, y = estimate)) +
     geom_line() +
     geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
     labs(y = "Predicted Bat Calls", x = "Moonlight", title = "Marginal Effect of Moonlight") +
     theme_minimal() + 
     scale_color_viridis(discrete = TRUE, option = "D") +  # Use a colorblind-friendly palette
     theme(legend.title = element_blank())  # Remove legend title for cleaner appearance
p2.1

# Marginal effect plot for wind speed.
pred2.2<- predictions(m1.4nb, 
                       newdata = datagrid(nit_avg_wspm.s_s = seq(min(bm2$nit_avg_wspm.s_s), max(bm2$nit_avg_wspm.s_s), 
                                                            length.out = 100)),
                      re.form = NA) # re.form = NA to get the marginal effect of wind speed without random effects.)

# Calculate non-standardized wind speed for graphics
pred2.2$wspm.s<-pred2.2$nit_avg_wspm.s_s*(2*sd(bm2$nit_avg_wspm.s))+mean(bm2$nit_avg_wspm.s) # recalculate the non-standardized wind speed for graphics.

p2.2 <- ggplot(pred2.2, aes(x = wspm.s, y = estimate, ymin= conf.low, ymax=conf.high)) +
        geom_line()+
        geom_ribbon( alpha = 0.2)+
        labs(y = "Bat Calls", x = "Wind Speed", title = "Marginal Effect of Wind Speed m1.2nb") +
        theme_minimal() +
        scale_color_viridis(discrete = TRUE, option = "D") +  # Use a colorblind-friendly palette
        theme(legend.title = element_blank())  # Remove legend title for cleaner appearance
p2.2


# Marginal effect plot for year. 

t<-plot_predictions(m1.4nb,                # The model we fit
                       type = "response",     # We'd like the predictions on the scale of the response
                       conf_level = 0.95,     # With a 95% confidence interval
                       condition = c("trmt_bin", "yr_s"), # Plot predictions for "l.illum_s" while holding all others at their mean
                       vcov = TRUE) +         # Compute and display the variance-covariance matrix
  xlab("") +            # Labels for axes
  ylab("bat calls") +
  theme_bw()       

t
pred2.3 <- predictions(m1.4nb, 
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
  labs(y = "Bat Calls", x = "Ear/Arm Ratio", title = "Marginal Effect of Ear/Arm Ratio m1.2nb") +
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
pred2.5 <- predictions(m1.4nb, 
                       by = "t.lepidoptera_s",
                       newdata = datagrid(t.lepidoptera_s = t_lepidoptera_s_seq),
                       re.form = NA) # re.form = NA to get the marginal effect without random effects

# Calculate non-standardized total lepidoptera for graphics

pred2.5$t.lepidoptera<-pred2.5$t.lepidoptera_s*(2*sd(bm2$t.lepidoptera,na.rm = T))+mean(bm2$t.lepidoptera, na.rm=T) # recalculate the non-standardized total lepidoptera for graphics.


p2.5 <- ggplot(pred2.5, aes(x = t.lepidoptera, y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(y = "Bat Calls", x = "Total Lepidoptera", title = "Marginal Effect of Total Lepidoptera m1.2nb") +
  theme_minimal() +
  scale_color_viridis(discrete = TRUE, option = "D") +  # Use a colorblind-friendly palette
  theme(legend.title = element_blank())  # Remove legend title for cleaner appearance
p2.5

# marginal effect for interaction treatment and yday

pred2.6 <- predictions(m1.4nb,
                       re.form = NA)
#adding treatment labels
pred2.6$tmt<-ifelse(pred2.6$trmt_bin==1, "light", "dark")
#calculating non-standardized jday for graphics
pred2.6$jday<-pred2.6$jday_s*(2*sd(bm2$jday))+mean(bm2$jday) # recalculate the non-standardized jday for graphics.


p2.6<-ggplot(pred2.6, aes(x = jday, y = estimate, colour =tmt )) +
  geom_point() +
  geom_smooth(method = 'lm', colour="black") +
  facet_wrap(~tmt) +
  labs(y = "Bat Calls", x = "Julian Day", title = "Marginal Effect of Julian m1.2nb") +
  theme_minimal() +
  scale_color_viridis(discrete = T, option = "H") +  # Use a colorblind-friendly palette
  theme(legend.title = element_blank())  # Remove legend title for cleaner appearance

p2.6
# now the random effects. 

# treatment by sp.

pred2.7 <- predictions(m1.4nb) # re.form = NULL to get the marginal effect without random effects

pred2.7 <- predictions(m1.4nb,
                     newdata = datagrid(sp = bm2$sp,
                                        trmt_bin = unique(bm2$trmt_bin)))
# Calculate non-standardized treatment for graphics
pred2.7$trmt<-ifelse(pred2.7$trmt_bin==1, "light", "dark")
# add sp names with actual species.
pred2.7<-left_join(pred2.7, btrait, by = c("sp" = "six_sp"))


p2.7 <- ggplot(pred2.7, aes(x = trmt, y = estimate, color = cname, group = sp)) +
  geom_point() +
  geom_line() +
  labs(y = "Bat Calls", x = "Treatment", title = "Marginal Effect of Treatment by Species") +
  theme_minimal() +
  scale_color_viridis(discrete = TRUE, option = "D") +  # Use a colorblind-friendly palette
  theme(legend.title = element_blank())  # Remove legend title for cleaner appearance

p2.7
p2.7i<-p2.7+ facet_wrap(~cname)
p2.7i


# marginal effect for treatment by year




# Create a list of plots and their corresponding filenames
figures_dir <- "figures/glmm_v1/marginleffects_out"


plots <- list(p2.1, p2.2, p2.3, p2.4, p2.5, p2.6, p2.7, p2.7i)
filenames <- c("p2.1.png", "p2.2.png", "p2.3.png", "p3.png", "p4.png", "p5.1.png", "p5.png", "p6.png")

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


# Marginal effect plots m1.5nb --------------------------------------------

pred1.5<- predictions(m1.5nb, 
                       newdata = datagrid(sp = bm2$sp,
                                          trmt_bin = unique(bm2$trmt_bin)),
                       conf_level = 0.95)

# Calculate non-standardized treatment for graphics
pred1.5$trmt<-ifelse(pred1.5$trmt_bin==1, "light", "dark")

# create the plot for treatment by species
p1.5.1 <- ggplot(pred1.5, aes(x = trmt, y = estimate, color = sp, group = sp)) +
  geom_point() +
  geom_line() +
  labs(y = "Bat Calls", x = "Treatment", title = "Marginal Effect of Treatment by Species m1.5nb") +
  theme_minimal() +
  scale_color_viridis(discrete = TRUE, option = "D") +  # Use a colorblind-friendly palette
  theme(legend.title = element_blank())  # Remove legend title for cleaner appearance

p1.5.1

# Marginal effect plots m1.6nb





# marginal effect plots m1.6nb --------------------------------------------

plot_model(m1.6nb)
plot_model(m1.6nb, type = "re")


# community response to treatment. 

pred1 <- predictions(m1.6nb,
                     newdata = datagrid(sp  = NA, 
                                        lux_s = seq( min(bm2$lux_s), max(bm2$lux_s), length.out = 100),
                                        site = unique(bm2$site)
                                        ),
                     re.form = NA)

# Calculate non-standardized lux for graphics
mean_lux <- mean(bm2$lux, na.rm = TRUE)

sd_lux   <- sd(bm2$lux, na.rm = TRUE)

# add to the pred1 data frame
pred1 <- pred1 %>%
  mutate(lux = (lux_s * sd_lux) + mean_lux)

# ad lit vs dark sites to plot two lines

lit<-c("long01", "long03", "iron05", "iron01", "iron03")
pred1 <- pred1 %>%
  mutate(trmt = ifelse(site %in% lit, "lit", "dark")) # add dark and lit labels to sites

# add treatment. 

p1<-ggplot(pred1, aes(x = lux, y = estimate)) +
  geom_line()+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(
    x = "Light (lux)",
    y = "Predicted Bat Calls") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),   # remove major grid lines
    panel.grid.minor = element_blank())
p1


# species response to treatment. 

pred2 <- predictions(m1.6nb,,
                     newdata = datagrid(
                     sp = (bm2$sp),
                     lux_s = seq( min(bm2$lux_s), max(bm2$lux_s), length.out = 100)
                     )
)

pred2$sp <- factor(pred2$species,  # order species 
                        levels = pred2 %>% 
                          group_by(species) %>% 
                          summarise(mean_est = mean(estimate)) %>%
                          arrange(desc(mean_est)) %>% 
                          pull(species))

#simplify species names so the labels are just the initial of genus and the species like C.brachyrhynchos
pred2 <- pred2 %>%
  mutate(species_short = str_replace(species, 
                                     "^([A-Za-z])[a-z]+\\s+", 
                                     "\\1. "))

#as data set 

pred2<-as_tibble(pred2)

p2<-ggplot(pred2, aes(x = lux_s, y = estimate )) +
  geom_line(alpha = 0.5, color = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  scale_fill_manual(values = c("black", "white")) +
  facet_wrap(~ sp, scales = "free_y") +
  theme_minimal()
p2

ggplot(pred2, aes(x = species_short, y = estimate, color = factor(trmt_bin))) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = "Species", y = "Predicted Bird Calls", color = "Treatment") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(pred2, aes(x = factor(trmt_bin), y = estimate, group = species)) +
  geom_point(aes(color = factor(trmt_bin)), size = 3) +
  geom_line(aes(group = species)) +
  facet_wrap(~ species, scales= "free_y") +
  labs(x = "Treatment", y = "Predicted Bird Calls", color = "Treatment") +
  theme_minimal()



p2<-ggplot(pred2, aes(x = species_short, y = estimate, 
                      color = factor(trmt_bin), 
                      shape = factor(trmt_bin))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("darkblue", "orange"), 
                     labels = c("Dark", "Lit")) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("Dark", "Lit")) +
  labs(x = "Species", y = "Predicted Bird Calls", 
       color = "Treatment", shape = "Treatment") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

p2
ggsave(p2, file = "images/glmm_birdnet_v2/species_treatment_effects.png", width = 10, height = 6)

# summary to describe the results. 
species_means <- pred2 %>%
  group_by(species) %>%
  summarise(
    mean_dark = estimate[trmt_bin == -1],
    mean_lit  = estimate[trmt_bin == 1]
  ) %>%
  arrange(desc(mean_dark))


diff_data <- pred2 %>%
  select(species, trmt_bin, estimate) %>%
  pivot_wider(names_from = trmt_bin, values_from = estimate) %>%
  mutate(Difference = `1` - `-1`)

ggplot(diff_data, aes(x = species, y = Difference)) +
  geom_col(fill = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Species", y = "Treatment Effect (Light - Dark)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# effect sizes plots -------------------------------------------------------------------

plot_model(m1.5nb, type = "est", show.values = TRUE, value.offset = 0.3)

plot_model(m1.5nb, type = "re", terms = "trmt_bin", show.values = TRUE, value.offset = 0.3)

# save models -------------------------------------------------------------

#save image 
save.image(file = "working_env/glmm_v2.RData")

save(m1.5nb, file = "models/m1.5nb.RData") # best model up to 8/1/2025

load("models/my_models.RData")




# trash -------------------------------------------------------------------


# the function below is not working. 
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


p2 <- plot_marginal_effects(m1.1nb, "jday_s", "sp", "Marginal effect of jday_s by species")




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

