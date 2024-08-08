# ---------------------------
##
## Script name:  glmm_v1
##
## Purpose of script: Start running models for the 2021 data with glmm
##
## Author: Carlos Linares, Jen Cruz (collaborator)
##
## Date Created: 05/21/2024
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes: need to annotate what script produced each data set imputed. 
##   
##
## ---------------------------
## # inputs ------------------------------------------------------------------
#   data_for_analysis/prep_for_glmm/bm.csv
#   data_for_analysis/Bat_trait.csv
# outputs ----------------------

# 

# this should be a database ready to analyze with the glmm_v1 script. 



# libraries  --------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(lme4)
library(sjPlot)
library(ggeffects)
library(car)
library(glmmTMB)
library(corrplot)


#load environment 
#last worked 08/01/2024
load(file = "working_env/glmm_v1.RData")

# data --------------------------------------------------------------------

bm<-read.csv('data_for_analysis/prep_for_glmm/bm.csv', header = T) 
filtered_bm <- bm[!(bm$AUTO.ID. %in% c("Noise","NoID")), ] # if noise is not filtered there is more records for dark sites. 
colnames(filtered_bm)[2]<-"sp"

# bm<-read.csv('data_for_analysis/data_glmm/bat_counts.csv', header = T) # bat counts by jday for 2021
# bm$datefromyday<-as_date(bm$jday, origin= "2021-01-01")
# filtered_bm <- bm[bm$AUTO.ID. != "Noise", ] # if noise is not filtered there is more records for dark sites. 
# colnames(filtered_bm)[3]<-"sp" # change Auto.ID to sp

btrait<-read.csv('data_for_analysis/Bat_trait.csv', header = T)
btrait$Species<-toupper(btrait$Species) # makes all caps
batnames<-read.csv('data_for_analysis/Species_bats.csv')
colnames(batnames)[3]<-"sp"
colnames(btrait)[2]<-"sp"

elev<-read.csv('data_for_analysis/elev/elevation.csv', header = T)
elev<-elev %>% rename("site" = "name")
elev$site <-tolower(elev$site)
#Use gsub to replace 'viz' with 'vizc' in the 'site' column of df1
elev$site <- gsub("viz(\\d{2})", "vizc\\1", elev$site)


ndvi<-read.csv('data_for_analysis/NDVI/NDVI_of_rip2021.csv', header = T)
ndvi<- ndvi %>% dplyr::select(-c("ele", "time", "magvar", "geoidheigh", "dgpsid") ) # remove unncessary cols
ndvi<-ndvi %>% rename("site" = "name")
ndvi<-ndvi %>% rename("ndvi_mean" = "X_mean")
ndvi$site<-tolower(ndvi$site)
#Use gsub to replace 'viz' with 'vizc' in the 'site' column of df1
ndvi$site <- gsub("viz(\\d{2})", "vizc\\1", ndvi$site)

weather<-read.csv('data_for_analysis/weather/nigh_averages.csv', header = T) #load nightly averages
weather$date<-as_date(weather$date)


# moon clean-up
moon<-read.csv('data_for_analysis/moon_pred.csv')
moon$date<- as_date(moon$date)
moon.adj<-moon %>% mutate(
  phase = ifelse(above_horizon==FALSE,0,phase),
  fraction= ifelse(above_horizon==FALSE,0,fraction),
  l.illum= ifelse(above_horizon==FALSE,0,l.illum)
)

# effort_hrs<-read.csv('data_for_analysis/data_glmm/effort_hrs.csv')# we might not need effort because we are not doing an offset model. 


# merge -------------------------------------------------------------------



namestraits<-left_join(btrait, batnames, by= "sp")
colnames(namestraits)[2]<-"sp4"
traits<- namestraits %>% dplyr::select(c("sp4","Six.letter.species.code","Mass","Aspect","PeakFreq","Loading" ) )
colnames(traits)[2]<-"sp"
rm(namestraits) 

filtered_bm<- left_join(filtered_bm, traits, by="sp")

summary(filtered_bm)
# rows_with_na <- filtered_bm[rowSums(is.na(filtered_bm)) > 0, ]

bm2<-left_join(filtered_bm, elev, by="site" )
bm2 <- left_join(bm2, ndvi, by = "site") %>%
  select(-c( "time", "buff_area.x", "buff_area.y"))

colnames(bm2)[1]<-"date" # change name noche to date for merging

bm2$date<-lubridate::as_date(bm2$date)
bm2<- left_join(bm2, weather, by="date")


bm2<- left_join(bm2, moon.adj, by=c("date"))
# bm2<- left_join(bm2, effort_hrs, by=c("jday","site"))


# Identify rows with NA values
rows_with_na <- bm2[rowSums(is.na(bm2)) > 0, ] # some wind and temp have NAs because we are missing august 2021 weather data.

summary(bm2)

# correlation  ------------------------------------------------------------


# check for correlation 
numeric_cols<- sapply(bm2, is.numeric) # separate all the num col
cor1<-bm2[,numeric_cols] #keeps just the numeric
cor1<-cor1 %>% select(-c("X_min","X_max","altitude","parallacticAngle","angle","lat", "lon", "yr" ))

c1<- cor(cor1,use="pairwise.complete.obs")
corrplot(c1, order= 'AOE')

# make jday 

bm2$jday<-yday(bm2$date)



# List of variable names to be scaled
variables_to_scale <- c(
  "elev_mean",
  "ndvi_mean",
  "percent",
  "jday",
  "PeakFreq",
  "l.illum",
  "avg_wind_speed",
  "avg_temperature",
  "phase",
  "fraction"
)

# Loop over each variable, scale it, and assign it back to the data frame with a new name
for (var in variables_to_scale) {
  bm2[[paste0(var, "_s")]] <- scale(bm2[[var]], center = TRUE, scale = TRUE)
}

# 
# # Loop over each variable, scale it by dividing by two standard deviations, and assign it back to the data frame with a new name
# for (var in variables_to_scale) {
#   # Check if the column exists and is numeric
#   if (var %in% names(bm2) && is.numeric(bm2[[var]])) {
#     mean_val <- mean(bm2[[var]], na.rm = TRUE)
#     sd_val <- sd(bm2[[var]], na.rm = TRUE)
#     bm2[[paste0(var, "_s")]] <- (bm2[[var]] - mean_val) / (2 * sd_val)
#   } else {
#     warning(paste("Column", var, "is not numeric or does not exist in the data frame."))
#   }
# }

# make treatment -1 to 1

bm2$trmt_bin <- ifelse(bm2$trmt_bin == 1, 1, -1)

# make year between -1:1
bm2 <- bm2 %>%
  mutate(yr_s = case_when(
    yr == 2021 ~ -1,
    yr == 2022 ~ 0,
    TRUE ~ 1
  ))



# explore data ------------------------------------------------------------

summary(bm2)
# there are some bat 

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

# models -------------------------------------------------------------------


# glm no random effects ---------------------------------------------------


# the model below is severely over-disperse. 
m1.4 <- glmmTMB(
  n ~ trmt_bin + jday_s + I(jday_s ^ 2) + percent_s + l.illum_s + elev_mean_s+
    avg_wind_speed_s + avg_temperature_s+ yr_s,
  data = bm2,
  nbinom2(link = "log"))
          

summary(m1.4)

# calculate c-hat
# Calculate residual deviance and residual degrees of freedom
residual_deviance <- deviance(m1.4)
residual_df <- df.residual(m1.4)

# Calculate c-hat using residual deviance
c_hat_deviance <- residual_deviance / residual_df
print(c_hat_deviance)

# Calculate Pearson's chi-square statistic
pearson_chisq <- sum(residuals(m1.4, type = "pearson")^2)

# Calculate c-hat using Pearson's chi-square statistic
c_hat_pearson <- pearson_chisq / residual_df
print(c_hat_pearson)

plot_model(m1.4, vline.color = "grey", transform = "exp")

plot_model(m1.4, type = "pred", terms = c("trmt_bin[exp]"), ci.lvl = NA)
plot_residuals(m1.4)
get_model_data(m1.4)

#----#----#---------------------------------------------------------------#



# Negative Binomial -------------------------------------------------------

#this one needs a long time to run 
library(MASS)
m1.4_nb <- glmer.nb(
  n ~ trmt_bin + jday_s + I(jday_s^2) + ndvi_mean_s + percent_s + PeakFreq_s + l.illum_s +
    avg_wind_speed_s + avg_temperature_s + (1 | site) + (1 + trmt_bin + ndvi_mean_s + jday_s + I(jday_s^2) | sp),
  data = bm2
)
summary(m1.4_nb)


# Calculate residual deviance and residual degrees of freedom
residual_deviance <- deviance(m1.4_nb)
residual_df <- df.residual(m1.4_nb)

# Calculate c-hat using residual deviance
c_hat_deviance <- residual_deviance / residual_df
print(c_hat_deviance)

model_convergence <- m1.4_nb@optinfo$conv$opt
print(model_convergence)

#coefficients
mcoef<-coef(m1.4_nb)

?vif()



# m1.5nb ------------------------------------------------------------------
m1.5anb <- glmmTMB(
  n ~ trmt_bin + jday + I(jday_s ^ 2) + percent_s  + l.illum_s + 
    avg_wind_speed_s + avg_temperature_s + yr + elev_mean_s +
    (1 |site) + (1 | sp),
  data = bm2,
  nbinom2(link = "log"))
 
m1.5nb <- glmmTMB(
  n ~ trmt_bin + jday + I(jday_s ^ 2) + percent_s  + l.illum_s + 
    avg_wind_speed_s + avg_temperature_s + yr + elev_mean_s +
  (1 |site) + (1 + trmt_bin + ndvi_mean_s + jday_s + I(jday_s ^ 2) | sp),
  data = bm2,
  nbinom2(link = "log")
)


summary(m1.5nb)
confint(m1.5nb)

# Calculate residual deviance and residual degrees of freedom
residual_deviance <- deviance(m1.5nb)
residual_df <- df.residual(m1.5nb)

# Calculate c-hat using residual deviance
c_hat_deviance <- residual_deviance / residual_df
print(c_hat_deviance)

saveRDS(m1.5nb,"models/m1.5nb")

plot_model(m1.5nb)


# model with offset effort  -------------------------------------------------------


m1.4_off <- glmer.nb(
  n ~ trmt_bin + jday_s + I(jday_s ^ 2) + ndvi_mean_s + percent_s  + l.illum_s +
    avg_wind_speed_s + avg_temperature_s + (1 | site) + (1 + trmt_bin + ndvi_mean_s + jday_s + I(jday_s^2) | sp),
   offset = eff.hrs_s,
  data = bm2,
)

summary(m1.4_off)

# Calculate residual deviance and residual degrees of freedom
residual_deviance <- deviance(m1.4_off)
residual_df <- df.residual(m1.4_off)

# Calculate c-hat using residual deviance
c_hat_deviance <- residual_deviance / residual_df
print(c_hat_deviance)

model_convergence <- m1.4_nb@optinfo$conv$opt
print(model_convergence)



#rstan
library(rstanarm)

# rstan_models ------------------------------------------------------------


# Define the formula
formula <- n ~ trmt_bin + jday_s + I(jday_s^2) + ndvi_mean_s + percent_s + l.illum_s +
  avg_wind_speed_s + avg_temperature_s + (1 | site) + (1 + trmt_bin + ndvi_mean_s + jday_s + I(jday_s^2) | sp)

# Fit the negative binomial model using rstanarm
m1.4_off <- stan_glmer(
  formula,
  family = neg_binomial_2(),  # Negative binomial family with log-link
  # offset = log(scale(bm2$eff.hrs)),    # Offset term
  data = bm2,
  chains = 4,                 # Number of Markov chains
  cores = 4                    # Number of cores to use (adjust as needed)
)


loo_m1.4_off <- loo(m1.4_off)
print(loo_m1.4_off)

plot_model(m1.4_off)
# plots -------------------------------------------------------------------


# plots with SJplot to see what I need to get manally
# model m1.4_nb
sjPlot::plot_model(m1.5nb, type = "re")

plot_jday <- effect("jday_s", m1.4_nb, xlevels = list(jday_s = seq(
  min(bm2$jday_s), max(bm2$jday_s), length.out = 100
)))

plot_trmt_bin <- effect("trmt_bin", m1.4_nb, xlevels = list(trmt_bin = seq(
  min(bm2$trmt_bin), max(bm2$trmt_bin), length.out = 100
)))

plot_ndvi <- effect("ndvi_mean_s", m1.4_nb, xlevels = list(ndvi_mean_s = seq(
  min(bm2$ndvi_mean_s), max(bm2$ndvi_mean_s), length.out = 100
)))

plot_percent_s <- effect("percent_s", m1.4_nb, xlevels = list(percent_s = seq(
  min(bm2$percent_s), max(bm2$percent_s), length.out = 100
)))

plot_PeakFreq_s  <- effect("PeakFreq_s ", m1.4_nb, xlevels = list(PeakFreq_s  = seq(
  min(bm2$PeakFreq_s,na.rm = TRUE), max(bm2$PeakFreq_s,na.rm = TRUE ), length.out = 100
)))

plot_l.illum_s  <- effect("l.illum_s ", m1.4_nb, xlevels = list(l.illum_s  = seq(
  min(bm2$l.illum_s,na.rm = TRUE), max(bm2$l.illum_s,na.rm = TRUE ), length.out = 100
)))

plot_avg_wind_speed_s   <- effect("avg_wind_speed_s  ", m1.4_nb, xlevels = list(avg_wind_speed_s   = seq(
  min(bm2$avg_wind_speed_s ,na.rm = TRUE), max(bm2$avg_wind_speed_s ,na.rm = TRUE ), length.out = 100
)))

plot_avg_temperature_s     <- effect("avg_temperature_s    ", m1.4_nb, xlevels = list(avg_temperature_s     = seq(
  min(bm2$avg_temperature_s   ,na.rm = TRUE), max(bm2$avg_temperature_s   ,na.rm = TRUE ), length.out = 100
)))

plot(plot_jday, multiline = TRUE)
plot(plot_trmt_bin, multiline = TRUE)
plot(plot_ndvi, multiline = TRUE)
plot(plot_percent_s, multiline = T)
plot(plot_PeakFreq_s, multiline = T)
plot(plot_l.illum_s, multiline = T)
plot(plot_avg_wind_speed_s, multiline = T)
plot(plot_avg_temperature_s, multiline = T)





# Sjplots  ----------------------------------------------------------------
#Here we plot the negative binomial model using sjplot
# coeffsicients
b<-plot_model(m1.5nb, type =c("re"), show.p = T)
a<-plot_model(m1.5nb, type =c("est"), se=T, show.p = T)

# marginal effectrs
c<-plot_model(m1.5nb, type = "pred")
c <- plot_model(m1.5nb, type = "pred", terms = c("trmt_bin [all]", "jday_s [all]", "ndvi_mean_s [all]"))


#diagnostics 

d<-plot_model(m1.5nb, type = "diag")




# partial predictor manual jday -------------------------------------------


# plotting partial predictors manually

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
f.jday <- fixef(m1.4_nb)[c(1, 3,4)] # extract the intercept and the coeficient for jday

#predicted ab

predabund <- exp(f.jday %*% t(cbind(int, jday.s, jday.sqr)))

# mean abu

mabund <- apply(predabund, MARGIN = 2, FUN = mean)

#95% CI

CIabund <- apply(
  predabund,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.975)
)

#Data

abunddf <- data.frame(mabund, t(CIabund), jday.s, ord.day) # make a df with all the above

colnames(abunddf)[1:3] <- c("Mean", "lowCI", "highCI")

ggplot(abunddf, aes(x = ord.day, y = Mean)) +
  theme_classic(base_size = 17) +
  ylab("bat calls") +
  xlab("jday") +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.3, aes(ymin = lowCI, ymax = highCI))





# treatment partial prfedictor

# obs. values
trmt_bin<-c("dark", "light") # how do you plot the obs. values when you have
trmt_bin_s<- c(0,1) # this should be -1 and 1  and remember to rerun the model for that 
f.trmt <- fixef(m1.4_nb)[c(1,2 )]

t<-confint(m1.4_nb)

mean.pred <- c(exp( trmt_bin_s[1]*f.trmt[1]), exp(trmt_bin_s[2]*(f.trmt[1]+f.trmt[2]) ))
lowCI <-  c(exp( trmt_bin_s[1]*t[1,1]), exp(trmt_bin_s[2]*(t[1,1]+t[2,1]) ))
highCI <- c(exp( trmt_bin_s[1]*f.trmt[1]), exp(trmt_bin_s[2]*(f.trmt[1]+f.trmt[2]) ))
abunddf <- data.frame(mean.pred, lowCI, highCI, trmt_bin)

colnames(abunddf )[1:3] <- c(  "Mean", "lowCI", "highCI" )

ggplot( abunddf, aes( x = trmt_bin, y = Mean) ) +
  theme_classic( base_size = 17) +
  ylab( "bat calls" ) +
  xlab( "treatment" ) +
  geom_point()+
  geom_errorbar( aes( ymin = lowCI, ymax = highCI ) )




# marginal effects  -------------------------------------------------------

day<-seq(min(bm2[,"jday"]), max(bm2[,"jday"]), length.out=n) # how do you plot the obs. values when you have
day.s<- scale(day)
f.day <- fixef(m1.4_nb)[c(1,2 )]
predabund <- exp( f.trmt %*% t( cbind( int, trmt_bin_s) ) )
mabund <- apply( predabund, MARGIN = 2, FUN = mean )
CIabund <- apply( predabund, MARGIN = 2, FUN = quantile, 
                  probs = c(0.025, 0.975) )
abunddf <- data.frame(mabund, t(CIabund), trmt_bin_s, trmt_bin)

colnames(abunddf )[1:3] <- c(  "Mean", "lowCI", "highCI" )

ggplot( abunddf, aes( x = trmt_bin_s, y = Mean) ) +
  theme_classic( base_size = 17) +
  ylab( "bat calls" ) +
  xlab( "treatment" ) +
  geom_line( size = 1.5) +
  geom_ribbon( alpha = 0.3, aes( ymin = lowCI, ymax = highCI ) )




# random effects plot

re_site <- ranef(m1.4_nb)$site
re_sp <-ranef(m1.4_nb)$sp
# Plot random intercepts for site

fixefs<-fixef(m1.4_nb)

rss.light<- exp(re_sp[,1] + fixefs[1]+  re_sp[,2]+fixefs[2])
rss.nolight <- exp(re_sp[,1]+ fixefs[1])


plot(re_site$'(Intercept)', xlab = "Site", ylab = "Random Intercept", main = "Random Intercept Variability: Site")
abline(h = 0, col = "red", lty = 2)  # Add a reference line at zero






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
jday_s_values <- seq(from = min(bm2$jday_s), to = max(bm2$jday_s), length.out = 100)

new_data <- data.frame(
  jday_s = jday_s_values,
  jday_s2 = jday_s_values^2,
  # Set other variables to their mean or reference values
  trmt_bin = -1,  # or the most common value
  ndvi_mean_s = mean(bm2$ndvi_mean_s),
  percent_s = mean(bm2$percent_s),
  PeakFreq_s = 0 ,
  l.illum_s = mean(bm2$l.illum_s),
  avg_wind_speed_s = 0,
  avg_temperature_s = 0,
  site = factor("long04"),  # or an appropriate reference level
  sp = factor("NoID")  # or an appropriate reference level
)

predictions <- predict(m1.5nb, newdata = new_data, type = "response", se.fit = TRUE)
new_data$fit <- predictions$fit
new_data$se.fit <- predictions$se.fit

ggplot(new_data, aes(x = jday_s, y = fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = fit - 1.96 * se.fit, ymax = fit + 1.96 * se.fit), alpha = 0.2) +
  labs(title = "Partial Predictor Plot for jday_s and I(jday_s^2)",
       x = "Standardized Julian Day (jday_s)",
       y = "Predicted Response") +
  theme_minimal()



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





# partial predictor trmt_bin

trmt_levels<- unique(bm2$trmt_bin)


new_data_trmt <- data.frame(
  trmt_bin = trmt_levels,
  jday_s = mean(bm2$jday_s),
  jday_s2 = mean(bm2$jday_s)^2,
  ndvi_mean_s = mean(bm2$ndvi_mean_s),
  percent_s = mean(bm2$percent_s),
  PeakFreq_s = 0,
  l.illum_s = mean(bm2$l.illum_s),
  avg_wind_speed_s = 0,
  avg_temperature_s = 0,
  site = factor("long04"),  # Use an appropriate reference level or the most common level
  sp = factor("NoID")  # Use an appropriate reference level or the most common level
)

predictions_trmt <- predict(m1.5nb, newdata = new_data_trmt, type = "response", se.fit = TRUE)

# Combine predictions with the new data
new_data_trmt$fit <- predictions_trmt$fit
new_data_trmt$se.fit <- predictions_trmt$se.fit

# Plot the predictions

ggplot(new_data_trmt, aes(x = factor(trmt_bin), y = fit)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = fit - 1.96 * se.fit, ymax = fit + 1.96 * se.fit), width = 0.1) +
  labs(title = "Partial Predictor Plot for trmt_bin",
       x = "Treatment (trmt_bin)",
       y = "Predicted Response") +
  theme_minimal()


#chat gpt code
# Extract random effects for sp (intercept only)
re_sp <- coef(m1.5nb)$cond$sp

# Plot random intercepts for sp
plot(re_sp$'(Intercept)', xlab = "Species", ylab = "Random Intercept", main = "Random Intercept Variability: Species")
abline(h = 0, col = "red", lty = 2)  # Add a reference line at zero

# Extract random effects for sp (slopes)
re_sp_slopes <- coef(m1.5nb)$cond$sp

# Plot random slopes for trmt_bin, jday_s, I(jday_s^2), ndvi_mean_s
par(mfrow = c(2, 2))  # Set up a 2x2 plot layout

plot(re_sp_slopes$trmt_bin, main = "Random Slope: trmt_bin", xlab = "Species")
abline(h = 0, col = "red", lty = 2)

plot(re_sp_slopes$jday_s, main = "Random Slope: jday_s", xlab = "Species")
abline(h = 0, col = "red", lty = 2)

plot(re_sp_slopes$I(jday_s^2), main = "Random Slope: I(jday_s^2)", xlab = "Species")
abline(h = 0, col = "red", lty = 2)

plot(re_sp_slopes$ndvi_mean_s, main = "Random Slope: ndvi_mean_s", xlab = "Species")
abline(h = 0, col = "red", lty = 2)

par(mfrow = c(1, 1))  # Reset plotting layout to default



# extracting marginal effects 
## Generate marginal effects for each species
library(ggeffects)
marginal_effects <- ggpredict(m1.5nb, terms = c("trmt_bin", "sp")) # this did not worked...

ggplot(marginal_effects, aes(x = x, y = predicted, color = group)) +
  geom_line() +
  facet_wrap(~ group, scales = "free") +  # Facet by species
  labs(title = "Marginal Effects by Species",
       x = "Predictor Variable",
       y = "Predicted Value") +
  theme_minimal()

plot_model(m1.5nb, type = "pred", terms = c("trmt_bin [all]", "sp"))

mydf<-predict_response(m1.5nb, terms = "trmt_bin[all]")
ggplot(mydf, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)

# save models -------------------------------------------------------------

#save image 
save.image(file = "working_env/glmm_v1.RData")

save(c(m1.2,m1.3,m1.4,m1.4_nb,m1.4_quasi, m1.5nb), file = 'models/glmm_v1.RData')
load("models/glmm_v1.RData")





# trash -------------------------------------------------------------------

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
