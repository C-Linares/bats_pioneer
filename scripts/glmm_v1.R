# ---------------------------
##
## Script name:  glmm_v1
##
## Purpose of script: Start running models for the 2021 data with glmm
##
## Author: Carlos Linares
##
## Date Created: 05/21/2024
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes: need to anotate what script produced each data set inputed. 
##   
##
## ---------------------------
## 


# libraries  --------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(lme4)
library(sjPlot)
library(ggeffects)
library(car)
library(glmmTMB)


# data --------------------------------------------------------------------


bm<-read.csv('data_for_analysis/data_glmm/bat_counts.csv', header = T) # bat counts by jday
bm$datefromyday<-as_date(bm$jday, origin= "2021-01-01")
filtered_bm <- bm[bm$AUTO.ID. != "Noise", ] # if noise is not filtered there is more records for dark sites. 
colnames(filtered_bm)[3]<-"sp" # change Auto.ID to sp

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

weather<-read.csv('data_for_analysis/weather/dailyavg.csv', header = T)
weather$date<-as_date(weather$date)

# moon clean-up
moon<-read.csv('data_for_analysis/moon_pred.csv')
moon$date<- as_date(moon$date)
moon.adj<-moon %>% mutate(
  phase = ifelse(above_horizon==FALSE,0,phase),
  fraction= ifelse(above_horizon==FALSE,0,fraction),
  l.illum= ifelse(above_horizon==FALSE,0,l.illum)
)


# merge -------------------------------------------------------------------



namestraits<-left_join(btrait, batnames, by= "sp")
colnames(namestraits)[2]<-"sp4"
traits<- namestraits %>% dplyr::select(c("sp4","Six.letter.species.code","Mass","Aspect","PeakFreq","Loading" ) )
colnames(traits)[2]<-"sp"
rm(namestraits) 

filtered_bm<- left_join(filtered_bm, traits, by="sp")

bm2<-left_join(filtered_bm, elev, by="site" )
bm2 <- left_join(bm2, ndvi, by = "site") %>%
  select(-c("X", "time", "buff_area.x", "buff_area.y"))

colnames(bm2)[7]<-"date"

bm2$date<-as_date(bm2$date)
bm2<- left_join(bm2, weather, by="date")


bm2<- left_join(bm2, moon.adj, by=c("date"))


# Identify rows with NA values
rows_with_na <- bm2[rowSums(is.na(bm2)) > 0, ] # some wind and temp have NAs because we are missing august 2021 weather data.

# correlation  ------------------------------------------------------------


# check for correlation 
numeric_cols<- sapply(bm2, is.numeric)
cor1<-bm2[,numeric_cols] #keeps just the numeric
cor1<-cor1 %>% select(-c("X_min","X_max","altitude","parallacticAngle","angle","lat", "lon" ))

c1<- cor(cor1,use="pairwise.complete.obs")
corrplot(c1, order= 'AOE', tl.)


# scale data 

# bm2$elevm_s<- scale(bm2$elev_mean, center = T, scale = T) Use the loop instead. 
# bm2$ndvim_s<- scale(bm2$ndvi_mean)
# bm2$percent_s<-scale(bm2$percent)
# bm2$jday_s<-scale(bm2$jday)
# bm2$PeakFreq_s<-scale(bm2$PeakFreq)
# bm2$l.illum_s<-scale(bm2$l.illum)
# bm2$wind_s<scale(bm2$avg_wind_speed)

# List of variable names to be scaled
variables_to_scale <- c("elev_mean", "ndvi_mean", "percent", "jday", 
                        "PeakFreq", "l.illum", "avg_wind_speed","avg_temperature","phase","fraction" )

# Loop over each variable, scale it, and assign it back to the data frame with a new name
for (var in variables_to_scale) {
  bm2[[paste0(var, "_s")]] <- scale(bm2[[var]], center = TRUE, scale = TRUE)
}

# explore data ------------------------------------------------------------

# Plot the distribution of the count data
ggplot(filtered_bm, aes(x = n)) +
  geom_histogram(binwidth = 40, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Bat Calls", x = "Number of Calls", y = "Frequency")

ggplot(filtered_bm, aes(x=jday, y=n, col=treatmt))+
  geom_point()+
  facet_wrap("site")


variables_to_plot <- c("elev_mean", "ndvi_mean", "percent", "jday", 
                        "PeakFreq", "l.illum", "avg_wind_speed","avg_temperature","phase","fraction" )
# Generate the plots
plots <- lapply(variables_to_plot, function(var) {
  if (is.numeric(bm2[[var]])) {
    ggplot(bm2, aes_string(x = var)) +
      geom_histogram(fill = "blue", color = "black") +
      theme_minimal() +
      ggtitle(paste("Histogram of", var))
  } else {
    message(paste("Skipping non-numeric variable:", var))
    NULL
  }
})
# Combine the plots using patchwork
plots <- Filter(Negate(is.null), plots)# removes all the null or NAs from the plot. We have somefor missing data. in avg temp, wind and peak freq
combined_plot <- wrap_plots(plots)
# Print the combined plot
print(combined_plot)

# models -------------------------------------------------------------------

# Fit a Poisson model intercept model. 
m1 <- glmer(n ~ (1|site),
                     data = bm2, 
                     family = poisson)
exp(m1@beta)# base line for bats at all sites. 

summary(m1)



m1.2<- glmer(n ~ trmt_bin * jday_scaled + percent_s +elev_mean_s+  (1|site),
                   data = bm2, 
                   family = poisson)

plot_model(m1.2)

exp(m1.2@beta)
summary(m1.2)

plot(fitted(m1.2), residuals(m1.2), main = "Residuals vs Fitted", 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(residuals(m1.2), main = "Normal Q-Q Plot")
qqline(residuals(m1.2), col = "red")

# Histogram of residuals
hist(residuals(m1.2), breaks = 30, main = "Histogram of Residuals", 
     xlab = "Residuals")

# Scale-Location Plot
plot(fitted(m1.2), sqrt(abs(residuals(m1.2))), main = "Scale-Location Plot",
     xlab = "Fitted values", ylab = "Square Root of |Residuals|")
abline(h = 0, col = "red")

# Calculate VIF
vif(m1.2)



# quasipoisson -------------------------------------------------------

#trying to make a quasi poisson.

m1.2 <- glmer(
  n ~ trmt_bin * jday_s + ndvi_mean_s + percent_s + PeakFreq_s + l.illum_s +
    avg_wind_speed_s + avg_temperature_s + (1+sp|site),
  data = bm2,
  family = 'poisson',
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

summary(m1.2)

exp(coef(m1.2))

plot_model(m1.2,)

model_convergence <- m1.2@optinfo$conv$opt
print(model_convergence)

# Check if the Hessian matrix is positive definite
is_positive_definite <- all(eigen(m1.2@optinfo$derivs$Hessian)$values > 0)
print(is_positive_definite)
