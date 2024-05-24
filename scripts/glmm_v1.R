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
## Notes: 
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


# data --------------------------------------------------------------------


bm<-read.csv('data_for_analysis/data_glmm/bat_counts.csv', header = T) # bat counts by jday
bm$jday_scaled<-scale(bm$jday)
bm$dateytyt<-as_date(bm$jday, origin= "2021-01-01")

btrait<-read.csv('data_for_analysis/Bat_trait.csv', header = T)
btrait$Species<-toupper(btrait$Species)
batnames<-read.csv('data_for_analysis/Species_bats.csv')
colnames(batnames)[3]<-"Species"

elev<-read.csv('data_for_analysis/elev/elevation.csv', header = T)
elev<-elev %>% rename("site" = "name")

ndvi<-read.csv('data_for_analysis/NDVI/NDVI_of_rip2021.csv', header = T)
weather<-read.csv('data_for_analysis/weather/dailyavg.csv', header = T)
moon<-read.csv('data_for_analysis/moon_pred.csv')

  

#merge 
namestraits<-left_join(btrait, batnames, by= "Species")
traits<- namestraits %>% dplyr::select(c("Species","Six.letter.species.code","Mass","Aspect","PeakFreq","Loading" ) )

t<-left_join(bm, elev, by="site", )

# explore data ------------------------------------------------------------

# Plot the distribution of the count data
ggplot(bm, aes(x = n)) +
  geom_histogram(binwidth = 40, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Bat Calls", x = "Number of Calls", y = "Frequency")

ggplot(bm, aes(x=jday, y=n, col=treatmt))+
  geom_point()+
  facet_wrap("site")


# model -------------------------------------------------------------------

# Fit a Poisson model without random effects to check goodness-of-fit
m1 <- glmer(n ~ (1|site),
                     data = bm, 
                     family = poisson)
exp(m1@beta)# base line for bats at all sites. 

summary(poisson_model)

m1.1<- glmer(n ~ trmt_bin +(1|site),
           data = bm, 
           family = poisson)
exp(m1.1@beta)
summary(m1.1)


m1.2<- glmer(n ~ trmt_bin + jday_scaled+ (1|site),
                   data = bm, 
                   family = poisson)


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



# negative binomial -------------------------------------------------------


library(MASS)
m1.2_nb <- glmer.nb(n ~ trmt_bin + jday_scaled + (1|site), data = bm)
summary(m1.2_nb)

plot(fitted(m1.2_nb), residuals(m1.2_nb), main = "Residuals vs Fitted", 
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
