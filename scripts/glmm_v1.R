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
library(effects)
library(reshape2)


#load environment 
#last worked 08/16/2024
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
   n ~ trmt_bin + jday_s + I(jday_s ^ 2) + percent_s  + l.illum_s +  # less complex model 
     avg_wind_speed_s + avg_temperature_s + yr_s + elev_mean_s +
     (1 |site) + (1 | sp),
   data = bm2,
   nbinom2(link = "log"))

 
m1.5nb <- glmmTMB(
  n ~ trmt_bin + jday_s + I(jday_s ^ 2) + percent_s  + l.illum_s +
    avg_wind_speed_s + avg_temperature_s + yr_s + elev_mean_s +
    (1 |site) + (1 + trmt_bin + jday_s + I(jday_s ^ 2) | sp),
  data = bm2,   nbinom2(link = "log")
)

m1.6nb <- glmmTMB(
  n ~ trmt_bin + jday_s + I(jday_s ^ 2) + percent_s  + l.illum_s + # this one didn't even converge. 
    avg_wind_speed_s + avg_temperature_s + yr_s + elev_mean_s +
    (1 |site) + (1 + trmt_bin + jday_s + I(jday_s ^ 2) | sp),
  data = bm2,   nbinom1(link = "log")
)
save(m1.5nb, 'models/')
AIC( m1.5nb,m1.5anb,m1.6nb) 
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
load('models/m1.5nb')

tab_model(m1.5nb)

emmeans::emmeans(m1.5nb, "percent_s",type="response")

plot(ggeffects::ggpredict(m1.5nb, "percent_s [all]"))
plot(ggeffects::ggpredict(m1.5nb, "l.illum_s [all]"))
plot(ggeffects::ggpredict(m1.5nb, "avg_wind_speed_s [all]"))
plot(ggeffects::ggpredict(m1.5nb, "avg_temperature_s [all]"))




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
  ylab("bat calls") +
  xlab("jday") +
  geom_line(size = 1.5) 


# Create the melted data for plotting# Create the melted data for plotting# Create the melted data for plotting all columns
abunddf_melted <- melt(abunddf, id.vars = "ord.day",variable.name = "sp", value.name ="predicted calls" , measure.vars = 1:14)
unique(abunddf_melted$sp)

abunddf_long<- pivot_longer(abunddf, cols= 1:14, names_to = "sp", values_to = "predicted calls")
unique(abunddf_long$sp)
head(abunddf_long)

# Create the ggplot object
p6<-ggplot(abunddf_melted, aes(x = ord.day, y = `predicted calls`, color = sp)) +
  scale_color_viridis_d()+
  # Add geom_point to plot points for each variable
  scale_shape_manual(values = 1:6)+ 
  geom_point(size = 2) +
    # Use different shapes
  # Set labels and title
  labs(title = "Abundance by Day", x = "Day", y = "bat calls") +
  # Set theme for better visuals (optional)
  theme_classic()
p6
ggplot(abunddf_long, aes(x = ord.day, y = `predicted calls`, color = sp)) +
  # Add geom_point to plot points for each variable
  geom_point(size = 1) +
  facet_wrap( ~ sp, scales = "free_y") +
  # Set labels and title
  labs(title = "Bat calls by julian day", x = "julian day", y = "bat calls") +
  # Set theme for better visuals (optional)
  theme_classic()
rm(abunddf)



p6
ggsave(filename = "jday_raneff_v1.png",plot = p6,device = "png", path = 'figures/glmm_v1/' )


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
  scale_color_viridis_d(direction = -1, option = "A")+
  scale_shape_manual(values = 1:14) +  # Customize shapes for each 'sp', adjust as needed
  ylab( "bat calls" ) +
  xlab( "year" ) +
  geom_point(size=3)+
  geom_line()+
  # labels
  ylab("Bat calls") +
  xlab("") 
p7


ggsave(filename = "trmt_raneff_v1.png",plot = p7,device = "png", path = 'figures/glmm_v1/' )



# save models -------------------------------------------------------------

#save image 
save.image(file = "working_env/glmm_v1.RData")

save(m1.5nb, file = "models/my_models.RData")
load("models/my_models.RData")




# trash -------------------------------------------------------------------

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

