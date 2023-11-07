

##     created by Carlos Linares            ##
##     
##     Purpose: Replicate analysis from the 
##     https://github.com/quantitativeconservationlab/AppPopnEco/blob/master/CountBayesAnalysis.R
##     
##     and https://github.com/quantitativeconservationlab/AppPopnEco/blob/master/CountBayesPlotting.R
##     

#libraries 
library(tidyverse)
library(jagsUI)

# load data
# this is the data from bat_pop_analysis script and was cleaned in the cleanup_script

js21<-read.csv('data_analysis/bat_pop_analysis/bat_js.csv')

#sites
I<-length(unique(js21$site))

#visits
j<-10 # thast the number of weeks 

# what varibles we scale?




