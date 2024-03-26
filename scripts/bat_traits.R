#----------------- Bat specie trait prep
# ------- purpose : get the traits for the specie of bats present in the sites. 
# 
# author: Carlos Linares
# date: Enero 2024
# 

#------------  load libraries
library(tidyverse)
library(magrittr)
library(waldo)#compare objects in R

#---- load data ---
# here we load the 2021 bat data to get the species names.

bat2021<-read.csv('data_for_analysis/bat2021_v3.csv',check.names = T)

species<-unique(bat2021$SppAccp)

# now we load bat species traits data bases Dylans and Michela (if available)

dylan<-read.csv('data_for_analysis/Bat_trait.csv') # this trait list has data from Zou et al. 2022 (doi: 10.3389/fevo.2022.1031548)

identical(species, dylan$Species)
all(sort(species) == sort(dylan$Species)) # species in dylans paper and mine are the same.

