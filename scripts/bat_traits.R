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
library(data.table)



# now we cant to caculate the ratio earlenght forarm length provided by Dr.Corsini. 

# Load necessary libraries
library(ggplot2)

# Create data frame
bat_trait <- data.frame(
  Species = c("Silver-haired", "Hoary bat", "Little brown bat", "Big brown bat", 
              "Long-legged Myotis", "Yuma myotis", "Myotis californicus", 
              "Long-eared Myotis", "Pallid bat"),
  Ear_Forearm_Ratio = c(0.22, 0.25, 0.28, 0.30, 0.34, 0.36, 0.38, 0.42, 0.50)
)

# Create the plot
ggplot(bat_trait, aes(x = reorder(Species, Ear_Forearm_Ratio), y = Ear_Forearm_Ratio)) +
  geom_point(size = 5, fill = "black", color = "black") +  # Shapes with fill
  theme_minimal() +
  labs(x = "Species", y = "Ear length / Forearm length", title = "Ear Length to Forearm Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fwrite(bat_trait, 'data_for_analysis/bat_eararm.csv')





# trash

# 
# #---- load data ---
# # here we load the 2021 bat data to get the species names.
# 
# bat2021<-read.csv('data_for_analysis/bat2021_v3.csv',check.names = T)
# 
# species<-unique(bat2021$SppAccp)
# 
# # now we load bat species traits data bases Dylans and Michela (if available)
# 
# dylan<-read.csv('data_for_analysis/Bat_trait.csv') # this trait list has data from Zou et al. 2022 (doi: 10.3389/fevo.2022.1031548)
# 
# identical(species, dylan$Species)
# all(sort(species) == sort(dylan$Species)) # species in dylans paper and mine are the same.
# 
