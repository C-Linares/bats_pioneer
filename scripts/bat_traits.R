#----------------- Bat specie trait prep
# ------- purpose : get the traits for the specie of bats present in the sites. 
# 
# author: Carlos Linares
# date: Enero 2024


# Notes: there are some species missing so we will have to add them. 

#------------  load libraries

pacman::p_load(
  tidyverse, # for data manipulation
  magrittr, # for the pipe operator %>%
  waldo, # compare objects in R
  data.table, # for fread and fwrite
  ggplot2 # for plotting
)


# now we cant to calculate the ratio ear-lenght to forearm length provided by Dr.Corsini. 

# Load necessary libraries

# Create data frame
bat_trait <- data.frame(
  Species = c("Silver-haired bat", "Hoary bat", "Little brown bat", "Big brown bat", 
              "Long-legged Myotis", "Yuma myotis", "Myotis californicus", 
              "Long-eared Myotis", "Pallid bat", "Silver-haired bat","Western small-footed myotis","Long-legged myotis",""),
  Ear_Forearm_Ratio = c(0.22, 0.25, 0.28, 0.30, 0.34, 0.36, 0.38, 0.42, 0.50)
)

# load Species bat data
# this file contains the four letter species codes.

spcodes<-fread("data_for_analysis/Species_bats.csv")


# Merge the bat trait and spcodes by species

bat_trait <- merge(bat_trait, spcodes, by.x = "Species", by.y = "Common name")

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
