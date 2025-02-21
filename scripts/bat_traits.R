#----------------- Bat specie trait prep
# ------- purpose : get the traits for the specie of bats present in the sites. 
# 
# author: Carlos Linares
# date: Enero 2024


# Notes: there are some species missing so we will have to add them. 
# 
# inputs: graph created by Dr.Corsini with the forearm/ear ratio.
# outputs: data_for_analysis/bat_traits/bat_eararm_v2.csv has the complete data species in the pionner project.

#------------  load libraries

pacman::p_load(
  tidyverse, # for data manipulation
  magrittr, # for the pipe operator %>%
  waldo, # compare objects in R
  data.table, # for fread and fwrite
  ggplot2 # for plotting
)

# load bats from pioneers.

bm <- read_csv('data_for_analysis/prep_for_glmm/bm.csv')

# get species pioneer project

bats.pioneer<- bm %>%
  select(AUTO.ID.) %>%
  distinct() %>%
  arrange(AUTO.ID.) %>% 
  filter(!AUTO.ID.==c("Noise", "NoID"))


# load Species bat data

spcodes<-fread("data_for_analysis/Species_bats.csv")

# merge spcodes and pioneer bats

bats.pioneer<-merge(spcodes, bats.pioneer, by.x = "Six-letter species code", by.y = "AUTO.ID.")
# change col names 
colnames(bats.pioneer)[1]<-"six_sp"
colnames(bats.pioneer)[2]<-"cname"
colnames(bats.pioneer)[3]<-"sname"
colnames(bats.pioneer)[4]<-"four_sp"



# now we can  calculate the ratio ear-length to forearm length provided by Dr.Corsini. 

# Create data frame
bat_trait <- data.frame(
  Species = c("Silver-haired bat", "Hoary bat", "Little brown bat", "Big brown bat", 
              "Long-legged Myotis", "Yuma myotis", "California myotis", 
              "Long-eared Myotis", "Pallid bat","Western small-footed myotis","Long-legged myotis"),
  Ear_Forearm_Ratio = c(
    0.22, 0.25, 0.28, 0.30, 0.34, 0.36, 0.38, 0.42, 0.50, NA, NA)
)



# Merge the bat trait and bat.pioneer

bat_trait<- left_join(bats.pioneer, bat_trait, by=c("cname" = "Species"))

colnames(bat_trait)[5]<-"ear.arm"

# Recreate Dr.Corsini plot

p1<-ggplot(bat_trait, aes(x = reorder(cname, ear.arm), y = ear.arm)) +
  geom_point(size = 4, fill = "black", color = "black") +  # Shapes with fill
  theme_minimal() +
  labs(x = "", y = "Ear length / Forearm length", title = "bat trait ear/arm ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave("figures/bat_traits/ear_arm_ratio.png", p1, width = 10, height = 5, units = "in", dpi = 300)


# the table below was manually modified and completed. bat_eararm_v2.csv
# fwrite(bat_trait, 'data_for_analysis/bat_traits/bat_eararm.csv') 




# complete the missing species with the mean of the ear/arm ratio.

# Create the dataset
bat_data <- data.frame(
  Species = c("Myotis evotis", "Myotis thysanodes", "Myotis volans", "Parastrellus hesperus"),
  Ear_Length_mm = c("18–25", "16–21", "11.8", "11–15"),
  Forearm_Length_mm = c("36–41", "39–46", "35.2", "28–34"),
  Ear_to_Forearm_Ratio = c("0.50–0.61", "0.35–0.54", "0.34", "0.32–0.54"),
  Reference = c(
    "https://centralcoastbatsurvey.org/our-central-central/",
    "https://cnhp.colostate.edu/download/documents/cbwg-pdfs/ConPlanRevisionFinal/Chapter%2013%20CBWG_BatConservationPlan_2ndEdition_2018_SpeciesAccounts.pdf",
    "https://animaldiversity.org/accounts/Myotis_volans/",
    "https://www.batsnorthwest.org/meet_washingtons_bats.html"
  )
)



# update Dr. Corsini plot 

fread('data_for_analysis/bat_traits/bat_eararm_v2.csv') -> bat_trait

p2<-ggplot(bat_trait, aes(x = reorder(cname, ear.arm), y = ear.arm)) +
  geom_point(size = 4, fill = "black", color = "black") +  # Shapes with fill
  theme_minimal() +
  labs(x = "", y = "Ear length / Forearm length", title = "bat trait ear/arm ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2

p3<-p2+theme_sjplot(base_size = 20, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "figures/bat_traits/ear_arm_ratio_v3.png",
  p3,
  device = "tiff",
  width = 10,
  height = 5,
  units = "in",
  dpi = 300
)


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
