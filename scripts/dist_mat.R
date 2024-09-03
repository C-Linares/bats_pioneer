
# Script name: Prep_for_glmm
# 
# Purpose of script: Combine the 2021, 2022, 2023, data sets and get them ready to go throuhg glmm_v1
# 
# Author: Carlos Linares
# 
# Date Created: 09/03/2024
# 
# Email: carlosgarcialina@u.boisestate.edu
# 
# ---------------------------
#   
#   Notes: short script to add distance to lights as a variable instead of dark or light. the output should be added in the glmm_v1 as a coovariat/ predictor. 


# inputs ------------------------------------------------------------------
#  Home>R>bats_pioneer>data_for_analysis>site_dist_matrix.csv
# 
# outputs ----------------------

# tbd 


# libraries

library(tidyverse)
library(beepr)



# data --------------------------------------------------------------------

dst.mat<- read.csv('data_for_analysis/site_dist_matrix.csv')

# make all site labels small caps

dst.mat$InputID<-tolower(dst.mat$InputID)
dst.mat$TargetID<-tolower(dst.mat$TargetID)
head(dst.mat)

# add the two distances. 

sum.dst <- dst.mat %>%
  group_by(InputID) %>%
  summarize(
    TargetIDs = str_c(TargetID, collapse = ","),
    Sum_Distance = sum(Distance),
    .groups = 'drop'
  )

colnames(sum.dst)[1]<-"site"

sum.dst <- sum.dst %>%
  mutate(site = if_else(str_starts(site, "viz"), 
                        str_replace(site, "^viz", "vizc"), 
                        site))

sum.dst

write.csv(sum.dst, file = 'data_for_analysis/sum.dst.csv', row.names = F)

