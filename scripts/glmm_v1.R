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


# data --------------------------------------------------------------------


bm<-read.csv('data_for_analysis/data_glmm/bat_counts.csv', header = T)


# model -------------------------------------------------------------------

null <- glmer(n ~ (1 | site), data = bm,
                      family = poisson(link = "log"))
