# ---------------------------
##
## Script name:  guano_v1.R
##
## Purpose of script: Access GUANO metadata in audio files. 
##
## Author: Carlos Linares
##
## Date Created: 05/21/2024
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes: Trying to developing versions for to guano metadata processing. 
##   
##
## ---------------------------
## # A bat file with data that has been run trhough sonobat for buzz. detection. 
##


# outputs ----------------------

# 

# this should be a database ready to analyze with the glmm_v1 script. 


#set wd
setwd("Z:/PioneerLights_2021/")


# libraries  --------------------------------------------------------------

library(httr)
set_config(config(ssl_verifypeer=0L))

install.packages("devtools")
devtools::install_github("riggsd/guano-r", subdir="guano")

library("guano")
dirname<-"data_for_analysis/IRON01__1__20230706_203741.wav"

a<-read.guano(dirname)

# the library above gave me trouble installing and was not able to produce the output I wanted.

# install.packages("devtools")
# devtools::install_github("vulpes-vulpes/batr")

library(batr)

# path to the z drive new suggestion of code
pp<-"data_for_analysis/IRON01__1__20230706_203741.wav"
pp2<-"E:/pioneer2023/robomoth_2023_all"

GUANO_reader(pp)
batr::read.guano('data_for_analysis/sample_wavs/IRON01_20210729_232359-Mylu.wav') # this function did not work...

GUANO_reader(pp2, "robomoth_2023") # function to be retired soon. # this worked
import_GUANO("New",pp2,site_col = "pioneer",site_col="user","E:/pioneer2023/robomoth_2023_all/sonobat_output")
