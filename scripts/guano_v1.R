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




# another try at guano function to work. This one works. It is a little tricky to install guano. 
# 
library(guano)
dt1 <- read.guano.dir(dirname = 'data_for_analysis/sample_wavs/' ,
                      pattern = "*.wav" ,
                      recursive = T) # this code from guano in R worked.

# now I want to see if the files that have been ID with kp and sono bat show or if running Kp erase sonobat. 

# I ran the same files but with kaleidoscope to see if the process changed the metadata. 

execution_time <- system.time({  #processing 882 files took 83 secs. Maybe I can do all of them...
  dt2 <- read.guano.dir(dirname = 'Z:/PioneerLights_2021/', recursive = FALSE)
})

# Print the execution time
print(execution_time)

# now let's try with the bioacustics package
# the following code work but does not provide an easy way to save the metadata into R. 
# I also learnt that Kaleidoscope metadata is not recovered this way unless it writes audio files. 



install.packages("bioacoustics")

library(bioacoustics)
t<-read_audio('data_for_analysis/sample_wavs/IRON01_20210729_232359-Mylu.wav')
class(t)
a<-metadata("data_for_analysis/sample_wavs/IRON01_20210729_232359-Mylu.wav")
a<-metadata(t)


# conclusion: seems like a working function is the read.guano.dir() rom the guano package. 





# trash -------------------------------------------------------------------

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



# Another package that did not work. 
library(batr)

# path to the z drive new suggestion of code
pp<-"data_for_analysis/IRON01__1__20230706_203741.wav"
pp2<-"E:/pioneer2023/robomoth_2023_all"

GUANO_reader(pp)
batr::read.guano('data_for_analysis/sample_wavs/IRON01_20210729_232359-Mylu.wav') # this function did not work...

<<<<<<< HEAD
GUANO_reader(pp, "pioneer") # function to be retired soon. # this worked
import_GUANO("New","data_for_analysis/sample_wavs/", 
             "C:/Users/Carlos/Documents/R/bats_pioneer/data_for_analysis/sample_wavs/batr_out")


# guano from python


library(reticulate)
use_python("C:/Users/Carlos/AppData/Local/Programs/Python/Python38/python.exe")  # Specify your Python path if needed
guano <- import("guano")

# Load a .WAV file with GUANO metadata
g <- guano$GuanoFile('C:/Users/Carlos/Documents/R/bats_pioneer/data_for_analysis/sample_wavs/IRON01_20210729_232359-Mylu.wav')

# Access metadata
version <- g['GUANO|Version']
make <- g['Make']
model <- g['Model']

print(version)
print(make)
print(model)
=======
  GUANO_reader(pp2, "robomoth_2023") # function to be retired soon. # this worked
import_GUANO("New",pp2,site_col = "pioneer",site_col="user","E:/pioneer2023/robomoth_2023_all/sonobat_output")
>>>>>>> 3e24bdc4bc8d35895cb77181dd2b80757808ad45
