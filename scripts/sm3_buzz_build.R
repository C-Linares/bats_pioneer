#=============================================================================
# Script Name:    sm3_buzz_build.R
# Purpose:        build buzz files for all the sm3 data across all sites. 
#
# Version:        1
#
# Author:         Carlos Linares
# Collaborator:   
# Created:        2026-06-30
# Contact:        carlosgarcialina@u.boisestate.edu
#
# Notes:
#   - original data stored in the bk_linares and PhD_bk hard drives
#   - we ran all the files through the sonobat detector but the IDs are missing in several files. we will match the ids from past databases. Possibly using the files from the prep_for_glmm_v2.R
#
#
# Inputs:
# - 2021 
# 
#
# Outputs:
#   - data bases cleaned for analysis
# =============================================================================
# 
# What I want to do with the data for the buzz is to merge it all into a database of buzz only and then pair ir with the ID from Kpro taht we have been using for the full analysis. 

# load packages

if(!require("pacman")) install.packages("pacman")
pacman::p_load(
  "tidyverse",
  "lubridate",
  "stringr",
  "ggplot2",
  "corrplot",
  "glmmTMB",
  "performance",
  "marginaleffects",
  "DHARMa",
  "sjPlot",
  "patchwork",
  "janitor"
)

# data --------------------------------------------------------------------

# first robomoth data paths

path_2021 <- "G:/PioneerLights_2021/buzz_out/"   # single file
path_2022 <- "G:/PioneerLights_2022/sm3/buzz_out/" # folder
path_2023 <- "G:/PioneerLights_2023/sm3/sm3/buzz/"             # folder

# lets list the files inside the paths. 

# function -----------------------------------------------------------------


#function to read files 
# when mergin files the parent directory col was being read as character or number depending on the file.Thus we modified the function to always read that column as character 

read_buzz_file <- function(file, year) {
  read_tsv(
    file,
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  ) %>%
    clean_names() %>%
    mutate(
      year = as.integer(year),
      source_file = basename(file)
    )
}



# robomoth data -----------------------------------------------------------


# read 2021 

sm3_buzz_2021 <- read_buzz_file(path_2021, 2021)



# read 2022 and 2023

files_2021 <- list.files(path_2021, pattern = "\\.txt$", full.names = TRUE)
files_2022 <- list.files(path_2022, pattern = "\\.txt$", full.names = TRUE)
files_2023 <- list.files(path_2023, pattern = "\\.txt$", full.names = TRUE)

# read and combine files

sm3_buzz_2021 <- map_dfr(files_2021, read_buzz_file, year = 2021)
sm3_buzz_2022 <- map_dfr(files_2022, read_buzz_file, year = 2022)
sm3_buzz_2023 <- map_dfr(files_2023, read_buzz_file, year = 2023)


#combine robomoth --------------------------------------------------------

# combine all years into one data frame
sm3_buzz_all <- bind_rows(sm3_buzz_2021, sm3_buzz_2022, sm3_buzz_2023)

str(sm3_buzz_all)


