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



# sm3 buzz data -----------------------------------------------------------


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


# remove year buzz files to free memory
rm(sm3_buzz_2021, sm3_buzz_2022, sm3_buzz_2023)
gc()

# acoustic out 
# here we remove the acoustic files marked with a _0_ in the file name 

sm3_buzz_all<- sm3_buzz_all %>%
  filter(!str_detect(filename, "__0__"))

# columns to keep

keep <- c(
  "filename",
  "source_file",
  "spp",
  "spp_accp",
  "species_manual_id",
  "auto_buzz_count",
  "manual_buzz_count",
  "hi_f",
  "lo_f",
  "year",
  "monitoring_night"
)

sm3_buzz_all<- sm3_buzz_all %>%
  select(all_of(keep))

glimpse(sm3_buzz_all)



# check for duplicate filename rows 
# it seems theres no duplicate rows. 

buzz_dup_summary <- sm3_buzz_all %>%
  group_by(filename) %>%
  summarise(
    n_rows = n(),
    n_source_files = n_distinct(source_file),
    source_files = paste(sort(unique(source_file)), collapse = "; "),
    years = paste(sort(unique(year)), collapse = "; "),
    monitoring_nights = paste(sort(unique(monitoring_night)), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(n_rows > 1) %>%
  arrange(desc(n_rows))




# site  -------------------------------------------------------------------


# we create a site column based on the file name. 

sm3_buzz_all <- sm3_buzz_all %>%
  mutate(
    site = str_to_lower(str_extract(filename, "^[A-Za-z]{3,4}[0-9]{2}"))
  )
head(sm3_buzz_all$site)


# check sites 
unique(sm3_buzz_all$site)
sum(is.na(sm3_buzz_all$site)) 

sm3_buzz_all %>%
  filter(is.na(site)) %>%
  count(filename, source_file, sort = TRUE)


# check if we have any NAs in the site column. We have 11 why? these files are files that show some problem with the file name tag being improperly written. 
# We will filter them out as these are files that we didn't use in the analysis.

# remove NA for site column 

sm3_buzz_all <- sm3_buzz_all %>%
  filter(!is.na(site))


# standardize site names

site_key <- c(
  "lon01" = "long01",
  "lon02" = "long02",
  "lon04" = "long04",
  "lon05" = "long05",
  "vis01" = "vizc01",
  "vis02" = "vizc02",
  "vis03" = "vizc03",
  "vis04" = "vizc04",
  "viz01" = "vizc01",
  "viz03" = "vizc03",
  "viz04" = "vizc04"
)

sm3_buzz_all <- sm3_buzz_all %>%
  mutate(
    site = recode(site, !!!site_key)
  )

# check sites again

unique(sm3_buzz_all$site)

sum(is.na(sm3_buzz_all$site)) # check if we have any NAs in the site column. We have 0 so we are good.

# now we correct the buzz column and create a single column. We give priority to the manual buzz count if it exists, otherwise we use the auto buzz count.
# will also use the manual buzz count to assess how many files we have vetted manually. 

sm3_buzz_all <- sm3_buzz_all %>%
  mutate(
    auto_buzz_count = parse_number(auto_buzz_count),
    manual_buzz_count = parse_number(manual_buzz_count),
    c_buzz = coalesce(manual_buzz_count, auto_buzz_count)
  )

summary(sm3_buzz_all)

# we have some NAs in the c_buzz and the auto_buzz_count columns. We will check where are these coming from. It could be possible we haven't ran all the files through the buzz tool in sonobat.


t<- sm3_buzz_all %>%
  filter(is.na(c_buzz)) %>%
  select(
    filename,
    source_file,
    auto_buzz_count,
    manual_buzz_count,
    c_buzz,
    spp,
    spp_accp,
    species_manual_id,
    year,
    monitoring_night,
    site
  )
table(t$source_file)

# iron01 has 11342; files that can be removed becuase these are 0kb file so no information on them.
# iron0_2022 has one file with 0kb
# iron03_2021 has 2038 files with 0kb to remove. 
# long01_2023 has 116 files with 0kb to remove. 
# vizc01_2022 has 71 fiels that can be removed
# vizc03_2021 has 61 files that have data but can be removed 
# I took the decision to remove the ones that have one file or three. 

sm3_buzz_all <- sm3_buzz_all %>%
  filter(!is.na(c_buzz))

# How many rows have manual buzz count and manual species manual ID. 2599 taht is very little. 

count_manual <- sm3_buzz_all %>%
  filter(!is.na(manual_buzz_count) & !is.na(species_manual_id)) %>%
  nrow()
count_manual

# date time

sm3_buzz_all <- sm3_buzz_all %>%
  mutate(
    date_time = ymd_hms(str_extract(filename, "\\d{8}_\\d{6}"), tz = "America/Denver")
  )

sum(is.na(sm3_buzz_all$date_time))  # there is two NAs but they come from an error in the file name so these can be removed.

t<- sm3_buzz_all %>%
  filter(is.na(date_time)) %>%
  select(
    filename,
    source_file,
    auto_buzz_count,
    manual_buzz_count,
    c_buzz,
    spp,
    spp_accp,
    species_manual_id,
    year,
    monitoring_night,
    site
  )

sm3_buzz_all <- sm3_buzz_all %>%
  filter(!is.na(date_time))

names(sm3_buzz_all)


# tere are no missing buzz auto count files. 
missing_auto <- sm3_buzz_all %>%
  filter(
    is.na(auto_buzz_count) |
      auto_buzz_count == "" |
      auto_buzz_count == "NA"
  ) %>%
  select(filename, auto_buzz_count, source_file )

table(missing_auto$source_file)

# I believe this is ready to be saved as a csv file. Then we will bind the cbuzz to the other data used for the glmms.
# 

write_csv(sm3_buzz_all, file = 'data_for_analysis/sm3_buzz_build/sm3_buzz_all.csv') # raw combine data 

file.exists("data_for_analysis/sm3_buzz_build/sm3_buzz_all.csv")

# Create a README file with information about the script
readme_content <- "Carlos Linares 7/07/2026 

this folder contains the sm3 data outputs for buzz counts. All files have been analyzed using the buzz tool in the sonobat software. The data has been cleaned and combined into a single csv file called sm3_buzz_all.csv.

glimpse(sm3_buzz_all)
Rows: 1,094,812
Columns: 14
$ filename          the file name for the recording 
$ source_file       the output file from the vetting process. 
$ spp               the possible species ID from the sonobat software.
$ spp_accp          the species accepted by sonobat algorithm
$ species_manual_id the species assigned after manual vetting of the files.
$ auto_buzz_count   the number of buzzes detected by the sonobat software.
$ manual_buzz_count <dbl> corrected number of buzzes manually 
$ hi_f              <chr> the numbe of high frequency buzzes detected by the sonobat software.
$ lo_f              <chr> the number of low frequency buzzes detected by the sonobat software.
$ year              <int> year
$ monitoring_night  <chr> the monitoring night when a recording has times past midnight assigns the date of the previous night. 
$ site              <chr> the site where the recording was made.
$ c_buzz            <dbl> corrected number of buzzes using the manual and auto columns. 
$ date_time         <dttm> 2021-08-07 22:50:52, 2021-08-07 23:01:48, 2021-08-11 23:00:27, 2021-08-07 23:27:13, 2021-08-11 22:33:47, 2021-08-07 23:27:29, 2021-08-07 22:45:32, 2021-08-07 22:57:00, 202…


"

# Write the README content to a file
writeLines(readme_content, "data_for_analysis/sm3_buzz_build/README.txt")




