
#=============================================================================
# Script Name:    robomoth_build.R
# Purpose:        build rbomoth and speaker database
#
# Version:        1
#
# Author:         Carlos Linares
# Collaborator:   
# Created:        2026-03-17
# Contact:        carlosgarcialina@u.boisestate.edu
#
# Notes:
#   - original data stored in the bk_linares and PhD_bk hard drives
#
# Inputs:
# - paths to data
#
# Outputs:
#   - data bases cleaned for analysis
# =============================================================================
# 

# libraries ---------------------------------------------------------------

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

path_2021 <- "F:/pioneer_2021/sm4/sb_buzz_2021_vetted.txt"   # single file
path_2022 <- "G:/PioneerLights_2022/part2/robomoth/buzz_out" # folder
path_2023 <- "F:/pioneer_2023/robomoth_2023_all"             # folder

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

# read 2021 

robomoth_2021 <- read_buzz_file(path_2021, 2021)

# needs some fixes first create a new column for sp that inputs the spp_accp but if missing puts the species_manual_id or the spp column 
# we will do this at the end when we have all the files together and I am cleaning data. 



# read 2022 and 2023


files_2022 <- list.files(path_2022, pattern = "\\.txt$", full.names = TRUE)
files_2023 <- list.files(path_2023, pattern = "\\.txt$", full.names = TRUE)

# read and combine files

robomoth_2022 <- map_dfr(files_2022, read_buzz_file, year = 2022)
robomoth_2023 <- map_dfr(files_2023, read_buzz_file, year = 2023)



# combine robomoth --------------------------------------------------------

# combine all years into one data frame
buzz_all <- bind_rows(robomoth_2021, robomoth_2022, robomoth_2023)



# sp column ---------------------------------------------------------------


# create sp column
# we create this column because some observations have buzz count but not manual id we want to recover those IDs to not losse data. We build a new one combining the spp_accp, the manual id and the spp column and the 1st, 2nd, 3nd columns 



buzz_all.1 <- buzz_all %>%
  mutate(
    sp = case_when(
      !is.na(spp_accp) ~ spp_accp,
      is.na(spp_accp) & !is.na(species_manual_id) ~ species_manual_id,
      is.na(spp_accp) & is.na(species_manual_id) & !is.na(spp) ~ spp,
      TRUE ~ NA_character_
    )
  )

# site column 

# we create a site column based on the file name. 

buzz_all.1 <- buzz_all.1 %>%
  mutate(
    site = str_to_lower(str_extract(filename, "^[A-Za-z]{3,4}[0-9]{2}"))
  )

# check sites 
unique(buzz_all.1$site)
sum(is.na(buzz_all.1$site)) # check if we have any NAs in the site column. We have 0 so we are good.


# standardize site names

site_key <- c(
  "viz01" = "vizc01",
  "viz02" = "vizc02",
  "viz03" = "vizc03",
  "viz04" = "vizc04"
)

buzz_all.1 <- buzz_all.1 %>%
  mutate(
    site = recode(site, !!!site_key)
  )

# check sites again

unique(buzz_all.1$site)

sum(is.na(buzz_all.1$site)) # check if we have any NAs in the site column. We have 0 so we are good.


# create the date column
buzz_all.1 <- buzz_all.1 %>%
  mutate(
    date_time = ymd_hms(str_extract(filename, "\\d{8}_\\d{6}"), tz = "America/Denver")
  )

# we use the monitoring night for the date_s (standardize date column)

# after checking the difference I realize that my method creates some observations where the date is not corrected properly specially if close or around 12 am. So I will just use the monitoring night column as the date column for the analysis

buzz_all.1 <- buzz_all.1 %>%
  rename(noche = monitoring_night)

# time column 

buzz_all.1 <- buzz_all.1 %>%
  mutate(
    time = format(date_time, "%H:%M:%S")
  )

# checking for NAs in the auto buzz. this would indicate that that file has not been analyzed with the buzz tools from Sonobat. 

missing_auto <- buzz_all.1 %>%
  filter(
    is.na(auto_buzz_count) |
      auto_buzz_count == "" |
      auto_buzz_count == "NA"
  ) %>%
  select(filename, auto_buzz_count)


missing_files <- missing_auto %>%
  mutate(
    site = str_to_lower(str_extract(filename, "^[A-Za-z]{3,4}[0-9]{2}")),
    date = ymd(str_extract(filename, "\\d{8}"))
  )

files_to_check <- missing_files %>% # there are some files that need to be analyzed with Sonobat yet
  distinct(site, date)

# I want to see what columns don't have a species ID so I go back and check those. 

# check for missing species ID

missing_sp <- buzz_all.1 %>%
  filter(is.na(sp))

# the majority are zeros for buzz counts meaning there was nothing to ID. We might eliminate these

# check species column. 

unique(buzz_all.1$sp)

# first we make all small caps

buzz_all.1 <- buzz_all.1 %>%
  mutate(
    sp = str_to_lower(sp)
  )


# clean out ---------------------------------------------------------------

# make buzz counts numeric

buzz_all.1 <- buzz_all.1 %>%
  mutate(
    auto_buzz_count = as.numeric(auto_buzz_count),
    manual_buzz_count = as.numeric(manual_buzz_count),
  ) 

str(buzz_all.1)


# I realized I should not remove rows with NAs in the autobuzz counts. these are files that need to be analyzed with sonobat. 

# lets see first what we would be removing in auto buzz 

buzz_all.1 %>%
  summarise(
    n_total = n(),
    n_na_na = sum(is.na(auto_buzz_count) & is.na(sp)),
    n_na_zero = sum(is.na(sp) & auto_buzz_count == 0, na.rm = TRUE)
  )

# i will remove the nas and zeros for sp and autobuzzcouts.
buzz_all.1 <- buzz_all.1 %>%
    filter(!(is.na(sp) & auto_buzz_count == 0))

# now we remove the instances where sp is na and autobuzz is na. 
buzz_all.1 <- buzz_all.1 %>%
  filter(!(is.na(sp) & is.na(auto_buzz_count)))


# rules for simplifying multi IDs -----------------------------------------


# we have tags in the sp column that correspond to instances where the software id several calls as mixed species. I will first count those and then I will assing the first option. For example if it is myvo/myci/myca then I will assign it as mysp. If it is lano/epfu then I will assign it as lof.any other combination it will be assinged as the first species ID for example coto/epfu/myth/anpa will be assigned as coto.

buzz_all.1 %>%
  summarise(
    n_total = n(),
    n_multi_id = sum(str_detect(sp, "/"), na.rm = TRUE),
    prop_multi_id = n_multi_id / n_total
  )

# the proportion above indicates there's not that many instances where we have a multi ID. So we can drop them. 

# I will turn the myotis ones into mysp and the lano/epfu into lof
# 

# vector of tags you want treated as Myotis spp.
myotis_tags <- c("mylu", "myvo", "myci", "myca", "myev", "myyu")


buzz_all.1 <- buzz_all.1 %>%
  mutate(
    sp_clean = case_when(
      is.na(sp) ~ NA_character_,
      !str_detect(sp, "/") ~ str_to_lower(sp),
      map_lgl(str_split(str_to_lower(sp), "/"), ~ all(.x %in% myotis_tags)) ~ "mysp",
      str_detect(str_to_lower(sp), "lano") & str_detect(str_to_lower(sp), "epfu") ~ "lof",
      TRUE ~ str_extract(str_to_lower(sp), "^[^/]+")
    )
  )


buzz_all.1 %>%
  count(sp, sp_clean, sort = TRUE)

str(buzz_all.1)

unique(buzz_all.1$sp_clean)

## i looks like the rules for multi species work so we replace sp with sp_clean
buzz_all.1 <- buzz_all.1 %>%
  mutate(sp = sp_clean) %>%
  select(-sp_clean)


# sp standardized
# some species names are in four letter code and others in 6. I will have them all in 4 letter code.



sp_key <- c(
  "myocil" = "myci",
  "lasnoc" = "lano",
  "myoluc" = "mylu",
  "eptfus" = "epfu",
  "myoevo" = "myev",
  "myovol" = "myvo",
  "antpal" = "anpa",
  "myocal" = "myca",
  "lascin" = "laci",
  "parhes" = "pahe",
  "cortow" = "coto",
  "myoyum" = "myyu",
  "myothy" = "myth"
)



buzz_all.1 <- buzz_all.1 %>%
  mutate(
    sp_clean = recode(sp, !!!sp_key)
  )



buzz_all.1 %>%
  count(sp, sp_clean, sort = TRUE)

unique(buzz_all.1$sp)
# now we create a treatment column based on the site


# Define the experimental light treatment sites
litsites <- c("iron01", "iron03", "iron05", "long01", "long03")

# Create the treatment and binary columns in one step
buzz_all.1 <- buzz_all.1 %>%
  mutate(
    treatmt = ifelse(site %in% litsites, "lit", "dark"),
    trmt_bin = ifelse(treatmt == "lit", 1, 0)
  )
summary(buzz_all.1)

glimpse(buzz_all.1)

# keep certain cols date, time, sp, date_time, site, treatmt, trmt_bin, auto_buzz_count, manual_buzz_count, year, noche, hi_f, lo_f. 

keep<- c("date_time", "time", "sp", "site", "treatmt", "trmt_bin", "auto_buzz_count", "manual_buzz_count", "year", "noche", "hi_f", "lo_f")

buzz_all.1 <- buzz_all.1 %>%
  select(all_of(keep))

# clean the manual buzz counts
# some buzz counts are zeros and I did not manually check all of those. So if there is a zero in the auto_buzz_count column and a species id in the sp column then the manual buzz counts should get a zero. 

buzz_all.1 <- buzz_all.1 %>%
  mutate(
    manual_buzz_count = if_else(
      auto_buzz_count == 0 & !is.na(sp) & is.na(manual_buzz_count),
      0,
      manual_buzz_count
    )
  )




# Now lets work with the speaker data. We only have 2022 and 2023.


# speaker data ------------------------------------------------------------

path_2022 <- "F:/pioneer_2022/bat_speakers/speaker_organized" # folder
path_2023 <- "G:/PioneerLights_2023/speaker"             # folder

files_2022 <- list.files(path_2022, pattern = "\\.txt$", full.names = TRUE)
files_2023 <- list.files(path_2023, pattern = "\\.txt$", full.names = TRUE)

# read and combine files

speaker_2022 <- map_dfr(files_2022,read_buzz_file, year = 2022)
speaker_2023 <- map_dfr(files_2023, read_buzz_file, year = 2023)


# combine all years into one data frame
spkr_all <- bind_rows(speaker_2022, speaker_2023)
names(spkr_all)


# site column 

# we create a site column based on the file name. 

spkr_all <- spkr_all %>%
  mutate(
    site = str_remove(source_file, "\\.txt$"),
    site = str_replace(site, "^viz", "vizc")
  )

sum(is.na(spkr_all$site)) # check if we have any NAs in the site column. We have 0 so we are good.)



# # create the date column
spkr_all <- spkr_all %>%
  mutate(
    date_time = ymd_hms(str_extract(filename, "\\d{8}_\\d{6}"), tz = "America/Denver")
  )

unique(spkr_all$date_time) # check it work 
sum(is.na(spkr_all$date_time)) # check if we have any NAs in the date_time column. We have 0 so we are good.

# we use the monitoring night for the date_s (standardize date column)


spkr_all <- spkr_all %>%
  rename(noche = monitoring_night)

sum(is.na(spkr_all$noche)) # check if we have any NAs in the noche column. we have 8552 nas that might go away as soon as we remove nas. 

# time column 

spkr_all <- spkr_all %>%
  mutate(
    time = format(date_time, "%H:%M:%S")
  )

sum(is.na(spkr_all$date_time)) # check if we have any NAs in the date_time column. We have 0 so we are good.
# create the sp column 

spkr_all <- spkr_all %>%
  mutate(
    sp = case_when(
      !is.na(spp_accp) ~ spp_accp,
      is.na(spp_accp) & !is.na(species_manual_id) ~ species_manual_id,
      is.na(spp_accp) & is.na(species_manual_id) & !is.na(spp) ~ spp,
      TRUE ~ NA_character_
    )
  )


# rule for multi species --------------------------------------------------



# I think I created a rule to eliminate unecessary zeros. We remove rows where the buzz is zero but there was no species ID suggested in our sp columm.

# I want to see what columns don't have a species ID so I go back and check those. 

# check for missing species ID

missing_sp <- spkr_all %>%
  filter(is.na(sp))
# the majority are zeros for buzz counts meaning there was nothing to ID. So we are good. 

# check species column. 

unique(spkr_all$sp)

# first we make all small caps

spkr_all <- spkr_all %>%
  mutate(
    sp = str_to_lower(sp)
  )

# when we have multiple species separated by a / we will make them into mysp (myotis sp) low for lano/epfu

spkr_all %>%
  summarise(
    n_total = n(),
    n_multi_id = sum(str_detect(sp, "/"), na.rm = TRUE),
    prop_multi_id = n_multi_id / n_total
  )

# the proportion above indicates there's not that many instances where we have a multi ID. So we can drop them. 

# I will turn the myotis ones into mysp and the lano/epfu into lof
# below I make multi ID into single tags using the first species proposed as in robomoth. 

# vector of tags you want treated as Myotis spp.
myotis_tags <- c("mylu", "myvo", "myci", "myca", "myev", "myyu")


spkr_all <- spkr_all %>%
  mutate(
    sp_clean = case_when(
      is.na(sp) ~ NA_character_,
      !str_detect(sp, "/") ~ str_to_lower(sp),
      map_lgl(str_split(str_to_lower(sp), "/"), ~ all(.x %in% myotis_tags)) ~ "mysp",
      str_detect(str_to_lower(sp), "lano") & str_detect(str_to_lower(sp), "epfu") ~ "lof",
      TRUE ~ str_extract(str_to_lower(sp), "^[^/]+")
    )
  )


spkr_all %>%
  count(sp, sp_clean, sort = TRUE)

str(spkr_all)

unique(spkr_all$sp_clean)

## i looks like the rules for multi species work so we replace sp with sp_clean
spkr_all <- spkr_all %>%
  mutate(sp = sp_clean) %>%
  select(-sp_clean)

unique(spkr_all$sp)


# now we add the treatment column.


# Define the experimental light treatment sites
litsites <- c("iron01", "iron03", "iron05", "long01", "long03")

# Create the treatment and binary columns in one step
spkr_all <- spkr_all %>%
  mutate(
    treatmt = ifelse(site %in% litsites, "lit", "dark"),
    trmt_bin = ifelse(treatmt == "lit", 1, 0)
  )
summary(spkr_all)

# clean the manual buzz counts
# some buzz counts are zeros and I did not manually chechk all of those. So if there is a zero in the auto_buzz_count column and a species id in the sp column then the manual buzz counts should get a zero. 

spkr_all <- spkr_all %>%
  mutate(
    manual_buzz_count = if_else(
      auto_buzz_count == "0" & !is.na(sp_clean),
      "0",
      manual_buzz_count
    )
  )



# make buzz counts numeric

spkr_all <- spkr_all %>%
  mutate(
    auto_buzz_count = as.numeric(auto_buzz_count),
    manual_buzz_count = as.numeric(manual_buzz_count),
  ) 

str(spkr_all)

# lets remove rows where with NAs for auto buzz count and sp column because these are not meaningful biologically. 

# lets see first what we would be removing in auto buzz 

spkr_all %>%
  summarise(
    n_total = n(),
    n_na_na = sum(is.na(auto_buzz_count) & is.na(sp)),
    n_na_zero = sum(is.na(sp) & auto_buzz_count == 0, na.rm = TRUE)
  )

# i will remove the nas and zeros for sp and autobuzzcouts about 
spkr_all <- spkr_all %>%
  filter(!(is.na(sp) & auto_buzz_count == 0))


# # now we filter out those rows that have NAs in manual buzz count.NO! we want to filter where there is a zero in auto buzz count and no sp id in the sp column. 
# 
# spkr_all.1 <- spkr_all %>%
#   filter(!is.na(manual_buzz_count))


# now there are still some rows where noche is NA.
sum(is.na(spkr_all$noche)) 

# several of those have inf in call_sec column. I would say we need to remove those becasue it indicates some kind of error with the recording. so we remove rows where both calls_sec is inf and noche is NA. 

spkr_all <- spkr_all %>%
  filter(!(calls_sec == "Inf" & is.na(noche)))


# noche as date
spkr_all <- spkr_all %>%
  mutate(
    noche = ymd(noche, tz = "America/Denver")
  )

sum(is.na(spkr_all$noche)) 
 
# final checks. 

summary(spkr_all)

# keep just necessary columns

keep<- c("date_time", "time", "sp", "site", "treatmt", "trmt_bin", "auto_buzz_count", "manual_buzz_count", "year", "noche", "hi_f", "lo_f")

spkr_all <- spkr_all %>%
  select(all_of(keep))


unique(spkr_all$site)



# explore data ------------------------------------------------------------


# I want a bar graph of species in the x axis and buzz counts by treatment in the y axis.
# first we need to make the manual buzz count numeric.

spkr_all <- spkr_all %>%
  mutate(
    manual_buzz_count = as.integer(manual_buzz_count)
  )

# now ggplot to graph it 

sp_summary <- spkr_all %>%
  group_by(sp, treatmt) %>%
  summarise(
    total_buzz = sum(as.numeric(manual_buzz_count), na.rm = TRUE),
    .groups = "drop"
  )

ggplot(sp_summary, aes(x = sp, y = total_buzz, fill = treatmt)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Species", y = "Total Manual Buzz Count spkrs", fill = "Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# what about for robomoths in the buzz_all.1 data frame.

# first we need to make the manual buzz count numeric.

buzz_all.1 <- buzz_all.1 %>%
  mutate(
    manual_buzz_count = as.integer(manual_buzz_count)
  )


sp_summary <- buzz_all.1 %>%
  group_by(sp, treatmt) %>%
  summarise(
    total_buzz = sum(as.numeric(manual_buzz_count), na.rm = TRUE),
    .groups = "drop"
  )

ggplot(sp_summary, aes(x = sp, y = total_buzz, fill = treatmt)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Species", y = "Total Manual Buzz Count spkrs", fill = "Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






















# junk --------------------------------------------------------------------


buzz_all.1 <- buzz_all.1 %>%
  mutate(
    date_s = if_else(hour(date_time) < 9, date(date_time) - days(1), date(date_time))
  )

all(buzz_all.1$date_s == buzz_all.1$monitoring_night, na.rm = TRUE)

buzz_all.1 %>%
  filter(date_s != monitoring_night | 
           (is.na(date_s) != is.na(monitoring_night))) %>%
  select(filename, date_s, monitoring_night)


files_2022 <- list.files(path_2022, pattern = "\\.txt$", full.names = TRUE)
files_2023 <- list.files(path_2023, pattern = "\\.txt$", full.names = TRUE)

# read and combine files

robomoth_2022 <- map_dfr(files_2022, read_buzz_file, year = 2022)
robomoth_2023 <- map_dfr(files_2023, read_buzz_file, year = 2023)



#remove nas for the column manual buzz count, but first we create a sp column that combines the spp_accp, the manual id and the spp column. We do this because some observations have buzz count but not manual id as I was not always including them so we build a new one combining the spp_accp, the manual id and the spp column.
#this is because we will focus only on the ones that we have manually looked around 11k records.

# buzz_all.1 <- buzz_all %>%
#   filter(!is.na(manual_buzz_count))
#   
#   
#   # # create sp column
# we create this column because some observatins have buzz count but not manual id as I was not always including them so we build a new one combining the spp_accp, the manual id and the spp column.

spkr_all <- spkr_all %>%
  mutate(
    sp = case_when(
      !is.na(spp_accp) ~ spp_accp,
      is.na(spp_accp) & !is.na(species_manual_id) ~ species_manual_id,
      is.na(spp_accp) & is.na(species_manual_id) & !is.na(spp) ~ spp,
      TRUE ~ NA_character_
    )
  )
