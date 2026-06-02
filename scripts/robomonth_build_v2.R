
#=============================================================================
# Script Name:    robomoth_build_2.R
# Purpose:        build rbomoth and speaker database
#
# Version:        2
#
# Author:         Carlos Linares
# Collaborator:   
# Created:        2026-03-17
# Contact:        carlosgarcialina@u.boisestate.edu
#
# Notes:
#   - original data stored in the bk_linares and PhD_bk hard drives
#   - version 2 has script added to change  the sp column to a new one that combines the spp_accp, the manual id and the spp column. We do this because some observations have buzz count but not manual id as I was not always including them so we build a new one combining the spp_accp, the manual id and the spp column.
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

path_2021 <- "F:/pioneer_2021/sm4/sb_buzz_2021_vetted_v4.txt"   # single file
path_2022 <- "G:/PioneerLights_2022/part2/robomoth/buzz_out" # folder
path_2023 <- "F:/pioneer_2023/robomoth_2023_all/buzz_out/"             # folder

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

robomoth_2021 <- read_buzz_file(path_2021, 2021)



# read 2022 and 2023


files_2022 <- list.files(path_2022, pattern = "\\.txt$", full.names = TRUE)
files_2023 <- list.files(path_2023, pattern = "\\.txt$", full.names = TRUE)

# read and combine files

robomoth_2022 <- map_dfr(files_2022, read_buzz_file, year = 2022)
robomoth_2023 <- map_dfr(files_2023, read_buzz_file, year = 2023)



#check col names --------------------------------------------------------

# Combine all robomoth buzz files into one file list
# 2021 is a single file, 2022 and 2023 are folders
buzz_files_all <- c(
  path_2021,
  files_2022,
  files_2023
)

# Function to read only the column names from one file
check_kpro_column <- function(file) {
  
  # Read only the header / first row structure
  # n_max = 0 means we read no data rows, only column names
  col_names <- read_tsv(
    file,
    n_max = 0,
    show_col_types = FALSE
  ) %>%
    clean_names() %>%
    names()
  
  tibble(
    source_file = basename(file),
    full_path = file,
    has_kpro_autoid = "kpro_autoid" %in% col_names,
    has_wa_kaleidoscope_auto_id = "wa_kaleidoscope_auto_id" %in% col_names,
    n_columns = length(col_names),
    column_names = paste(col_names, collapse = "; ")
  )
}

# Apply the check to every file
kpro_column_check <- map_dfr(buzz_files_all, check_kpro_column)

# View summary
kpro_column_check %>%
  count(has_kpro_autoid, has_wa_kaleidoscope_auto_id)

# after checking 2021, 2022, and 2023 data I see only 2022 has been run through kpro and thus we can use the autoID col for all the years. 

#combine robomoth --------------------------------------------------------

# combine all years into one data frame
buzz_all <- bind_rows(robomoth_2021, robomoth_2022, robomoth_2023)

str(buzz_all)



# rename the kaleidoscope column to kpro_autoid

buzz_all<-rename(buzz_all, kpro_autoid = wa_kaleidoscope_auto_id) # we rename the kaleidoscope column to kpro_autoid to avoid confusion with the spp column that also has kaleidoscope IDs.)

glimpse(buzz_all)



# sp column ---------------------------------------------------------------


# create sp column
# we create this column because some observations have buzz count but not manual id we want to recover those IDs to not loss data. We build a new one combining the spp_accp, the manual id and the spp column and the 1st, 2nd, 3nd columns 



buzz_all.1 <- buzz_all %>%
  mutate(
    sp = case_when(
      !is.na(species_manual_id) ~ species_manual_id, # this gives priority to manual ID.
      is.na(species_manual_id) & !is.na(spp_accp) ~ spp_accp,
      is.na(species_manual_id) & is.na(spp_accp) & !is.na(spp) ~ spp,
      is.na(species_manual_id) & is.na(spp_accp) & is.na(spp) & !is.na(kpro_autoid) ~ kpro_autoid,
      TRUE ~ NA_character_ # this makes sure sp is a character column. 
    )
  )

sum(is.na(buzz_all.1$sp)) # we have 26 223 NAs in the sp column.

# site column 

# we create a site column based on the file name. 

buzz_all.1 <- buzz_all.1 %>%
  mutate(
    site = str_to_lower(str_extract(filename, "^[A-Za-z]{3,4}[0-9]{2}"))
  )


# check sites 
unique(buzz_all.1$site)
sum(is.na(buzz_all.1$site)) # check if we have any NAs in the site column. We have 0 so we are good.but the names need work. 


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

sum(is.na(buzz_all.1$date_time)) # check if we have any NAs in the date_time column. We have 0 so we are good.)

# we use the monitoring night for the date_s (standardize date column)

# after checking the difference I realize that my method creates some observations where the date is not corrected properly specially if close or around 12 am. So I will just use the monitoring night column as the date column for the analysis

buzz_all.1 <- buzz_all.1 %>%
  rename(noche = monitoring_night)

sum(is.na(buzz_all.1$noche)) # we have 146 missing monitoring nights. many have no data so I will recheck after clean up 

# time column 

buzz_all.1 <- buzz_all.1 %>%
  mutate(
    time = format(date_time, "%H:%M:%S")
  )

sum(is.na(buzz_all.1$date_time)) # check if we have any NAs in the date_time column. We have 0 so we are good.



# checking for NAs in the auto buzz. this would indicate that that file has not been analyzed with the buzz tools from Sonobat. therw are some missing might be from the old table. checking 
#

missing_auto <- buzz_all.1 %>%
  filter(
    is.na(auto_buzz_count) |
      auto_buzz_count == "" |
      auto_buzz_count == "NA"
  ) %>%
  select(filename, auto_buzz_count, kpro_autoid, source_file )


# several vizc03 2022 files are mising autobuzz we should check. 

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
  filter(is.na(sp)) %>% 
  select(filename, sp, auto_buzz_count, manual_buzz_count, source_file)

# the majority are zeros for buzz counts meaning there was nothing to ID. We might eliminate these. In addition it seems these also have inf for call/sec column possibly indicating these are noise. se we might remove all inf rows for call/sec and NAs for sp. 

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

# lets see first what we would be removing in auto buzz 

buzz_all.1 %>%
  summarise(
    n_total = n(),
    n_na_na = sum(is.na(auto_buzz_count) & is.na(sp)),
    n_na_zero = sum(is.na(sp) & auto_buzz_count == 0, na.rm = TRUE)
  )

# # i will remove the nas and zeros for sp and autobuzzcouts.
# buzz_all.1 <- buzz_all.1 %>%
#     filter(!(is.na(sp) & auto_buzz_count == 0))
# 
# # now we remove the instances where sp is na and autobuzz is na. 
# buzz_all.1 <- buzz_all.1 %>%
#   filter(!(is.na(sp) & is.na(auto_buzz_count)))

# remove rows with noise, noid and inf for the call/sec column  
# the code below removes rows where sp=NA and calls_sec is inf. 

buzz_all.1 <- buzz_all.1 %>%
  filter(
    !(is.na(sp) & calls_sec == "Inf"),
    !str_to_lower(sp) %in% c("noise", "noid") 
  )

# we are going to remove NAs for sp about 116 rows
# these rows show as having buzz counts but no sp. these might need to be check in the future.  

sum(is.na(buzz_all.1$sp))

buzz_all.1 <- buzz_all.1 %>%
  filter(!is.na(sp))

# rules for simplifying multi IDs -----------------------------------------


# we have tags in the sp column that correspond to instances where the software id several calls as mixed species. I will first count those and then I will assign the first option. For example if it is myvo/myci/myca then I will assign it as mysp. If it is lano/epfu then I will assign it as lof.any other combination it will be assigned as the first species ID for example coto/epfu/myth/anpa will be assigned as coto.

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
      # if there's no / lower case it.
      !str_detect(sp, "/") ~ str_to_lower(sp),
      # if all the tags in the multi ID are myotis then we assign it as mysp
      map_lgl(str_split(str_to_lower(sp), "/"), ~ all(.x %in% myotis_tags)) ~ "mysp",
      # if the multi ID contains lano and epfu then we assign it as lof
      str_detect(str_to_lower(sp), "lano") &
        str_detect(str_to_lower(sp), "epfu") ~ "lof",
      
      TRUE ~ str_extract(str_to_lower(sp), "^[^/]+")
    )
  )

sum(is.na(buzz_all.1$sp_clean)) # we have 0 NAs in the sp clean column. 

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

unique(buzz_all.1$sp)


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
  "myothy" = "myth",
  "hifrag" = "hif",
  "lofrag" = "lof"
  )



buzz_all.1 <- buzz_all.1 %>%
  mutate(
    sp = recode(sp, !!!sp_key)
  )



buzz_all.1 %>%
  count(sp, sort = TRUE)

unique(buzz_all.1$sp)


# recode myth, myca, myyu as my sp while we code euma as lof. 

buzz_all.1 <- buzz_all.1 %>%
  mutate(
    sp = case_when(
      sp %in% c("myth", "myca", "myyu") ~ "mysp",
      sp == "euma" ~ "lof",
      TRUE ~ sp
    )
  )



# now we correct the counts by giving priority to the manual buzz counts but if these are missing then we input the auto buzz counts. 

buzz_all.1 <- buzz_all.1 %>%
  mutate(
    c_buzz = if_else( # corrected buzz counts 
      is.na(manual_buzz_count) & !is.na(auto_buzz_count),
      auto_buzz_count,
      manual_buzz_count
    )
  )


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

# keep certain cols date, time, sp, date_time, site, treatmt, trmt_bin, auto_buzz_count, manual_buzz_count, year, noche, hi_f, lo_f and file.

keep <- c(
  "filename",
  "source_file",
  "date_time",
  "time",
  "sp",
  "site",
  "treatmt",
  "trmt_bin",
  "auto_buzz_count",
  "manual_buzz_count",
  "c_buzz",
  "year",
  "noche",
  "hi_f",
  "lo_f"
)

buzz_all.1<- buzz_all.1 %>%
  select(all_of(keep))

glimpse(buzz_all.1)


# I need to check what files are duplicated in the buzz_all.1 data frame because I will be merging it with the amplitude data frame and I want to make sure I am not merging multiple rows with the same file name.

# cehck for duplicate files 
buzz_dup <- buzz_all.1 %>%
  count(filename, name = "n_buzz_rows") %>%
  filter(n_buzz_rows > 1) %>%
  arrange(desc(n_buzz_rows))

# make all letters small caps in the file name and remove any leading or trailing white space to check for duplicates again.
buzz_all.1 <- buzz_all.1 %>%
  mutate(
    filename_norm = str_to_lower(str_trim(basename(filename)))
  )

buzz_dup_norm <- buzz_all.1 %>%
  count(filename_norm, name = "n_buzz_rows") %>%
  filter(n_buzz_rows > 1) %>%
  arrange(desc(n_buzz_rows))

buzz_dup_norm


buzz_exact_dups <- buzz_all.1 %>%
  group_by(across(everything())) %>%
  filter(n() > 1) %>%
  ungroup()

buzz_exact_dups

# the buzz files don't have any duplicated files. 
# # -------------------------------------------------------------------------
# Handle duplicated audio files
# Rule:
# 1. If same filename + same species appears more than once, keep one row.
# 2. If same filename has different species, keep all species rows.
# -------------------------------------------------------------------------

buzz_all.1 <- buzz_all.1 %>%
  mutate(
    filename_norm = str_to_lower(str_trim(basename(filename)))
  )

buzz_all.2 <- buzz_all.1 %>%
  distinct(
    filename_norm,
    sp,
    .keep_all = TRUE
  )

# amplitude  --------------------------------------------------------------

# paths to folders with amplitude csv files
amp_path_2022_robo <- "G:/PioneerLights_2022/part2/robomoth/buzz_out"
amp_path_2023_robo <- "F:/pioneer_2023/robomoth_2023_all/buzz_out/"

# list amplitude files from both folders
amp_files_2022_robo <- list.files(
  path = amp_path_2022_robo,
  pattern = "_Amp\\.csv$",
  full.names = TRUE
)

amp_files_2023_robo <- list.files(
  path = amp_path_2023_robo,
  pattern = "_Amp\\.csv$",
  full.names = TRUE
)

# combine both file lists
amp_files_robo <- c(
  amp_files_2022_robo,
  amp_files_2023_robo
)

# check how many files were found and see them
length(amp_files_robo)
amp_files_robo


# Check which amplitude files have RMS column
# -------------------------------------------------------------------------

check_amp_columns <- function(file) {
  
  col_names <- read_csv(
    file,
    n_max = 0,
    show_col_types = FALSE
  ) %>%
    clean_names() %>%
    names()
  
  tibble(
    source_amp_file = basename(file),
    full_path = file,
    has_filename = "filename" %in% col_names,
    has_max_d_bfs = "max_d_bfs" %in% col_names,
    has_avg_d_bfs = "avg_d_bfs" %in% col_names,
    has_rms_value = "rms_value" %in% col_names,
    has_rms = "rms" %in% col_names,
    has_rms_val = "rms_val" %in% col_names,
    has_any_rms = any(col_names %in% c("rms_value", "rms", "rms_val")),
    column_names = paste(col_names, collapse = "; ")
  )
}

amp_column_check <- map_dfr(amp_files_robo, check_amp_columns)

# Summary of RMS availability
amp_column_check %>%
  count(has_any_rms)

# Files that DO have RMS
amp_column_check %>%
  filter(has_any_rms) %>%
  select(source_amp_file, has_rms_value, has_rms, has_rms_val)

# Files that DO NOT have RMS
amp_files_missing_rms <- amp_column_check %>%
  filter(!has_any_rms)

amp_files_missing_rms %>%
  select(source_amp_file, full_path, column_names) %>%
  print(n = Inf)




# now we read them and combine them 

amp_robo_all <- amp_files_robo %>%
  map_dfr(function(file) {
    
    read_csv(
      file,
      show_col_types = FALSE,
      col_types = cols(.default = col_character())
    ) %>%
      clean_names() %>%
      select(
        filename,
        max_d_bfs,
        avg_d_bfs
      ) %>%
      mutate(
        source_amp_file = basename(file),
        amp_folder = dirname(file),
        # Extract only the site code, even if filename has year
        site = str_to_lower(str_extract(source_amp_file, "^[A-Za-z]{3,4}[0-9]{2}"))
      )
  }) %>%
  mutate(
    max_dbfs = as.numeric(max_d_bfs),
    avg_dbfs = as.numeric(avg_d_bfs),
    filename_norm = str_to_lower(basename(filename)),
    
    # extract year from the filename, e.g. iron01-20220707_230159.wav
    year = as.integer(str_extract(filename_norm, "20\\d{2}"))
  ) %>%
  select(
    filename,
    filename_norm,
    max_dbfs,
    avg_dbfs,
    site,
    year,
    source_amp_file,
    amp_folder
  )


# let's check the output

glimpse(amp_robo_all)
unique(amp_robo_all$site)

amp_robo_all %>%
  count(year, site) #There's a mistake with the year on iron 04 22 or 23 files. 

# check misssing values of dbfs. 
amp_robo_all %>%
  summarise(
    n_rows = n(),
    n_files = n_distinct(filename_norm),
    n_missing_avg_dbfs = sum(is.na(avg_dbfs)),
    n_missing_max_dbfs = sum(is.na(max_dbfs))
  )

# check duplicated rows
amp_dups <- amp_robo_all %>%
  mutate(
    filename_norm = str_to_lower(str_trim(basename(filename)))
  ) %>%
  count(filename_norm, name = "n_amp_rows") %>%
  filter(n_amp_rows > 1) %>%
  arrange(desc(n_amp_rows))

amp_dups

# now I want to merge buzz_all.1 with the amplitude data frame. I will merge by filename 

# I have more than 29K observations that I need to add amplitude data by merging it with the amp_robo_all data frame. I will do a left join to keep all the observations in buzz_all.1 and add the amplitude data where available.

head(buzz_all.2$filename) #buzz.all.2 has the duplicated files removed 
head(amp_robo_all$filename)

# filenames differ so we need to standardize to merge properly.

buzz_all.2 <- buzz_all.2 %>%
  mutate(
    filename_norm = filename %>%
      str_to_lower() %>%
      str_replace("_", "-") %>%              # convert IRON03_20210729 to iron03-20210729
      str_remove("-[a-z]+(?=\\.wav$)") %>%   # remove species suffix like -myvo.wav
      str_trim()
  )
glimpse(buzz_all.2)

amp_robo_all <- amp_robo_all %>%
  mutate(
    filename_norm = filename %>%
      str_to_lower() %>%
      str_replace("_", "-") %>%
      str_trim()
  )

glimpse(amp_robo_all)

# remove nas and duplicates from the amplitude data frame before merging.
# there is no duplicates but just to be sure
amp_robo_clean <- amp_robo_all %>%
  group_by(filename_norm) %>%
  summarise(
    max_dbfs = mean(max_dbfs, na.rm = TRUE),
    avg_dbfs = mean(avg_dbfs, na.rm = TRUE),
    n_amp_rows = n(),
    amp_sources = paste(unique(source_amp_file), collapse = "; "),
    .groups = "drop"
  )

# merge buzz data and amplitude (dbfs)
buzz_all.2_amp <- buzz_all.2 %>%
  left_join(
    amp_robo_clean,
    by = "filename_norm",
    relationship = "many-to-one"
  )


# how many matched 
buzz_all.2_amp %>%
  summarise(
    n_rows = n(),
    n_with_amp = sum(!is.na(avg_dbfs)),
    n_missing_amp = sum(is.na(avg_dbfs)),
    prop_with_amp = mean(!is.na(avg_dbfs))
  )

# missing dbfs
# there are several missing amp data but we will check them later. 
missing<-buzz_all.2_amp %>%
  filter(year %in% c(2022, 2023), is.na(avg_dbfs)) %>%
  select(filename, filename_norm, year, site, sp, c_buzz) %>%
  arrange(year, site, filename_norm) 

# the file buzz_all.2_amp is the one that should be exported and analyzed. 
# graph the number of observations at each dbfs bin and sp dat for 2022 and 2023

ggplot(buzz_all.2_amp, aes(x = avg_dbfs)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  labs(title = "Distribution of Average dBFS", x = "Average dBFS", y = "Count") +
  theme_minimal()


buzz_all.2_amp %>%
  filter(
    year %in% c(2022, 2023),
    !is.na(avg_dbfs),
    !is.na(sp)
  ) %>%
  ggplot(aes(x = avg_dbfs)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  facet_wrap(~ sp, scales = "free_y") +
  labs(
    title = "Distribution of Average dBFS by Species",
    x = "Average dB",
    y = "Count"
  ) +
  theme_minimal()


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
    site = str_replace(site, "^viz(?!c)", "vizc")
  )

sum(is.na(spkr_all$site)) # check if we have any NAs in the site column. We have 0 so we are good.)

unique(spkr_all$site) # check sites


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


# sp column ---------------------------------------------------------------

# create the sp column 

spkr_all <- spkr_all %>%
  mutate(
    sp = case_when(
      !is.na(species_manual_id) ~ species_manual_id, # this gives priority to manual ID.
      is.na(species_manual_id) & !is.na(spp_accp) ~ spp_accp,
      is.na(species_manual_id) & is.na(spp_accp) & !is.na(spp) ~ spp,
      TRUE ~ NA_character_
    )
  )

sum(is.na(spkr_all$sp)) # we have 91192 NAs



# rule for multi species --------------------------------------------------

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

## i looks like the rules for multi species work so we replace sp_clean with sp
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



# -------------------------------------------------------------------------

# make buzz counts numeric

spkr_all <- spkr_all %>%
  mutate(
    auto_buzz_count = as.numeric(auto_buzz_count),
    manual_buzz_count = as.numeric(manual_buzz_count),
  ) 

str(spkr_all)

# lets see first what we would be removing in auto buzz 

spkr_all %>%
  summarise(
    n_total = n(),
    n_na_na = sum(is.na(auto_buzz_count) & is.na(sp)),
    n_na_zero = sum(is.na(sp) & auto_buzz_count == 0, na.rm = TRUE)
  )

# remove rows with noise, noid and inf for the call/sec column  

spkr_all <- spkr_all %>%
  filter(
    !(is.na(sp) & calls_sec == "Inf"),
    !str_to_lower(sp) %in% c("noise", "noid") 
  )

# we are going to remove NAs for sp about 381 rows 

sum(is.na(spkr_all$sp))

spkr_all <- spkr_all %>%
  filter(!is.na(sp))


# clean the manual buzz counts
# some buzz counts are zeros and I did not manually check all of those. So if there is a zero in the auto_buzz_count column and a species id in the sp column then the manual buzz counts should get a zero. 

spkr_all <- spkr_all %>%
  mutate(
    manual_buzz_count = if_else(
      auto_buzz_count == 0 & !is.na(sp),
      0,
      manual_buzz_count
    )
  )

spkr_all <- spkr_all %>%
  mutate(
    c_buzz = if_else( # corrected buzz counts 
      is.na(manual_buzz_count) & !is.na(auto_buzz_count),
      auto_buzz_count,
      manual_buzz_count
    )
  )

summary(spkr_all)




# noche as date
spkr_all <- spkr_all %>%
  mutate(
    noche = ymd(noche, tz = "America/Denver")
  )

sum(is.na(spkr_all$noche)) 
 
# checks. 

summary(spkr_all)

# keep just necessary columns

keep<- c("date_time", "time", "sp", "site", "treatmt", "trmt_bin", "auto_buzz_count", "manual_buzz_count", "c_buzz","year", "noche", "hi_f", "lo_f")

spkr_all <- spkr_all %>%
  select(all_of(keep))


unique(spkr_all$site)
summary(spkr_all)


# explore data ------------------------------------------------------------


# I want a bar graph of species in the x axis and buzz counts by treatment in the y axis.
# first we need to make the manual buzz count numeric.


# now ggplot to graph it 

sp_summary <- spkr_all %>%
  group_by(sp, treatmt,year) %>%
  summarise(
    total_buzz = sum(as.numeric(c_buzz), na.rm = TRUE),
    .groups = "drop"
  )

ggplot(sp_summary, aes(x = sp, y = total_buzz, fill = treatmt)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(~ year) +
  labs(x = "Species", y = "Total Manual Buzz Count spkrs", fill = "Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# what about for robomoths in the buzz_all.1 data frame.

# first we need to make the manual buzz count numeric.


sp_summary <- buzz_all.1 %>%
  group_by(sp, treatmt, year) %>%
  summarise(
    total_buzz = sum(as.numeric(manual_buzz_count), na.rm = TRUE),
    .groups = "drop"
  )

ggplot(sp_summary, aes(x = sp, y = total_buzz, fill = treatmt)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(~ year) +
  labs(x = "Species", y = "Total Manual Buzz Count spkrs", fill = "Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# after looking at the graphs I chose to analyse both 2021 and 2022 for robomoth and speakers as the data is consisten

# What about the zeros? is there more calls with no buzz in lit vs dark sites?



buzz_all.1 <- buzz_all.1 %>%
  mutate(
    zero_buzz = c_buzz == 0
  )

prop_zero <- buzz_all.1 %>%
  group_by(treatmt) %>%
  summarise(
    prop_zero = mean(zero_buzz, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

ggplot(prop_zero, aes(x = treatmt, y = prop_zero, fill = treatmt)) +
  geom_col() +
  labs(
    x = "Treatment",
    y = "Proportion of zero feeding events",
    title = "Zero feeding events (c_buzz = 0) by treatment"
  ) +
  theme_minimal()




count_zero <- buzz_all.1 %>%
  filter(zero_buzz) %>%
  count(treatmt)

ggplot(count_zero, aes(x = treatmt, y = n, fill = treatmt)) +
  geom_col() +
  labs(
    x = "Treatment",
    y = "Number of zero feeding events"
  ) +
  theme_minimal()


buzz_all.1 %>%
  count(treatmt, zero_buzz) %>%
  ggplot(aes(x = treatmt, y = n, fill = zero_buzz)) +
  geom_col(position = "fill") +
  labs(
    x = "Treatment",
    y = "Proportion",
    fill = "Zero buzz",
    title = "Proportion of zero vs non-zero feeding events"
  ) +
  theme_minimal()


# I am a little counfuse about how to analyze this section. should I model the proportions of zeros too? in the model 


# dir.create("data_for_analysis/prep_for_glm", showWarnings = FALSE) # just run if the dir is abscent

write.csv(buzz_all.1, file = 'data_for_analysis/robomoth_build/bat_robomoth.csv', row.names = F) # raw combine data 
write.csv(spkr_all, file = 'data_for_analysis/robomoth_build/spkr_all.csv', row.names = F) # raw combine data


# Create a README file with information about the script
readme_content <- "Carlos Linares 3/22/2026 
this folde contains the raw combined data for the robomoth and speaker data. The data is cleaned and ready for analysis The data are as follows:

buzz_all.1: this is the combined data frame for the robomoth data. It contains data from 2021, 2022 and 2023. The data has been cleaned and ready for analysis. 

the columns are as follows
date_time: the date and time of the recording
time: the time of the recording
noche: the monitoring night of the recording
sp: the species ID of the recording composed of sp_accp, manual id and spp columns. We created this column to recover as much species ID as possible.
site: the site of the recording
treatmt: the treatment of the site (lit or dark)
c_buzz: the corrected buzz count. We give priority to the manual buzz count but if this is missing then we use the auto buzz count
zero_buzz: created for graph the zero proportion probably unneccessary for analysis.

spkr_all: this is the combined data frame for the speaker data. It contains data from 2022 and 2023. The data has been cleaned and ready for analysis.

the columns are as follows:
date_time: the date and time of the recording
time: the time of the recording
noche: the monitoring night of the recording
sp: the species ID of the recording composed of sp_accp, manual id and spp columns
site: the site of the recording
treatmt: the treatment of the site (lit or dark)
c_buzz: the corrected buzz count. We give priority to the manual, same as robomoth, if this is missing then we use the auto buzz count

"
# Write the README content to a file
writeLines(readme_content, "data_for_analysis/robomoth_build/README.txt")








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
