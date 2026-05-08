# 
# Script name: Prep_for_glmm
# 
# Purpose of script: Combine the 2021, 2022, 2023, data sets and get them ready to go throuhg glmm_v1
# 
# Author: Carlos Linares
# 
# Date Created: 07/29/2024
# 
# Email: carlosgarcialina@u.boisestate.edu
# 
# ---------------------------
#   
#   Notes: we also should check we have the same temporal window for all variables.  
# sessionInfo() at end of script


# inputs ------------------------------------------------------------------
# -data_for_analysis/2021_kpro_raw/bats2021_kpro_v1.csv
# -data_for_analysis/2022_kpro_raw/bat2022_kpr.csv
# -data_for_analysis/2023_kpro_raw/bat2023_kpr.csv

# outputs ----------------------

# bat_combined.csv - process data no counts
# bm.csv - counts of bat calls by day from 2021 to 2023 all sites
# bm.miller.day - number of minutes of activity by day  for 2021-2023 data all sites

# this should be a database ready to analyze with the glmm_v1 script. 

# libraries

library(tidyverse)
library(beepr)
library(lubridate)
library(data.table)
library(janitor)


# kpro_data --------------------------------------------------------------

# load
# 
load("working_env/prep_for_glm.RData")


kpro_2021_bat <- fread(file = 'data_for_analysis/2021_kpro_raw/bats2021_kpro_v1.csv', header = T)
kpro_2022_bat <- fread(file = 'data_for_analysis/2022_kpro_raw/bat2022_kpr.csv', header = T)
kpro_2023_bat <- fread(file = 'data_for_analysis/2023_kpro_raw/bat2023_kpr.csv', header = T)


bat_combined <- bind_rows(kpro_2021_bat, kpro_2022_bat, kpro_2023_bat)
bat_combined<- clean_names(bat_combined)
head(bat_combined)
summary(bat_combined) # check the output
names(bat_combined)
keep <- c(
  "id",
  "indir",
  "outdir",
  "folder",
  "in_file",
  "duration",
  "date",
  "time",
  "hour",
  "date_12",
  "time_12",
  "auto_id",
  "pulses",
  "alternate_1",
  "alternate_2",
  "n",
  "fc",
  "sc",
  "dur",
  "fmax",
  "fmin",
  "fmean",
  "tbc",
  "fk",
  "tk",
  "s1",
  "tc",
  "qual",
  "files",
  "manual_id"
) # cols to keep

bat_combined <- bat_combined %>% select(all_of(keep)) # keeps variables of interest


# site 

bat_combined$site<-str_extract(bat_combined$OUTDIR, "[A-Za-z]{3,4}\\d{2}")

unique(bat_combined$site) # site labels need correction

# here we fix some problems with the sites misspellings

bat_combined$site = ifelse(bat_combined$site %in% "Iron02","iron02", bat_combined$site)
bat_combined$site = ifelse(bat_combined$site %in% "viz01","vizc01", bat_combined$site)
bat_combined$site = ifelse(bat_combined$site %in% "viz02","vizc02", bat_combined$site)
bat_combined$site = ifelse(bat_combined$site %in% "viz03","vizc03", bat_combined$site)
bat_combined$site = ifelse(bat_combined$site %in% "viz04","vizc04", bat_combined$site)

unique(bat_combined$site) # site labels


# date 

bat_combined$date<-lubridate::ymd(bat_combined$date)
sum(is.na(bat_combined$date)) # check for NAs. 

# noche/night
# this variable assigns the same date to calls that belong to a single night. If a call happened on July 6 at 3 am it assigns it to July 5. 
# update on this. I realized that sometimes when I do this the correct date is not input. However, the date.12 column in Kaleidocope does the right thing. Sow we updated this sectio to use the date.12 insetead of the following lines. 

# bat_combined$noche <-
#   if_else(bat_combined$HOUR < 9, # if it is less than 9 put the date of the previous day
#           true =  (date(bat_combined$DATE) - ddays(1)),
#           false = date(bat_combined$DATE))
# 
# sum(is.na(bat_combined$noche)) # check for NAs. 

bat_combined$noche<- bat_combined$date_12
sum(is.na(bat_combined$noche)) # chceck for NAs


# date time col. 
bat_combined$date<-as.character(bat_combined$date)
bat_combined$time<-as.character(bat_combined$time)

datetime<-paste(bat_combined$date, bat_combined$time)#merge date and time 

bat_combined$datetime<- datetime

bat_combined <- bat_combined %>%
  mutate(
    datetime = ymd_hms(datetime, tz = "America/Denver"), # Parse as POSIXct
    # Add 10 seconds to midnight times
    datetime = if_else(
      format(datetime, "%H:%M:%S") == "00:00:00", 
      datetime + seconds(10), 
      datetime
    )
  )

# View result

midnight_times <- bat_combined$datetime %>%
  as_tibble() %>%
  filter(hour(value) == 0 & minute(value) == 0 & second(value) == 10)

# View the rows that match
midnight_rows


#year

bat_combined$yr<-year(bat_combined$date)
unique(bat_combined$yr) #we check the three years are present.


# treatment column 

litsites<-c("iron01","iron03","iron05","long01","long03")


bat_combined$treatmt <- ifelse(bat_combined$site %in% litsites , "lit", "dark") # this makes a treatment variable.

bat_combined$trmt_bin<- ifelse(bat_combined$treatmt== "lit", 1, 0)


bat_combined$jday<-lubridate::yday(bat_combined$noche) # julian day


summary(bat_combined)


# species rules  ----------------------------------------------------------

# now we are going to create a series of rules to make sure we don't have species that are unlikely to be present in the sites we were working on like PARHES and EUDMAC. Two species that are unlike to be present in the study area.
# 
# first we tally the species 

species_table <- bat_combined %>%
  count(auto_id, sort = TRUE)

species_table

# ANPA there are 3k rows.



# remove noise and NoID
bat_clean <- bat_clean %>%
  filter(
    !str_to_upper(auto_id) %in% c("NOISE", "NOID"))

# remove large objcets to keep memory free
rm(kpro_2021_bat)
rm(kpro_2022_bat)
rm(kpro_2023_bat)

gc()



# rules for anpa 


# make all species tags into small caps 
bat_clean <- bat_combined %>%
  mutate(
    auto_id = str_to_lower(auto_id),
    alternate_1 = str_to_lower(alternate_1),
    alternate_2 = str_to_lower(alternate_2)
  )

# create sp column a exact copy of auto_id
bat_clean <- bat_clean %>%
  mutate(sp = auto_id)

# standadize empty spaces as NAs for alternate columns. 
bat_clean <- bat_clean %>%
  mutate(
    alternate_1 = na_if(alternate_1, ""),
    alternate_2 = na_if(alternate_2, "")
  )

# step one ANPA make "anpaepfu" when empfu is first alternative. 
# we need to update it to make all anpa when alternate_1 = Myevo into myoevo 

bat_clean <- bat_clean %>%
  mutate(
    auto_id = str_to_lower(auto_id),
    alternate_1 = str_to_lower(alternate_1),
    alternate_2 = str_to_lower(alternate_2),
    
    sp = case_when(
      
      # Stronger ANTPAL-like rows
      auto_id == "antpal" &
        fc >= 24 & fc <= 30 &
        fmin >= 22 & fmin <= 28 &
        dur >= 5 &
        fmax >= 43 & fmax <= 60 &
        (alternate_1 == "eptfus" | is.na(alternate_1)) ~ "antpal",
      
      # ANTPAL label, but alternate suggests CORTOW and call is short
      auto_id == "antpal" &
        alternate_1 == "cortow" &
        dur < 5 &
        fc >= 24 & fc <= 35 ~ "anpa_coto",
      
      # ANTPAL label, but alternate suggests MYOEVO and call has higher frequency
      auto_id == "antpal" &
        alternate_1 == "myoevo" &
        fc >= 32 ~ "myoevo",
      
      # ANTPAL label, but alternative suggests MYOTHY we don't have myothy in sites so lof
      auto_id == "antpal" &
        alternate_1 == "myothy" ~ "lof",
      
      # ANTPAL with EPTFUS alternative but not strong ANTPAL metrics
      auto_id == "antpal" &
        alternate_1 == "eptfus" ~ "anpa_epfu",
      
      # All other ANTPAL rows are low-frequency uncertain
      auto_id == "antpal" ~ "lof",
      
      # Everything else remains as originally labeled
      TRUE ~ auto_id
    ),
    
    inspect = case_when(
      sp %in% c("anpa_coto", "myoevo", "anpa_epfu", "lof_check") ~ "yes",
      TRUE ~ "no"
    ),
    
    rule_used = case_when(
      
      sp == "antpal" ~ 
        "ANTPAL retained: fc, fmin, fmax, and duration consistent with ANTPAL reference range",
      
      sp == "anpa_coto" ~ 
        "ANTPAL auto_id but CORTOW alternate and short duration; inspect as ANTPAL/CORTOW ambiguity",
      
      sp == "myoevo_check" ~ 
        "ANTPAL auto_id but MYOEVO alternate and higher fc; inspect as possible MYOEVO",
      
      sp == "myothy_check" ~ 
        "ANTPAL auto_id but MYOTHY alternate; inspect because MYTH fragments can mimic low-frequency calls",
      
      sp == "anpa_epfu" ~ 
        "ANTPAL auto_id with EPTFUS alternate but not strong ANTPAL metric match",
      
      sp == "lof_check" ~ 
        "ANTPAL auto_id but not enough support for confident ANTPAL",
      
      TRUE ~ "original auto_id kept"
    )
  )

# check if things worked.

bat_clean %>%
  filter(auto_id == "antpal") %>%
  count(auto_id, sp, alternate_1, inspect, rule_used, sort = TRUE)

table(bat_clean$sp)

# rows to inspect 
antpal_inspect <- bat_clean %>%
  filter(auto_id == "antpal", inspect == "yes") %>%
  select(
    in_file, auto_id, alternate_1, alternate_2,
    sp, inspect, rule_used,
    pulses, fc, dur, fmax, fmin, fmean, qual
  )

# CORTOW rules ------------------------------------------------------------


# CORTOW is often confused with ANPA, LANO, LACI, MYEV, also with Myca and MYth but whit these two cases it should remain as COTO. 

# we have about 1116 calls. 

coto<-bat_clean %>% 
  filter(auto_id == "cortow")

bat_clean %>% 
  filter(auto_id == "cortow") %>% 
  count(auto_id, sp, alternate_1)

# rule

bat_clean <- bat_clean %>%
  mutate(
    sp = case_when(
      
      # ------------------------------------------------------------
      # 1. Strong CORTOW-like calls: keep as CORTOW
      # ------------------------------------------------------------
      auto_id == "cortow" &
        fc >= 24 & fc <= 32 &
        dur >= 1.5 & dur <= 5.0 &
        fmax >= 34 & fmax <= 52 &
        fmin >= 21 & fmin <= 31 ~ "cortow",
      
      # ------------------------------------------------------------
      # 2. CORTOW with ANTPAL alternative
      # If duration is long, consider ANTPAL or CORTOW/ANTPAL ambiguity
      # ------------------------------------------------------------
      auto_id == "cortow" &
        alternate_1 == "antpal" &
        dur >= 5.5 &
        fc >= 24 & fc <= 30 ~ "antpal",
      
      auto_id == "cortow" &
        alternate_1 == "antpal" ~ "cortow_antpal",
      
      # ------------------------------------------------------------
      # 3. CORTOW with LASNOC alternative
      # LASNOC is usually longer and flatter
      # ------------------------------------------------------------
      auto_id == "cortow" &
        alternate_1 == "lasnoc" &
        dur >= 7 &
        fc >= 24 & fc <= 29 ~ "lasnoc",
      
      auto_id == "cortow" &
        alternate_1 == "lasnoc" ~ "cortow_lasnoc",
      
      # ------------------------------------------------------------
      # 4. CORTOW with LASCIN alternative
      # LASCIN usually lower fc and longer duration
      # ------------------------------------------------------------
      auto_id == "cortow" &
        alternate_1 == "lascin" &
        fc <= 24 &
        dur >= 7 ~ "lascin",
      
      auto_id == "cortow" &
        alternate_1 == "lascin" ~ "cortow_lascin",
      
      # ------------------------------------------------------------
      # 5. CORTOW with MYOEVO alternative
      # MYOEVO usually has higher fc and higher bandwidth
      # ------------------------------------------------------------
      auto_id == "cortow" &
        alternate_1 == "myoevo" &
        fc >= 33 ~ "myoevo",
      
      auto_id == "cortow" &
        alternate_1 == "myoevo" ~ "cortow_myoevo",
      
      # ------------------------------------------------------------
      # 6. CORTOW with MYOCAL or MYOTHY alternative
      # Your decision: keep as CORTOW, but flag for inspection
      # ------------------------------------------------------------
      auto_id == "cortow" &
        alternate_1 %in% c("myocal", "myothy") ~ "cortow",
      
      # ------------------------------------------------------------
      # 7. CORTOW rows that do not strongly match expected CORTOW metrics
      # ------------------------------------------------------------
      auto_id == "cortow" &
        (
          dur > 5 |
            fc < 24 | fc > 32 |
            fmax < 34 | fmax > 52
        ) ~ "cortow",
      
      # ------------------------------------------------------------
      # 8. Everything else unchanged
      # ------------------------------------------------------------
      TRUE ~ sp
    ),
    
    inspect = case_when(
      auto_id == "cortow" & sp != "cortow" ~ "yes",
      auto_id == "cortow" & alternate_1 %in% c("antpal", "lasnoc", "lascin", "myoevo", "myocal", "myothy") ~ "yes",
      auto_id == "cortow" & qual < 5 ~ "yes",
      TRUE ~ inspect
    ),
    
    rule_used = case_when(
      
      auto_id == "cortow" &
        fc >= 24 & fc <= 32 &
        dur >= 1.5 & dur <= 5.0 &
        fmax >= 34 & fmax <= 52 &
        fmin >= 21 & fmin <= 31 ~
        "CORTOW retained: metrics within expected CORTOW range",
      
      auto_id == "cortow" &
        alternate_1 == "antpal" &
        dur >= 5.5 &
        fc >= 24 & fc <= 30 ~
        "CORTOW recoded as ANTPAL: ANTPAL alternative with longer duration",
      
      auto_id == "cortow" &
        alternate_1 == "antpal" ~
        "CORTOW/ANTPAL ambiguity: inspect",
      
      auto_id == "cortow" &
        alternate_1 == "lasnoc" &
        dur >= 7 &
        fc >= 24 & fc <= 29 ~
        "CORTOW recoded as LASNOC: LASNOC alternative with longer duration",
      
      auto_id == "cortow" &
        alternate_1 == "lasnoc" ~
        "CORTOW/LASNOC ambiguity: inspect",
      
      auto_id == "cortow" &
        alternate_1 == "lascin" &
        fc <= 24 &
        dur >= 7 ~
        "CORTOW recoded as LASCIN: low fc and long duration",
      
      auto_id == "cortow" &
        alternate_1 == "lascin" ~
        "CORTOW/LASCIN ambiguity: inspect",
      
      auto_id == "cortow" &
        alternate_1 == "myoevo" &
        fc >= 33 ~
        "CORTOW recoded as MYOEVO: MYOEVO alternative with higher fc",
      
      auto_id == "cortow" &
        alternate_1 == "myoevo" ~
        "CORTOW/MYOEVO ambiguity: inspect",
      
      auto_id == "cortow" &
        alternate_1 %in% c("myocal", "myothy") ~
        "CORTOW retained but flagged: MYOCAL/MYOTHY alternative",
      
      auto_id == "cortow" &
        (
          dur > 5 |
            fc < 24 | fc > 32 |
            fmax < 34 | fmax > 52
        ) ~
        "CORTOW metrics outside conservative range: inspect",
      
      TRUE ~ rule_used
    )
  )

bat_clean %>%
  filter(auto_id == "cortow") %>%
  count(auto_id,alternate_1, sp, inspect, rule_used, sort = TRUE)

cortow_inspect <- bat_clean %>%
  filter(auto_id == "cortow", inspect == "yes") %>%
  select(
    in_file, auto_id, alternate_1, alternate_2,
    sp, inspect, rule_used,
    pulses, fc, dur, fmax, fmin, fmean, qual
  )

cortow_inspect # this are all the rows that need to be inspected. 


# rule for EUDMAC
# this batspecies is not likely at our sites. It is in idaho just not likely at our site. 
# all need to be recheked so we need to tag those rows as as yes for the inspect column and lof for the sp 
# but if the alternatives indicate a species that is likely such as lascin we should recode it as lascin and still check. 

bat_clean <- bat_clean %>%
  mutate(
    
    sp = case_when(
      auto_id == "eudmac" & !is.na(alternate_1) ~ alternate_1,
      auto_id == "eudmac" & is.na(alternate_1) & !is.na(alternate_2) ~ alternate_2,
      auto_id == "eudmac" ~ "lof",
      TRUE ~ sp
    ),
    
    inspect = case_when(
      auto_id == "eudmac" ~ "yes",
      TRUE ~ inspect
    ),
    
    rule_used = case_when(
      auto_id == "eudmac" & !is.na(alternate_1) ~
        "EUDMAC unlikely at study sites; recoded to alternate_1 and flagged for inspection",
      
      auto_id == "eudmac" & is.na(alternate_1) & !is.na(alternate_2) ~
        "EUDMAC unlikely at study sites; recoded to alternate_2 and flagged for inspection",
      
      auto_id == "eudmac" ~
        "EUDMAC unlikely at study sites; no alternative available, recoded as lof and flagged for inspection",
      
      TRUE ~ rule_used
    )
  )
# 

bat_clean %>%
  filter(auto_id == "eudmac") %>%
  count(auto_id,sp, inspect, rule_used)

eudmac_inspect <- bat_clean %>%
  filter(auto_id == "eudmac", inspect == "yes") %>%
  select(
    in_file, auto_id, alternate_1, alternate_2,
    sp, inspect, rule_used,
    pulses, fc, dur, fmax, fmin, fmean, qual
  )

eudmac_inspect


# MYCA rule ---------------------------------------------------------------

# we have about 3447 observations of mycal

# now we do Myotis Californicus MYCA or myocal
# this species is most often confused with mycil, myvol, mylu
myca<-bat_clean %>% 
  filter(auto_id == "myocal")

# effort ------------------------------------------------------------------


effort_days <- bat_combined %>%
  group_by(site, yr) %>%
  summarise(
    stard = min(noche),
    endd = max(noche),
    eff.days = as.numeric(difftime(max(noche), min(noche), units = "days"))
  )

effort_hrs <- bat_combined %>%
  group_by(site, noche, jday, yr) %>%
  summarise(stard = min(datetime), endd = max(datetime)) %>%
  mutate(eff.hrs = time_length(endd - stard, unit = "hours"))

# merge effort with bat combined 

bat_combined<- left_join(bat_combined, effort_hrs, by=c("site", "jday", "yr", "noche"))


keep<- c("AUTO.ID.", "PULSES", "site","noche","datetime", "yr","treatmt","trmt_bin","jday","eff.hrs") # cols to keep

bat_combined <- bat_combined %>% select(all_of(keep))

bat_combined <- bat_combined %>% rename(sp = AUTO.ID.)# change the auto.id to sp 



summary(bat_combined)




# count matrix ------------------------------------------------------------
# this is a matrix where we create a n column that tells us how many calls for each bat are there.
# daily counts

bm <- bat_combined %>% # 
  group_by(noche, sp, site,yr, treatmt, trmt_bin, eff.hrs) %>% 
  summarise(n = n(), .groups = 'drop') 

summary(bm)
head(bm)
# here we calculate the total calls per species to report in the results 
bm_summary <- bm %>%
  group_by(sp) %>%
  summarise(total_calls = sum(n)) %>%
  arrange(desc(total_calls))
# miller matrix -----------------------------------------------------------



# minutes activity 
# in here we calculate the minutes of activity insipired by Miller 2001 paper. 

bat_combined$rmins<-round(bat_combined$date_time, units="mins") #rounds to the nearest min

bm.miller<-bat_combined %>% #min of activity 
  group_by(site, AUTO.ID., noche, rmins) %>% 
  summarize(activity_min= n()) %>%  #calculate the num of min activity
  ungroup()


bm.miller.day <- bm.miller %>% # number of minutes active  by night. 
  group_by(site, noche, AUTO.ID.) %>%
  summarize(activity_min = sum(activity_min))

summary(bm.miller.day)


bm.miller<- bat_combined %>%
  # Extract the relevant columns and round to minute level
  mutate(rmins = round(date_time, units = "mins")) %>%
  # Remove duplicate entries for the same site, date, and minute
  distinct(site, sp, noche, rmins, .keep_all = TRUE) %>%
  # Group by site, date, and minute
  group_by(site,sp, noche, rmins) %>%
  # Summarize to count the number of unique minutes
  summarize(activity_min = n_distinct(rmins), .groups = 'drop')

bm.miller.day <- bm.miller %>% # number of minutes active  by night. 
  group_by(site, noche, sp) %>%
  summarize(activity_min = sum(activity_min))

head(bm.miller.day)
summary(bm.miller.day)



# correcting for abundance ------------------------------------------------

# This code maps sites into defined pairs, groups data by night, site pair, and species, calculates normalized bat activity by dividing experimental activity (lit treatment) by the mean control activity (dark treatment), and filters out noise and unidentified species. It also adds a Julian day column to the resulting dataset.

# Define the pairs
site_pairs <- list(
  c("long01", "long02"),
  c("long03", "long04"),
  c("iron01", "iron02"),
  c("iron03", "iron04"),
  c("iron05", "iron06")
)

# Add a column to map pairs
bm <- bm %>%
  mutate(
    pair_group = case_when(
      site %in% c("long01", "long02") ~ "long01:long02",
      site %in% c("long03", "long04") ~ "long03:long04",
      # site %in% c("long05") ~ "long05",
      site %in% c("iron01", "iron02") ~ "iron01:iron02",
      site %in% c("iron03", "iron04") ~ "iron03:iron04",
      site %in% c("iron05", "iron06") ~ "iron05:iron06",
      TRUE ~ NA_character_
    )
  )

normalized_bm <- bm %>%
  group_by(noche, pair_group, sp) %>%  # Add species (sp) to the grouping
  summarize(
    control_mean = ifelse(all(is.na(n[treatmt == "dark"])), 
                          .01,  # Replace NaN with .01 when no data exists for control
                          mean(n[treatmt == "dark"], na.rm = TRUE)),  # Average control activity per species,  # Average control activity per species
    control_activity = sum(n[treatmt == "dark"], na.rm = TRUE),  # Total control activity per species
    experimental_activity = sum(n[treatmt == "lit"], na.rm = TRUE),  # Total experimental activity per species
    .groups = "drop"
  ) %>%
  mutate(
    normalized_activity = experimental_activity / (control_mean+ experimental_activity)  # Calculate normalized activity
  ) %>% 
  mutate(
    bin_act = ifelse(experimental_activity>=control_mean, 1, 0)
  ) %>% 
  mutate(
  j_diff= experimental_activity - control_activity) # difference between experimental and control activity.
  
  


# View the result
normalized_bm

summary(normalized_bm) # check for NAs
head(normalized_bm,100)

normalized_bm$jday<-yday(normalized_bm$noche)
normalized_bm$year<-year(normalized_bm$noche)

# View the result
normalized_bm
# there are some rows that have NAs. that come from control sites not paired with any experimental site.we remove them. 
rows_with_na <- normalized_bm %>% filter(if_any(everything(), is.na)) 
normalized_bm<- normalized_bm %>% drop_na()
normalized_bm <- normalized_bm[!(normalized_bm$sp %in% c("Noise","NoID")), ]# filter out Noise and NoID rows. 
#remove all rows for long05
normalized_bm <- normalized_bm %>% filter(!pair_group %in% "long05")

summary(normalized_bm)
head(normalized_bm,100)

hist(normalized_bm$j_diff)
range(normalized_bm$j_diff) 

# show rows with 20 largest j_diff 
rowlarge<- normalized_bm %>%
  arrange(desc(j_diff)) %>%
  slice_head(n = 20)

print(rowlarge)
# save rowlarge as .csv
rowlarge <- rowlarge %>% select(-c(control_mean, control_activity, experimental_activity)) # remove unwanted columns

write.csv(rowlarge, file = 'data_for_analysis/prep_for_glmm/rowlarge.csv', row.names = F) # raw combine data 



# outputs -----------------------------------------------------------------


# dir.create("data_for_analysis/prep_for_glm", showWarnings = FALSE) # just run if the dir is abscent

write.csv(bat_combined, file = 'data_for_analysis/prep_for_glmm/bat_combined.csv', row.names = F) # raw combine data 
write.csv(bm, file = 'data_for_analysis/prep_for_glmm/bm.csv', row.names = F) #daily counts
write.csv(bm.miller.day, file = "data_for_analysis/prep_for_glmm/bm.miller.day.csv") # miller Ai index data
write.csv(normalized_bm, file = "data_for_analysis/prep_for_glmm/normalized_bm.csv")


# Create a README file with information about the script
readme_content <- "Carlos Linares 8/01/2024, 12/12/2024 
This directory contains the bat_combined.csv file which was created using the script prep_for_glmm.R combines bat species call abundance data. This script merges 2021-23 data that was previously scanned with Kaleidoscope pro

bat_combined.csv - process data no counts (Update: 12/2/2024 some modifications to date time column.) 
bm.csv - counts of bat calls by day from 2021 to 2023 all sites
bm.miller.day - number of minutes of activity by day  for 2021-2023 data all sites (last update 9/23/2024)
normalized_bm.csv - normalized bat activity by dividing experimental activity (lit treatment) by the mean control activity (dark treatment) by species and noche(date) "
# Write the README content to a file
writeLines(readme_content, "data_for_analysis/prep_for_glmm/README.txt")







# plots -------------------------------------------------------------------

#-------------------------- calls by site and year


bat_summary <- bat_combined %>%
  filter(!AUTO.ID. %in% c("Noise", "NoID")) %>%  # Filter out Noise and NoID tags
  group_by(site, yr, AUTO.ID.) %>%
  summarise(count = n()) %>%
  ungroup()


ggplot(bat_summary, aes(x = yr, y = count, fill = AUTO.ID.)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ site, scales = "free_y") +
  labs(title = "species by sites and years",
       x = "Year",
       y = "Count",
       fill = "Tag (AUTO.ID.)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#-------------------------- calls by jday and treatment

filtered_bat_combined <- bat_combined %>%
  filter(!AUTO.ID. %in% c("Noise", "EUDMAC"))

# Summarize the number of pulses per Julian day and treatment
summary_data <- filtered_bat_combined %>%
  group_by(jday, treatmt, yr) %>%
  summarise(count = n()) %>%
  ungroup()

# Create the plot
ggplot(summary_data, aes(x = jday, y = count, col = treatmt)) +
  geom_point() +
  facet_wrap(~ yr+ treatmt , scales = "free_y") +  labs(title = "Call activity by Julian Day and Treatment",
       x = "Julian Day",
       y = "Number of calls") +
  geom_vline(xintercept = 180, linetype = "dashed", color = "red") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#-------------------------- species calls by Julian day and treatment

summary_data <- filtered_bat_combined %>%
  group_by(jday, treatmt, AUTO.ID.) %>%
  summarise(count = n()) %>%
  ungroup()


# Create the plot
ggplot(summary_data, aes(x = jday, y = count, col = treatmt)) +
  geom_point() +
  facet_wrap( ~  AUTO.ID. + treatmt, scales = "free_y") +  labs(title = "Call activity by Julian Day and Treatment", x = "Julian Day", y = "call counts") +
  geom_vline(xintercept = 180,
             linetype = "dashed",
             color = "red") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))





summary_data <- filtered_bat_combined %>%
  group_by(jday, treatmt, AUTO.ID., yr) %>%
  summarise(count = n()) %>%
  ungroup()

# Create the plot
ggplot(summary_data, aes(x = jday, y = count, col = treatmt)) +
  geom_point() +
  facet_wrap(~ AUTO.ID.+ yr , scales = "free_y") +  labs(title = "Call activity by Sp and year",
                                                               x = "Julian Day",
                                                               y = "Number of Pulses") +
  geom_vline(xintercept = 180, linetype = "dashed", color = "red") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


# pulses by treatment 

ggplot(bat_combined, aes(x = treatmt, y = scale(PULSES), fill = treatmt)) +
  geom_boxplot() +
  labs(title = "Distribution of Bat Pulses by Treatment",
       x = "Treatment",
       y = "Number of Pulses") +
  theme_minimal()


# plots_normalized_activity  ----------------------------------------------

ggplot(normalized_bm, aes(x = normalized_activity)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black") +
  facet_wrap(~ year) +  # Facet by year
  labs(
    title = "Histogram of Normalized Activity by Year",
    x = "Normalized Activity",
    y = "Frequency"
  ) +
  theme_minimal()

ggplot(normalized_bm, aes(x = normalized_activity, colour = factor(year))) + 
  geom_density() +
  labs(
    title = "Density Plot of Normalized Activity by Year",
    x = "Normalized Activity",
    y = "Density"
  ) +
  theme_minimal()

ggplot(normalized_bm, aes(x = factor(year), y = normalized_activity)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(
    title = "Boxplot of Normalized Activity by Year",
    x = "Year",
    y = "Normalized Activity"
  ) +
  theme_minimal()

# Print the plot


# save --------------------------------------------------------------------

save.image("working_env/prep_for_glm.RData")

# Session info ------------------------------------------------------------



# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 22631)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: America/Denver
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] beepr_2.0       lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5    
# [8] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.5       compiler_4.4.1     tidyselect_1.2.1   scales_1.3.0       R6_2.5.1           generics_0.1.3    
# [7] sjlabelled_1.2.0   knitr_1.46         sjPlot_2.8.16      insight_0.20.2     munsell_0.5.1      pillar_1.9.0      
# [13] tzdb_0.4.0         sjstats_0.19.0     rlang_1.1.3        utf8_1.2.4         stringi_1.8.4      performance_0.11.0
# [19] xfun_0.44          audio_0.1-11       ggeffects_1.6.0    datawizard_0.12.1  timechange_0.3.0   cli_3.6.2         
# [25] withr_3.0.0        magrittr_2.0.3     grid_4.4.1         rstudioapi_0.16.0  hms_1.1.3          lifecycle_1.0.4   
# [31] sjmisc_2.8.10      vctrs_0.6.5        glue_1.7.0         fansi_1.0.6        colorspace_2.1-0   tools_4.4.1       
# [37] pkgconfig_2.0.3   
# 
# 
# 




# --------------------------- trash ----------------



midnight_strings <- datetime[grepl("^\\d{4}-\\d{2}-\\d{2} 00:00:00$", datetime)] # make sure the midnight strings are there. 

t<-ymd_hms(head(midnight_strings), tz = "America/Denver")

formatted_times <- format(t, "%Y-%m-%d %H:%M:%S %Z")
print(formatted_times)\\\



bat_combined$date_time<-ymd_hms(datetime) # add to data. 
sum(is.na(bat_combined$date_time)) # check for NAs. 


# maybe too complex normalization

# Normalize activity
normalized_bm <- bm %>%
  group_by(noche, pair_group, sp) %>%  # Add species (sp) to the grouping
  summarize(
    control_mean = ifelse(all(is.na(n[treatmt == "dark"])), 
                          .01,  # Replace NaN with .01 when no data exists for control
                          mean(n[treatmt == "dark"], na.rm = TRUE)),  # Average control activity per species
    experimental_activity = sum(n[treatmt == "lit"], na.rm = TRUE),  # Total experimental activity per species
    .groups = "drop"
  ) %>%
  mutate(
    normalized_activity = ifelse(control_mean == 0, 
                                 .01,  # Assign a small value when control_mean is zero
                                 experimental_activity / control_mean)  # Otherwise, calculate normally
  )%>%
  group_by(sp) %>%
  mutate(
    # Min-Max normalization to scale between 0 and 1
    normalized_activity = (normalized_activity - min(normalized_activity)) / 
      (max(normalized_activity) - min(normalized_activity))
  ) %>%
  ungroup()
