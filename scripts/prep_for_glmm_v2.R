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
library(corrplot)


# functions ---------------------------------------------------------------

# funtion to scale. 
scale_by_2sd_tidy <- function(data, variables_to_scale) {
  # Keep only variables that exist and are numeric
  valid_vars <- variables_to_scale[variables_to_scale %in% names(data) & sapply(data[variables_to_scale], is.numeric)]
  
  if (length(valid_vars) == 0) {
    warning("No valid numeric variables found to scale.")
    return(data)
  }
  
  data <- data %>%
    mutate(across(all_of(valid_vars),
                  ~ (. - mean(., na.rm = TRUE)) / (2 * sd(., na.rm = TRUE)),
                  .names = "{.col}_s"))
  
  return(data)
}


# kpro_data --------------------------------------------------------------

# load
# 
# load("working_env/prep_for_glm.RData")


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




# site --------------------------------------------------------------------


bat_combined$site<-str_extract(bat_combined$outdir, "[A-Za-z]{3,4}\\d{2}")

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
# update on this. I realized that sometimes when I do this the correct date is not input. However, the date.12 column in Kaleidocope does the right thing. Sow we updated this section to use the date.12 instead of the following lines. 

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
midnight_times


#year

bat_combined$yr<-year(bat_combined$date)
unique(bat_combined$yr) #we check the three years are present.


# treatment column 

litsites<-c("iron01","iron03","iron05","long01","long03")


bat_combined$treatmt <- ifelse(bat_combined$site %in% litsites , "lit", "dark") # this makes a treatment variable.

bat_combined$trmt_bin<- ifelse(bat_combined$treatmt== "lit", 1, 0)


bat_combined$jday<-lubridate::yday(bat_combined$noche) # julian day


summary(bat_combined)

# effort -------------------------------------------------------------------

# now when load the effort from the actual files. 

effort_days<- read.csv("data_for_analysis/effort/effort_days.csv")
effort_hrs<- read.csv("data_for_analysis/effort/effort_hrs.csv")

# below is the old way of calculating effort. 
# # we have to calculate effort before fitering becaue we don't have the noise and noID calls. 
# 
# effort_days <- bat_combined %>%
#   group_by(site, yr) %>%
#   summarise(
#     stard = min(noche),
#     endd = max(noche),
#     eff.days = as.numeric(difftime(max(noche), min(noche), units = "days"))
#   )
# 
# effort_hrs <- bat_combined %>%
#   group_by(site, noche, jday, yr) %>%
#   summarise(stard = min(datetime), endd = max(datetime)) %>%
#   mutate(eff.hrs = time_length(endd - stard, unit = "hours"))

# species rules  ----------------------------------------------------------

# now we are going to create a series of rules to make sure we don't have species that are unlikely to be present in the sites we were working on like PARHES and EUDMAC. Two species that are unlike to be present in the study area.
# 
# first we tally the species 

species_table <- bat_combined %>%
  count(auto_id, sort = TRUE)

species_table

# ANPA there are 3k rows.





# remove large objcets to keep memory free
rm(kpro_2021_bat)
rm(kpro_2022_bat)
rm(kpro_2023_bat)

gc()



# rules for anpa 

# first remove noise and NoID
bat_clean <- bat_combined %>%
  filter(
    !str_to_upper(auto_id) %in% c("NOISE", "NOID")) 

# make all species tags into small caps 
bat_clean <- bat_clean %>%
  mutate(
    auto_id = str_to_lower(auto_id),
    alternate_1 = str_to_lower(alternate_1),
    alternate_2 = str_to_lower(alternate_2)
  )

# create sp column a exact copy of auto_id
bat_clean <- bat_clean %>%
  mutate(sp = auto_id)

# standardize empty spaces as NAs for alternate columns. 
bat_clean <- bat_clean %>%
  mutate(
    alternate_1 = na_if(alternate_1, ""),
    alternate_2 = na_if(alternate_2, "")
  )

# step one ANPA make "anpaepfu" when empfu is first alternative. 
# we need to update it to make all anpa when alternate_1 = Myevo into myoevo 

# first how many ANPA we have 
unique(bat_clean$auto_id)

antpal<-bat_clean %>% 
  filter(auto_id == "antpal") # we have 3821 antpal 

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
        fc >= 24 & fc <= 35 ~ "cortow",
      
      # ANTPAL label, but alternate suggests MYOEVO and call has higher frequency
      auto_id == "antpal" &
        alternate_1 == "myoevo" &
        fc >= 32 ~ "myoevo",
      
      # ANTPAL label, but alternative suggests MYOTHY we don't have myothy in sites so lof
      auto_id == "antpal" &
        alternate_1 == "myothy" ~ "antpal",
      
      # ANTPAL with EPTFUS alternative but not strong ANTPAL metrics
      auto_id == "antpal" &
        alternate_1 == "eptfus" ~ "eptfus",
      
      # All other ANTPAL rows are low-frequency uncertain, but for now we keep them as antpal 
      auto_id == "antpal" ~ "antpal",
      
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

table(bat_clean$sp) # after filtering it went from 3821 to 282 antpal. 

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

# we have about 1116 calls for coto

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
        alternate_1 == "antpal" ~ "cortow",
      
      # ------------------------------------------------------------
      # 3. CORTOW with LASNOC alternative
      # LASNOC is usually longer and flatter
      # ------------------------------------------------------------
      auto_id == "cortow" &
        alternate_1 == "lasnoc" &
        dur >= 7 &
        fc >= 24 & fc <= 29 ~ "lasnoc",
      
      auto_id == "cortow" &
        alternate_1 == "lasnoc" ~ "cortow",
      
      # ------------------------------------------------------------
      # 4. CORTOW with LASCIN alternative
      # LASCIN usually lower fc and longer duration
      # ------------------------------------------------------------
      auto_id == "cortow" &
        alternate_1 == "lascin" &
        fc <= 24 &
        dur >= 7 ~ "lascin",
      
      auto_id == "cortow" &
        alternate_1 == "lascin" ~ "cortow",
      
      # ------------------------------------------------------------
      # 5. CORTOW with MYOEVO alternative
      # MYOEVO usually has higher fc and higher bandwidth
      # ------------------------------------------------------------
      auto_id == "cortow" &
        alternate_1 == "myoevo" &
        fc >= 33 ~ "myoevo",
      
      auto_id == "cortow" &
        alternate_1 == "myoevo" ~ "cortow",
      
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

table(bat_clean$sp) # after filter we went from 1116 to 2271 is more I don't know why. 


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

table(bat_clean$sp) # eudmac was recoded or marke as lof

# MYCA rule ---------------------------------------------------------------

# we have about 3447 observations of mycal

# now we do Myotis Californicus MYCA or myocal
# this species is most often confused with mycil, myvol, mylu
myca<-bat_clean %>% 
  filter(auto_id == "myocal")


likely_myotis_alt <- c("myocil", "myovol", "myoluc")

bat_clean <- bat_clean %>%
  mutate(
    
    sp = case_when(
      
      # ------------------------------------------------------------
      # 1. MYOCAL auto_id but alternate_1 supports MYOLUC
      # MYOLUC tends to have lower fc and longer duration.
      # ------------------------------------------------------------
      auto_id == "myocal" &
        alternate_1 == "myoluc" &
        fc >= 39 & fc <= 44 &
        fmax >= 55 & fmax <= 85 &
        dur >= 4.5 ~ "myoluc",
      
      # ------------------------------------------------------------
      # 2. MYOCAL auto_id but alternate_2 supports MYOLUC
      # ------------------------------------------------------------
      auto_id == "myocal" &
        alternate_2 == "myoluc" &
        fc >= 39 & fc <= 44 &
        fmax >= 55 & fmax <= 85 &
        dur >= 4.5 ~ "myoluc",
      
      # ------------------------------------------------------------
      # 3. MYOCAL auto_id but alternate_1 supports MYOVOL
      # MYOVOL around 40-43 kHz, often shorter/steeper than MYOLUC.
      # ------------------------------------------------------------
      auto_id == "myocal" &
        alternate_1 == "myovol" &
        fc >= 40 & fc <= 46 &
        fmax >= 55 & fmax <= 95 &
        dur >= 3.0 & dur <= 6.0 ~ "myovol",
      
      # ------------------------------------------------------------
      # 4. MYOCAL auto_id but alternate_2 supports MYOVOL
      # This catches the row you mentioned: myoyum + myovol.
      # ------------------------------------------------------------
      auto_id == "myocal" &
        alternate_2 == "myovol" &
        fc >= 40 & fc <= 46 &
        fmax >= 55 & fmax <= 95 &
        dur >= 3.0 & dur <= 6.0 ~ "myovol",
      
      # ------------------------------------------------------------
      # 5. MYOCAL auto_id but alternate_1 supports MYOCIL
      # MYOCIL tends to have high fmax.
      # ------------------------------------------------------------
      auto_id == "myocal" &
        alternate_1 == "myocil" &
        fc >= 40 & fc <= 45 &
        fmax >= 75 &
        dur >= 2.5 & dur <= 5.5 ~ "myocil",
      
      # ------------------------------------------------------------
      # 6. MYOCAL auto_id but alternate_2 supports MYOCIL
      # ------------------------------------------------------------
      auto_id == "myocal" &
        alternate_2 == "myocil" &
        fc >= 40 & fc <= 45 &
        fmax >= 75 &
        dur >= 2.5 & dur <= 5.5 ~ "myocil",
      
      # ------------------------------------------------------------
      # 7. Strong MYOCAL-like calls, but only if no likely target
      # alternative is available.
      # ------------------------------------------------------------
      auto_id == "myocal" &
        !(alternate_1 %in% likely_myotis_alt) &
        !(alternate_2 %in% likely_myotis_alt) &
        fc >= 47 & fc <= 53 &
        fmax >= 70 &
        fmin >= 40 &
        dur >= 2.5 & dur <= 5.5 ~ "myocal",
      
      # ------------------------------------------------------------
      # 8. MYOYUM alternatives without a target alternate
      # ------------------------------------------------------------
      auto_id == "myocal" &
        (alternate_1 == "myoyum" | alternate_2 == "myoyum") &
        !(alternate_1 %in% likely_myotis_alt) &
        !(alternate_2 %in% likely_myotis_alt) ~ "myocal",
      
      # ------------------------------------------------------------
      # 9. Short Myotis calls are weak for species-level ID
      # ------------------------------------------------------------
      auto_id == "myocal" &
        dur < 3.5 ~ "myocal",
      
      # ------------------------------------------------------------
      # 10. Lower-frequency 40 kHz Myotis metrics without a supported
      # target alternative
      # ------------------------------------------------------------
      auto_id == "myocal" &
        fc >= 39 & fc <= 45 ~ "myocal",
      
      # ------------------------------------------------------------
      # 11. Any remaining MYOCAL rows
      # ------------------------------------------------------------
      auto_id == "myocal" ~ "myocal",
      
      TRUE ~ sp
    ),
    
    inspect = case_when(
      auto_id == "myocal" ~ "yes",
      TRUE ~ inspect
    ),
    
    rule_used = case_when(
      
      auto_id == "myocal" &
        alternate_1 == "myoluc" &
        fc >= 39 & fc <= 44 &
        fmax >= 55 & fmax <= 85 &
        dur >= 4.5 ~
        "MYOCAL recoded to MYOLUC: alternate_1 supports MYOLUC and metrics fit lower-frequency longer Myotis",
      
      auto_id == "myocal" &
        alternate_2 == "myoluc" &
        fc >= 39 & fc <= 44 &
        fmax >= 55 & fmax <= 85 &
        dur >= 4.5 ~
        "MYOCAL recoded to MYOLUC: alternate_2 supports MYOLUC and metrics fit lower-frequency longer Myotis",
      
      auto_id == "myocal" &
        alternate_1 == "myovol" &
        fc >= 40 & fc <= 46 &
        fmax >= 55 & fmax <= 95 &
        dur >= 3.0 & dur <= 6.0 ~
        "MYOCAL recoded to MYOVOL: alternate_1 supports MYOVOL and metrics fit 40 kHz Myotis",
      
      auto_id == "myocal" &
        alternate_2 == "myovol" &
        fc >= 40 & fc <= 46 &
        fmax >= 55 & fmax <= 95 &
        dur >= 3.0 & dur <= 6.0 ~
        "MYOCAL recoded to MYOVOL: alternate_2 supports MYOVOL and metrics fit 40 kHz Myotis",
      
      auto_id == "myocal" &
        alternate_1 == "myocil" &
        fc >= 40 & fc <= 45 &
        fmax >= 75 &
        dur >= 2.5 & dur <= 5.5 ~
        "MYOCAL recoded to MYOCIL: alternate_1 supports MYOCIL and high fmax supports MYOCIL-like call",
      
      auto_id == "myocal" &
        alternate_2 == "myocil" &
        fc >= 40 & fc <= 45 &
        fmax >= 75 &
        dur >= 2.5 & dur <= 5.5 ~
        "MYOCAL recoded to MYOCIL: alternate_2 supports MYOCIL and high fmax supports MYOCIL-like call",
      
      auto_id == "myocal" &
        !(alternate_1 %in% likely_myotis_alt) &
        !(alternate_2 %in% likely_myotis_alt) &
        fc >= 47 & fc <= 53 &
        fmax >= 70 &
        fmin >= 40 &
        dur >= 2.5 & dur <= 5.5 ~
        "MYOCAL retained: metrics broadly consistent with higher-frequency MYOCAL and no likely target alternative",
      
      auto_id == "myocal" &
        (alternate_1 == "myoyum" | alternate_2 == "myoyum") &
        !(alternate_1 %in% likely_myotis_alt) &
        !(alternate_2 %in% likely_myotis_alt) ~
        "MYOCAL with MYOYUM alternative only: recoded as MYSP because target likely alternatives are not supported",
      
      auto_id == "myocal" &
        dur < 3.5 ~
        "MYOCAL recoded as MYSP: short Myotis calls are weak for species-level ID",
      
      auto_id == "myocal" &
        fc >= 39 & fc <= 45 ~
        "MYOCAL recoded as MYSP: lower-frequency 40 kHz Myotis metrics but no supported target alternative",
      
      auto_id == "myocal" ~
        "MYOCAL recoded as MYSP: MYOCAL less likely at study sites and no strong supported alternative",
      
      TRUE ~ rule_used
    )
  )

bat_clean %>%
  filter(auto_id == "myocal") %>%
  count(auto_id, alternate_1, alternate_2, sp, inspect, rule_used, sort = TRUE)

bat_clean %>%
  filter(auto_id == "myocal") %>%
  group_by(sp) %>%
  summarise(
    n = n(),
    mean_fc = mean(fc, na.rm = TRUE),
    mean_fmax = mean(fmax, na.rm = TRUE),
    mean_fmin = mean(fmin, na.rm = TRUE),
    mean_dur = mean(dur, na.rm = TRUE),
    min_dur = min(dur, na.rm = TRUE),
    max_dur = max(dur, na.rm = TRUE),
    .groups = "drop"
  )

myocal_inspect <- bat_clean %>%
  filter(auto_id == "myocal", inspect == "yes") %>%
  select(
    in_file, auto_id, alternate_1, alternate_2,
    sp, inspect, rule_used,
    pulses, fc, dur, fmax, fmin, fmean, qual
  )

table(bat_clean$sp) #it just went from 3447 to 3029 after the filtering. needs extra work might need to be fiter out.  


# 

# rules for parhes --------------------------------------------------------
# this species hasn't been reported in the area according to the idaho department of fish and game. We should reconsider many of othe calls. It often gets confused with myotis but at the frequency range where it calls it is hard to distinguish. 


parhes<-bat_clean %>% 
  filter(auto_id == "parhes") # we have about 5541 observations. 

bat_clean <- bat_clean %>%
  mutate(
    
    sp = case_when(
      
      # ------------------------------------------------------------
      # 1. Stronger PARHES-like calls
      # PARHES tends to have higher fc and low frequency near/above 39–45 kHz.
      # ------------------------------------------------------------
      auto_id == "parhes" &
        fc >= 44 & fc <= 50 &
        fmin >= 39 &
        fmax >= 48 & fmax <= 75 &
        dur >= 3.5 & dur <= 7.5 ~ "parhes",
      
      # ------------------------------------------------------------
      # 2. PARHES auto_id but alternate supports MYOLUC
      # MYOLUC is lower-frequency 40 kHz Myotis and tends to have longer duration.
      # ------------------------------------------------------------
      auto_id == "parhes" &
        alternate_1 == "myoluc" &
        fc >= 38 & fc <= 44 &
        dur >= 5.0 ~ "myoluc",
      
      auto_id == "parhes" &
        alternate_2 == "myoluc" &
        fc >= 38 & fc <= 44 &
        dur >= 5.0 ~ "myoluc",
      
      # ------------------------------------------------------------
      # 3. PARHES auto_id but alternate supports MYOVOL
      # MYOVOL is also around 40–43 kHz but often shorter/steeper than MYOLUC.
      # ------------------------------------------------------------
      auto_id == "parhes" &
        alternate_1 == "myovol" &
        fc >= 40 & fc <= 45 &
        dur >= 3.0 & dur <= 5.5 ~ "myovol",
      
      auto_id == "parhes" &
        alternate_2 == "myovol" &
        fc >= 40 & fc <= 45 &
        dur >= 3.0 & dur <= 5.5 ~ "myovol",
      
      # ------------------------------------------------------------
      # 4. PARHES auto_id but alternate supports MYOCAL
      # MYOCAL is a higher-frequency Myotis; use only with higher fc/fmin.
      # ------------------------------------------------------------
      auto_id == "parhes" &
        alternate_1 == "myocal" &
        fc >= 47 & fc <= 53 &
        fmin >= 40 &
        fmax >= 60 ~ "myocal",
      
      auto_id == "parhes" &
        alternate_2 == "myocal" &
        fc >= 47 & fc <= 53 &
        fmin >= 40 &
        fmax >= 60 ~ "myocal",
      
      # ------------------------------------------------------------
      # 5. Low-frequency Myotis-like calls without strong species support
      # ------------------------------------------------------------
      auto_id == "parhes" &
        fc >= 38 & fc <= 44 &
        dur >= 4.0 ~ "mysp",
      
      # ------------------------------------------------------------
      # 6. Very short or weak PARHES rows: insufficient for species-level ID
      # ------------------------------------------------------------
      auto_id == "parhes" &
        dur < 3.5 ~ "hif",
      
      # ------------------------------------------------------------
      # 7. Any remaining PARHES rows
      # ------------------------------------------------------------
      auto_id == "parhes" ~ "hif",
      
      TRUE ~ sp
    ),
    
    inspect = case_when(
      auto_id == "parhes" ~ "yes",
      TRUE ~ inspect
    ),
    
    rule_used = case_when(
      
      auto_id == "parhes" &
        fc >= 44 & fc <= 50 &
        fmin >= 39 &
        fmax >= 48 & fmax <= 75 &
        dur >= 3.5 & dur <= 7.5 ~
        "PARHES retained: metrics broadly consistent with PARHES high-frequency call; flagged for inspection",
      
      auto_id == "parhes" &
        alternate_1 == "myoluc" &
        fc >= 38 & fc <= 44 &
        dur >= 5.0 ~
        "PARHES recoded to MYOLUC: alternate_1 supports MYOLUC and metrics fit lower-frequency longer Myotis",
      
      auto_id == "parhes" &
        alternate_2 == "myoluc" &
        fc >= 38 & fc <= 44 &
        dur >= 5.0 ~
        "PARHES recoded to MYOLUC: alternate_2 supports MYOLUC and metrics fit lower-frequency longer Myotis",
      
      auto_id == "parhes" &
        alternate_1 == "myovol" &
        fc >= 40 & fc <= 45 &
        dur >= 3.0 & dur <= 5.5 ~
        "PARHES recoded to MYOVOL: alternate_1 supports MYOVOL and metrics fit 40 kHz Myotis",
      
      auto_id == "parhes" &
        alternate_2 == "myovol" &
        fc >= 40 & fc <= 45 &
        dur >= 3.0 & dur <= 5.5 ~
        "PARHES recoded to MYOVOL: alternate_2 supports MYOVOL and metrics fit 40 kHz Myotis",
      
      auto_id == "parhes" &
        alternate_1 == "myocal" &
        fc >= 47 & fc <= 53 &
        fmin >= 40 &
        fmax >= 60 ~
        "PARHES recoded to MYOCAL: alternate_1 supports MYOCAL and metrics fit higher-frequency Myotis",
      
      auto_id == "parhes" &
        alternate_2 == "myocal" &
        fc >= 47 & fc <= 53 &
        fmin >= 40 &
        fmax >= 60 ~
        "PARHES recoded to MYOCAL: alternate_2 supports MYOCAL and metrics fit higher-frequency Myotis",
      
      auto_id == "parhes" &
        fc >= 38 & fc <= 44 &
        dur >= 4.0 ~
        "PARHES recoded as MYSP: lower-frequency 40 kHz Myotis-like metrics but no supported species-level alternative",
      
      auto_id == "parhes" &
        dur < 3.5 ~
        "PARHES recoded as HIF: short call, insufficient for species-level ID",
      
      auto_id == "parhes" ~
        "PARHES recoded as HIF: PARHES auto_id not strongly supported by metrics or alternatives",
      
      TRUE ~ rule_used
    )
  )

bat_clean %>%
  filter(auto_id == "parhes") %>%
  count( auto_id, alternate_1, alternate_2, sp, inspect, rule_used, sort = TRUE)

table(bat_clean$sp) # it went from 5000 to 638 taht we can easily filter out. 

# effort ------------------------------------------------------------------
# # we have to calculate effort before filtering because we don't have the noise and noID calls. te following section was sent up the script
# 
# effort_days <- bat_combined %>%
#   group_by(site, yr) %>%
#   summarise(
#     stard = min(noche),
#     endd = max(noche),
#     eff.days = as.numeric(difftime(max(noche), min(noche), units = "days"))
#   )
# 
# effort_hrs <- bat_combined %>%
#   group_by(site, noche, jday, yr) %>%
#   summarise(stard = min(datetime), endd = max(datetime)) %>%
#   mutate(eff.hrs = time_length(endd - stard, unit = "hours"))

# merge effort with bat combined 

# but first these need to be tibble or the same class of object. 

bat_clean <- bat_clean %>%
  as_tibble() %>%
  mutate(
    site = as.character(site),
    jday = as.integer(jday),
    yr = as.integer(yr),
    noche = as.Date(noche)
  )

effort_hrs <- effort_hrs %>%
  ungroup() %>%
  as_tibble() %>%
  mutate(
    site = as.character(site),
    jday = as.integer(jday),
    yr = as.integer(year),
    noche = as.Date(monitoring_night)
  )

bat_clean<- left_join(bat_clean, effort_hrs, by=c("site", "jday", "yr", "noche"))
bat_clean<- left_join(bat_clean, effort_days, by=c("site", "year"))
summary(bat_clean)

effort_na <- bat_clean %>%
  filter(is.na(eff.hrs) | is.na(eff.days)) # we see the missing information rows come from acoustic files that will be removed later.


bat_clean %>%
  filter(is.na(eff.hrs) | is.na(eff.days)) %>%
  count(site, noche, year, monitoring_night, sort = TRUE)

# bat_clean_v2 ------------------------------------------------------------

# here we filter out unnecessary columns and keep only the ones we need for the analysis.

keep<- c("in_file", "date", "time","noche","datetime","auto_id","manual_id","site","noche", "datetime", "yr", "treatmt","trmt_bin", "jday","sp","eff.hrs","eff.days" ) # cols to keep

bat_clean_v2 <- bat_clean %>% select(all_of(keep))

# it seems we have acoustic files in the bat_clean_v2 that might need to be removed

bat_clean_v2 %>%
  mutate(
    file_type = case_when(
      str_detect(in_file, "__0__") ~ "acoustic",
      str_detect(in_file, "__1__") ~ "ultrasonic",
      TRUE ~ "other_or_unknown"
    )
  ) %>%
  count(file_type, sort = TRUE) 

bat_clean_v2 %>%
  filter(str_detect(in_file, "__0__")) %>% # see acousticfiles. 
  count(in_file, sort = TRUE) %>%
  head(20)

# we filter and keep only the ultrasnic files. 
bat_clean_v2 <- bat_clean_v2 %>%
  filter(str_detect(in_file, "__1__"))

summary(bat_clean_v2)
glimpse(bat_clean_v2)


# now we want to see how many files are duplicated. 

x<-bat_clean_v2 %>%
  count(in_file, sort = TRUE) %>%
  filter(n > 1)

x

dup_rows <- bat_clean_v2 %>%
  semi_join(x, by = "in_file") %>%
  arrange(in_file, datetime)

View(dup_rows)

# after inspection the duplicate rows in the file. It seem we can safely merge
# them into one row as they share the species and everything else

bat_clean_v3 <- bat_clean_v2 %>%
  distinct(in_file, sp, .keep_all = TRUE)

# how many rows we removed from one to another version (about 6244)

nrow(bat_clean_v2)
nrow(bat_clean_v3)

nrow(bat_clean_v2) - nrow(bat_clean_v3)

bat_clean_v3 # this files has acoustic, duplicates, and species vetted as much as possible. 


# count matrix ------------------------------------------------------------
# # this is a matrix where we create a n column that tells us how many calls for each bat are there.
# # daily counts
# 
# bm <- bat_clean_v2 %>% # 
#   group_by(noche, sp, site,yr, treatmt, trmt_bin, eff.hrs) %>% 
#   summarise(n = n(), .groups = 'drop') 
# 
# summary(bm)
# head(bm)

# we updated the previews section because it did not accounted for zeros. in other words when the species was not detected at a site but the site was sampled.

# 1. Summarize observed bat detections
bat_counts <- bat_clean_v3 %>%
  group_by(site, yr, jday, noche, sp) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )

# 2. Build effort grid
# Use effort_hrs if it contains all sampled site-nights.

effort_grid <- effort_hrs %>%
  ungroup() %>%
  distinct(site, yr, jday, noche, eff.hrs)

# 3. Species list
species_list <- bat_clean_v3 %>%
  distinct(sp)
# edit effort list to remove hif, mysp, lof

species_list <- species_list %>%
  filter(!sp %in% c("hif", "mysp", "lof"))

# 4. Create full site-night-species grid
bat_zero_db <- effort_grid %>%
  crossing(species_list)

# 5. Join observed counts and replace missing counts with zero
bat_zero_db <- bat_zero_db %>%
  left_join(
    bat_counts,
    by = c("site", "yr", "jday", "noche", "sp")
  ) %>%
  mutate(
    n = replace_na(n, 0)
  )

summary(bat_zero_db)


# calculate week for later mergin with insects

bat_zero_db <- bat_zero_db %>%
  mutate(
    wk = lubridate::week(noche)
  )

summary(bat_zero_db)

# here we calculate the total calls per species to report in the results 
bm_summary <- bat_zero_db %>%
  group_by(sp) %>%
  summarise(total_calls = sum(n)) %>%
  arrange(desc(total_calls))


# predictors --------------------------------------------------------------
# merge with predictors and standardize and see correlations 

# Craters weather (night)
crmo.wet.night <- read_csv("data_for_analysis/weather/craters_weater/craters_night.csv") %>% 
  clean_names() # load craters night weather

# Moon
moon.int <- read_csv('data_for_analysis/moon_pred/moon.int.csv') %>%
  clean_names()
summary(moon.int) # check the structure of the moon data

# convert date times from UTC to Denver/America
moon.int$denver_time <- with_tz(moon.int$date, tzone = "America/Denver") # Convert to Denver time zone
attr(moon.int$denver_time, "tzone") # check the timezone is correct

# summarize moonlight by date but conditional moon_alt_degrees > 0

# Step 1: Filter data where the moon is above the horizon
moon_filtered <- moon.int %>%
  filter(moon_alt_degrees > 0)

# Step 2: Create a new 'noche' variable (just the date part of the time stamp)
# but also makes nights any time stamps that are less than 9 am.
moon_filtered <- moon_filtered %>%
  mutate(
    hour = hour(denver_time),  # Extract hour from datetime
    noche = if_else(hour < 9,
                    true = as_date(denver_time) - days(1),
                    false = as_date(denver_time))
  )

# Step 3: Group by 'noche', then summarize the average values
moon_daily_avg <- moon_filtered %>%
  group_by(noche) %>%
  summarise(
    avg_moonlight = mean(moonlight_model, na.rm = TRUE),
    avg_twilight = mean(twilight_model, na.rm = TRUE),
    avg_illumination = mean(illumination, na.rm = TRUE),
    n_obs = n()  # optional: number of observations per night
  )

# insects 
c_bugs <- read_csv("data_for_analysis/insect_wranglin/c_bugs.csv") %>%  # load insect data
  clean_names() %>%
  rename(yr = yrs) # safe rename
# add treatment. 

litsites<-c("iron01","iron03","iron05","long01","long03")


c_bugs$treatmt<-ifelse(c_bugs$site %in% litsites , "lit", "dark") # this makes a treatment variable.


# calculate mean by yr, trmt and site. I will use this to substitute the NA valeus from tha appeared when merging the bat data. 

c_bugs_mean <- c_bugs %>%
  group_by(yr,treatmt, site) %>%
  summarise(
    t_insect = mean(t_insect, na.rm = TRUE),
    t_lepidoptera = mean(t_lepidoptera, na.rm = TRUE)
  ) %>%
  ungroup()



# light
light <- read_csv("data_for_analysis/lights/lightspectra_pioneer.csv") %>%
  clean_names() %>%
  filter(vert_horiz == "Horizontal") %>%                          # Keep only horizontal measures
  mutate(mwatts = rowMeans(across(c(watts_m1, watts_m2, watts_m3)), na.rm = TRUE)) %>%  # Mean of watts
  select(site, lux, yr, mwatts)    # Select relevant columns  

# Merge datasets ------------------------------------------------------------


# merge weather
# Merge crmo.wet.night into filtered_bm by matching dates
bm2 <- bat_zero_db %>%
  left_join(crmo.wet.night, by = c("noche" = "date"))
summary(bm2)

# merge with moon
bm2 <- bm2 %>%
  left_join(moon_daily_avg, by = "noche") 

summary(bm2) #there is just 78 NA I can keep it that way for now because we have to recalculate the moon predictors 

missing_moon_nights <- bat_zero_db %>%
  distinct(noche) %>%
  anti_join(
    moon_daily_avg %>% distinct(noche),
    by = "noche"
  )

missing_moon_nights

# merge with insects

bm2 <- bm2 %>%
  left_join(c_bugs, by = c("site", "wk", "yr")) # merge

# check for NAs
summary(bm2)  #lots of NAs in t.insect and t.lepidoptera.

# replace NAs in t.insect and t.lepidoptera with the mean values from c_bugs_mean


bm2 <- bm2 %>%
  left_join(c_bugs_mean, 
            by = c("yr", "site"),
            suffix = c("", "_mean")) %>%  # rename mean columns directly
  mutate(
    t_insect = coalesce(t_insect, t_insect_mean),
    t_lepidoptera = coalesce(t_lepidoptera, t_lepidoptera_mean)
  ) %>%
  select(-t_insect_mean, -t_lepidoptera_mean)

summary(bm2)

# merge light 

bm2<- bm2 %>%
  left_join(light, by = c("site", "yr"))


# we need to add the treatment
litsites<-c("iron01","iron03","iron05","long01","long03")

bm2<- bm2 %>%
  mutate(
    treatmt = if_else(site %in% litsites, "lit", "dark"),
    trmt_bin = if_else(treatmt == "lit", 1, -1)
  )

summary(bm2)
# note:
# we are going to write this table to use for the predictors in the dbrda analysis. 
# write_csv(bm2, "data_for_analysis/dbrda/bm2.csv")

# correlation -------------------------------------------------------------


# 1. Select numeric columns, optionally drop unique id/group columns
numeric_data <- bm2 %>% 
  select(where(is.numeric)) %>% 
  select(-any_of(c("yr", "trmt_bin"))) # add/remove columns as needed

# 2. Remove zero-variance columns
numeric_data <- numeric_data %>% select(where(~sd(., na.rm = TRUE) > 0))

# 3. Compute correlation matrix
cor_mat <- cor(numeric_data, use = "pairwise.complete.obs")

# 4. Visualize correlation matrix
corrplot(cor_mat, order = 'AOE', method = 'color', tl.col = 'black', tl.cex = 0.8, addCoef.col = 'black', number.cex = 0.7)


# standardize -------------------------------------------------------------

# #calculate Julian day from 'noche' and scale variables
# bm2<-bm2 %>% 
#   mutate(
#     jday = lubridate::yday(noche),  # Calculate Julian day from 'noche'
#   )

variables_to_scale <- c(
  "avg_moonlight",
  "avg_twilight",
  "avg_illumination",
  "nit_avg_temp_c",
  "nit_avg_wspm_s",
  "t_lepidoptera",
  "t_insect",
  "lux",
  "mwatts",
  "jday"
)

bm2 <- bm2 %>%
  scale_by_2sd_tidy(variables_to_scale)

summary(bm2)

# year standardize. 

# make year between -1:1
bm2 <- bm2 %>%
  mutate(yr_s = case_when(
    yr == 2021 ~ -1,
    yr == 2022 ~ 0,
    TRUE ~ 1
  ))


# species names for graphs

species <- data.frame(
  sp = c("ANTPAL", "CORTOW", "EPTFUS", "EUDMAC", "LASCIN", "LASNOC",
         "MYOCAL", "MYOCIL", "MYOEVO", "MYOLUC", "MYOTHY", "MYOVOL",
         "MYOYUM", "PARHES"),
  species_name = c("Antrozous pallidus", "Corynorhinus townsendii", "Eptesicus fuscus", "Euderma maculatum",
                   "Lasiurus cinereus", "Lasiurus noctivagans", "Myotis californicus", "Myotis ciliolabrum",
                   "Myotis evotis", "Myotis lucifugus", "Myotis thysanodes", "Myotis volans",
                   "Myotis yumanensis", "Parastrellus hesperus")
)

species <- species %>%
  mutate(
    sp = tolower(sp),  # Ensure species codes are lowercase for consistency)
    genus = word(species_name, 1),
    species = word(species_name, 2),
    sp_label = paste0(substr(genus, 1, 1), ".", species)
  )

# we make treatment bin -1 for dark and 1 for lit. 
bm2 <- bm2 %>%
  left_join(species %>% select(sp, sp_label), by = "sp") # add species labels for plotting

glimpse(bm2)
summary(bm2)

# now we merge the bm_ai with the bm2 data to have both the Miller activity index and the echolocation count data.

t <- bm2 %>%
  left_join(
    bm_miller %>%
      select(site, noche, sp, activity_min),
    by = c("site", "noche", "sp")
  ) %>%
  mutate(
    activity_min = replace_na(activity_min, 0L)
  )

summary(bm2)
# last check before export to the glmm_v5 script. lets see if theres any NA in the activity_min column.and if we need to remove noise or hif observations 

unique(bm2$sp)
tbm2 <- table(bm2$sp)


summary(bm2)
# check the rows that were removed to see if it is right.

# nrow(bm2) - nrow(t) # we removed  1, 2, and 3 species. we romved abotu 1017 rows that adds to the mysp+hif+lof calls.

# miller matrix -----------------------------------------------------------



# minutes activity 
# in here we calculate the minutes of activity inspired by Miller 2001 paper. 

bat_clean_v3 <- bat_clean_v3 %>%
  mutate(
    rmins = floor_date(datetime, unit = "minute")
  ) #rounds to the nearest min

bm_miller <- bat_clean_v3 %>%
  distinct(site, sp, noche, rmins) %>%
  group_by(site, sp, noche) %>%
  summarise(
    activity_min = n(),
    .groups = "drop"
  )

# bm_miller.day <- bm_miller %>% # number of minutes active  by night.  redundant might erase later. 
#   group_by(site, noche, sp) %>%
#   summarize(activity_min = sum(activity_min))
# 
# summary(bm_miller.day)

# 
# bm.miller<- bat_combined %>%
#   # Extract the relevant columns and round to minute level
#   mutate(rmins = round(date_time, units = "mins")) %>%
#   # Remove duplicate entries for the same site, date, and minute
#   distinct(site, sp, noche, rmins, .keep_all = TRUE) %>%
#   # Group by site, date, and minute
#   group_by(site,sp, noche, rmins) %>%
#   # Summarize to count the number of unique minutes
#   summarize(activity_min = n_distinct(rmins), .groups = 'drop')
# 
# bm.miller.day <- bm.miller %>% # number of minutes active  by night. 
#   group_by(site, noche, sp) %>%
#   summarize(activity_min = sum(activity_min))
# 
# head(bm.miller.day)
# summary(bm.miller.day)



# sm3 _buzz  --------------------------------------------------------------

# here we create a section to merge the sm3 buzz data with the bat_clean_v2 data. 

# load the sm3_buzz data. 

sm3_buzz <- read.csv("data_for_analysis/sm3_buzz_build/sm3_buzz_all.csv") 
glimpse(sm3_buzz)  

# standardize site.

unique(sm3_buzz$site)

sm3_buzz$site = ifelse(sm3_buzz$site %in% "lon03","long03", sm3_buzz$site)
sm3_buzz$site = ifelse(sm3_buzz$site %in% "viz02","vizc02", sm3_buzz$site)
sm3_buzz$site = ifelse(sm3_buzz$site %in% "viz04","vizc04", sm3_buzz$site)

unique(sm3_buzz$site)
# we can't trust monitoring night because there are a many missing from the sm3_buzz file. I have to review the script that created it firts to figure it out. For now we are going to create the noche column using the date time column 

# first we make date_time into a date and time column 
sm3_buzz<-sm3_buzz %>% mutate(
datetime = ymd_hms(date_time, tz = "America/Denver")
) %>% 
select(-c(date_time, monitoring_night)) # we don't need the date_time column anymore.

# now we create using the datetime column 
# we will use a noon-to noon boundary for the day. 

sm3_buzz<- sm3_buzz %>% 
  mutate(
    noche = as.Date(datetime - hours(12))
  )


summary(sm3_buzz)
# I believe the best idea is to merge the bat data with the sm3_buzz data by file.

# Keep the filename and species information from bat_clean_v3
bat_species_lookup <- bat_clean_v3 %>%
  select(filename = in_file, sp) %>%
  distinct()

# Add sp to sm3_buzz by matching filenames
sm3_buzz <- sm3_buzz %>%
  left_join(
    bat_species_lookup,
    by = "filename",
    relationship = "many-to-one"
  )

summary(sm3_buzz)

# how many sp were succesfully added
sm3_buzz %>%
  summarise(
    total_rows = n(),
    sp_added = sum(!is.na(sp)),
    sp_not_added = sum(is.na(sp)),
    percent_added = round(mean(!is.na(sp)) * 100, 2)
  )

summary(sm3_buzz)



# how many column have species information now
names(sm3_buzz)
# spp has all the posibble id we won't use it 
# spp_accp has all species that were accepted by the sonobat 
# species_manual_id has the vetted ones 
# sp the species ID with Kpro that  were added to the sm3_buzz data in the previous step. 

# clean up and create un single species column sp_clen 

# we create the sp_clean column where the priority is as follows sp_clean = species_manual_id > sp_accp > sp
# # now there's still some rows that have sp id in the spp column but this is composed of several option separated by slash myvol/myoluc/myocal. I want to take the firs option only. 

sm3_buzz <- sm3_buzz %>%
  mutate(
    sp_clean = case_when(
      !is.na(species_manual_id) ~ species_manual_id,
      !is.na(spp_accp)          ~ spp_accp,
      !is.na(sp)                ~ sp,
      !is.na(spp)               ~ sub("/.*$", "", spp),
      TRUE                      ~ NA_character_
    )
  )


summary(sm3_buzz)

# I check and this worked. Now the sp_clean column has the species ID that we can use for analysis, but needs to be cleaned first from 6code ID to 4 code so myovol should become myvo.

sm3_buzz %>%
  select(species_manual_id, spp_accp, sp, spp, sp_clean) %>%
  head(30)

# fix the sp_name in sp_clean to be 4 letter code instead of 6 letter code.
sm3_buzz <- sm3_buzz %>%
  mutate(
    sp_clean = tolower(sp_clean),
    
    sp_clean = if_else(
      nchar(sp_clean) == 6,
      paste0(
        substr(sp_clean, 1, 2),  # first 2 genus letters
        substr(sp_clean, 4, 5)   # first 2 species letters
      ),
      sp_clean
    )
  )

# now we have to decide what observations we keep, zeros with sp ID or any zeros. I will keep zeros associated with species. I also want to see how many feeding buzz above zero have no species ID associated to them 

# first I want to know what species have buzzes but no ID. the opposite. 

buzz_id_summary <- sm3_buzz %>%
  mutate(
    category = case_when(
      is.na(sp_clean) & c_buzz > 0  ~ "Buzzes but no species ID",
      !is.na(sp_clean) & c_buzz == 0 ~ "Species ID but no buzzes",
      !is.na(sp_clean) & c_buzz > 0  ~ "Species ID and buzzes",
      is.na(sp_clean) & c_buzz == 0  ~ "No species ID and no buzzes",
      is.na(c_buzz)                  ~ "Missing buzz count"
    )
  ) %>%
  count(category, name = "observations") %>%
  mutate(
    percent = round(100 * observations / sum(observations), 2)
  )

buzz_id_summary

# we filter out the rows with NA in sp_clen and c_buzz. 
sm3_buzz<- sm3_buzz %>%
  filter(
    !is.na(sp_clean),
    sp_clean != "", # possibly unnecessary
    !is.na(c_buzz)
  )

# tally/count
sm3_buzz %>%
  count(sp_clean, sort = TRUE)

summary(sm3_buzz)
glimpse(sm3_buzz)

# now we summarize bat buzzes by species site and monitoring night and year. If we don't do this then we have several buzzes counted for the same species in a site on a given monitoring night. 
# 
sm3_buzz <- sm3_buzz %>%
  group_by(
    site,
    noche,# we can't use monitoring night thus noche
    year,
    sp_clean
  ) %>%
  summarise(
    t_buzzes = sum(c_buzz),
    identified_files = n(),
    .groups = "drop"
  )

summary(sm3_buzz)


# now we want to add the effort from the table 

# # join effort with sm3_buzz. effort is leaving several NAs and I think we have some issues with this. 
# 



sm3_buzz <- sm3_buzz %>% 
  left_join(
    effort_lookup,
    by = c("site", "year", "noche")  )
summary(sm3_buzz) # as it stands there's 32 sites that do not have effort data. I will need to check this ones for now, I will give them the minimum effort of 1 hour.

sm3_buzz <- sm3_buzz %>%
  mutate(
    eff.hrs = if_else(is.na(eff.hrs), 1, eff.hrs),
    eff.days = if_else(is.na(eff.days), 1, eff.days)
  )

summary(sm3_buzz) # now we have no NAs in the effort column.

# now we add the zeros 



# outputs -----------------------------------------------------------------


# dir.create("data_for_analysis/prep_for_glm", showWarnings = FALSE) # just run if the dir is abscent

write.csv(bat_combined, file = 'data_for_analysis/prep_for_glmm_v2//bat_combined.csv', row.names = F) # raw combine data 
write.csv(bm2, file = 'data_for_analysis/prep_for_glmm_v2/bm2.csv', row.names = F) #daily counts
write.csv(bm_miller, file = "data_for_analysis/prep_for_glmm_v2/bm_miller.csv") # miller Ai index data
# write.csv(sm3_buzz, file = "data_for_analysis/prep_for_glmm/sm3_buzz_sp.csv") # this is not ready to write 


# Create a README file with information about the script
readme_content <- "

Carlos Linares 8/01/2024, 12/12/2024, 7/8/2025 
This directory contains the bat_combined.csv file which was created using the script prep_for_glmm_v2.R combines bat species call abundance data. This script merges 2021-23 data that was previously scanned with Kaleidoscope pro. but in this new version we created a series of rules to re-code and flagg species that are not likely to be found in the study area (bat_clean_v2). 

- The script also creates a count matrix bat_zero_db taht has the summary of calls by species with zeros added. These represent times when a bat was not detected but the site was monitored. 

the script also produces a Miller 2001 index of activity matrix bm_miller.

bat_combined.csv - process data no counts (Update: 7/8/2026 rules for bat species) 
bm2.csv - bats counts ready for analysis with zero added, predictors and standardized 
bm_miller.csv - number of minutes of activity by day  for 2021-2023 data all sites (last update 7/8/2026)


"
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
