# =================================================
# = build database from Sonobat output files  =
# =================================================
# the objective is to build a data base with the individual files that Sonobat creates. 
# these are saved in different folders depending on the year, site, instrument. 
# 
# 
# 
# Carlos Linares  12/19/2023

# libraries ---------------
library(data.table)
library(purrr)
library(stringr)


# path to the z drive. 


setwd("Z:/PioneerLights_2021") #you have to set the wd() to the place where the files are stored

load("database_workspace.RData")

patt<- "Z:/PioneerLights_2021/"

# paths list

myfiles_long <- list.files(path = patt, recursive = TRUE,pattern = "*.txt", full.names = T) # thjs takes a long time. 
myfiles <-fs::dir_ls(path = patt, recurse = T, glob = "*.txt", ) # 695 this one runs faster than the above

# I just checked and both lines of code have more or less the same amount of paths.


# sm3 files

sm3dirs<- myfiles[grepl("sm3", myfiles)] # this does it for just the sm3 files. 

# Filter paths  

filtered_paths <- sm3dirs[grepl("CumulativeSonoBatch.*(_v4\\.4\\.5\\.txt|_v420\\.txt)$", sm3dirs)] # this grabs both sonofiles outputs v.4.4.5 and v420

paths_v4.4.5<- sm3dirs[grepl("CumulativeSonoBatch.*(_v4\\.4\\.5\\.txt)$",sm3dirs)] #adds to 6
pats_v420<- sm3dirs[grepl("CumulativeSonoBatch.*(_v420.txt)$",sm3dirs)] # adds to 82 both add to 88 the total in line 88



# t<-setdiff(as.vector(sonofiles1), as.vector(filtered_paths))



# Merging all 

dlist<- sapply(filtered_paths, read.delim, simplify = F, fill=T) #read all the files simplify assure they are returned as a list of dataframes. 


# Rbind efficiently, handling potential column differences
bat2021 <- data.table::rbindlist(dlist, fill = TRUE, idcol = T)


write.csv(x = bat2021, file ="bats2021_v5.csv")


# load old data base and compare

oldbats<-read.csv("bats2021.csv")



#compare old and new data base

nuevo<-bat2021$Path
viejo<-oldbats$Path

a<-setdiff(nuevo,viejo) #they are definetively different

# let's see if the new data base has all sites. 

bat2021$site<-str_extract(bat2021$Filename, "^[A-Za-z]{3,4}\\d{2}")

unique(bat2021$site)

#I want to know how many data points we have for each site. 




save.image(file = "scripts/database_workspace.RData")# this saves the workspace so I don't have to run this every time


# robomoth db -------------------------------------------------------------


sm4dirs<- myfiles[grepl("sm4",myfiles)] # 88 txt files

sm4sono<- sm4dirs[grepl("_CumulativeSonoBatch_v\\d+\\.\\d+\\.\\d+\\.txt", sm4dirs)]
sm4sono_v420<- sm4dirs[grepl("CumulativeSonoBatch_v420.txt", sm4dirs)] # files analyzed with the older version of sonobat

#list kpro files 
myfiles <-fs::dir_ls(path = patt, recurse = T, glob = "*.csv", ) # 695 this one runs faster than the above

kproout<- sm4dirs[grepl("^id.csv", sm4dirs)]


sm4sonomerge<- sapply(sm4sono, read.delim, simplify = F) #read all the files 
sm4sonomerge_v420<-sapply(sm4sono_v420, read.delim, simplify = F)
sm4_kpromerge<- sapply(kpro, read.delim,simplify = F)

robo2021 <- data.table::rbindlist(sm4sonomerge)
robo2021_v420 <- data.table::rbindlist(sm4sonomerge_v420)
robo2021_kpro<-data.table::rbindlist(sm4_kpromerge)

# check both databases have the same cols
columns_df1 <- colnames(robo2021)
columns_df2 <- colnames(robo2021_v420)

# Check if the column names are the same
if (identical(columns_df1, columns_df2)) {
  print("Both databases have the same columns.")
} else {
  print("The databases have different columns.")
}

# merge data bases of data analyzes with v4.4.5 and v420 of sonobat.

# Add a column to each dataframe indicating the source
robo2021 <- mutate(robo2021, source = "v.4.4.5")
robo2021_v420 <- mutate(robo2021_v420, source = "v420")

# keep cols
col_keep <-
  c(
    "Path",
    "Filename",
    "HiF",
    "LoF",
    "SppAccp",
    "Prob",
    "calls.sec",
    "X1st",
    "X.Spp",                                                                
    "X.Prob",
    "source"
  )




# Merge and create a new dataframe
merged_df <- full_join(robo2021, robo2021_v420, by=col_keep) %>% 
  select(all_of(col_keep))

write.csv(x = merged_df, file = "robo2021.csv",)




# Old code below----------------------

# files2process<-list.files(filessm3, recursive = T, pattern = "Data_SonoBatch_v420.txt", full.names = T)


#test with listing dir instead of files seems faster. 

t<- list.dirs(path = ".", recursive = T) # list all the files  in the z folders for pioneer lights 2021

sm3dirs<- t[grepl("sm3", t)] # this does it for just the sm3 files. 

# t2<-sm3dirs[1:100] # this was just a test with 100 rows. 

files2process<- list.files(sm3dirs, recursive = T, pattern = "Data_SonoBatch_v420.txt", full.names = T)

dlist<- sapply(files2process, read.delim, simplify = F) #read all the files 

bat2021 <- data.table::rbindlist(dlist)
bat2021_update<-bat2021

# now we save the data on a different file but be aware that this needs to be updated once the all the data has been ran on sonobat. In december 2023 I have been trying to update this file. 

write.csv(bat2021_update,file = 'bats2021_update.csv')

# bat2021 <- data.table::rbindlist(dlist)


# lets clean the data base in a different script.

write.csv(bat2021_update,file = 'bats2021_update.csv') # we need to check this update and how it compare to the first one.

# Now I want the data from 2022

# list the files in a directory

#test directory is iron 01 
tp <-
  list.files(
    'PioneerLights_2022/sm3/iron01/08052022/',
    recursive = T,
    pattern = ".txt$",
    full.names = T
  )

ls2023 <-
  list.files(
    'PioneerLights_2022/sm3/',
    recursive = T,
    pattern = ".txt$",
    full.names = T
  )

sbat<- ls2023[grepl("CumulativeSonoBatch_v\\d{3}\\.txt", ls2023)]

dlist<- sapply(sbat, read.delim, simplify = F)

bat2022 <- data.table::rbindlist(dlist)



# sonofiles1<- sm3dirs[grepl("CumulativeSonoBatch_", sm3dirs)] # this grabs more than we want. 
# sonofiles2<- sm3dirs[grepl("SonoBatch_v420.txt", sm3dirs)]