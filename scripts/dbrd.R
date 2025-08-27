# ---------------------------
##
## Script name:  bat_diversity_v2.R
##
## Purpose
#  We want to assess if the bat community changes with the light treatment. 
# 
## Author: Carlos Linares, 
## Date Created:8/25/2025
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes: Bat data product of the prep_for_glmm.R script
##        analysis following the methods from https://esajournals-onlinelibrary-wiley-com.libproxy.boisestate.edu/share/KSJBJ6ZJGCCGVC4TBJPR?target=10.1002/ecy.70128
##   
##
## ---------------------------
## # inputs ------------------------------------------------------------------
# - bm<- read_csv("data_for_analysis/prep_for_glmm/bm.csv")
#

# outputs ----------------------

# plots comparing activity between sites 


# libraries  --------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  "tidyverse",
  "lubridate",
  "here",# for reproducible file paths
  "janitor",
  "purrr",
  "patchwork",
  "vegan"
)



# data --------------------------------------------------------------------

# bat data 
bm<- read_csv("data_for_analysis/prep_for_glmm/bm.csv") %>% 
  clean_names()

bat_combined<-read_csv("data_for_analysis/prep_for_glmm/bat_combined.csv")

#effort

effort_days <- bat_combined %>%
  group_by(site, yr) %>%
  summarise(
    stard = min(noche),
    endd = max(noche),
    eff.days = as.numeric(difftime(max(noche), min(noche), units = "days"))
  )

total_effort <- effort_days %>%
  group_by(site) %>%
  summarise(total_effort = sum(eff.days, na.rm = TRUE))


# predictors --------------------------------------------------------------


# light data 

light<- read_csv("data_for_analysis/lights/lightspectra_pioneer.csv") %>% 
  clean_names() %>% 
  filter(vert_horiz == "Horizontal") %>%  # select just horizontal measures
  select("site","vert_horiz","watts_m1","watts_m2","watts_m3","lux","yr") %>% # keep this cols
  mutate(mwatts = rowMeans(across(c(watts_m1, watts_m2, watts_m3)), na.rm = TRUE)) # calculate mean watts. 

# calculate mean light by site

mlight <- light %>%
  group_by(site) %>%
  summarise(mean_mwatts = mean(mwatts, na.rm = TRUE))


# elevation

elev<-read_csv('data_for_analysis/elev/elevation.csv', name_repair = "universal")
elev<-rename(elev, site=name)

# insects

insect<-read_csv('data_for_analysis/insect_wranglin/c_bugs.csv', name_repair = "universal")

# we need to summarize total lepidoptera by site 

insect<-insect %>% 
  group_by(site) %>% 
  summarise(t.leps= sum(t.lepidoptera))

# check predictors for colinearity

cor(pred_matrix %>% select(mean_mwatts, t.leps, elev_mean), use = "pairwise.complete.obs")



# matrix build ------------------------------------------------------------



# We need a site × species matrix with abundances.

bat_comm <- bm %>%
  group_by(site, auto_id) %>%     # site = sampling unit, auto_id = species
  summarise(abundance = sum(n), .groups = "drop") %>%
  pivot_wider(names_from = auto_id, values_from = abundance, values_fill = 0) %>% 
  select(-c(NoID, Noise)) # remove Noise and NoID

# Check matrix
head(bat_comm)

comm_matrix <- bat_comm %>% # make column site into rownames 
  arrange(site) %>%
  column_to_rownames("site")

# standardize and divide by number of days sampled after 3 years. 

comm_matrix_std <- comm_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("site") %>%
  left_join(total_effort, by = "site") %>%
  mutate(across(-c(site, total_effort), ~ .x / total_effort)) %>%
  select(-total_effort) %>%
  tibble::column_to_rownames("site")

comm_matrix_std

# Create Bray-Curtis distance matrix
bray_dist <- vegdist(comm_matrix_std, method = "bray")



# predictor matrix

pred_matrix<-left_join(mlight,insect, by="site")
pred_matrix<-left_join(pred_matrix, elev, by="site")
pred_matrix<-pred_matrix %>% select(-c(time, buff_area))

# standardize predictors

pred_matrix_std <- pred_matrix %>%
  mutate(across(
    c(mean_mwatts, t.leps, elev_mean),
    ~ (.-mean(. , na.rm = TRUE)) / (2*sd(. , na.rm = TRUE))
  ))

# Run partial dbRDA---------------------------------------------------

# # just light predictor
# dbrda_full <- dbrda(bray_dist ~ mean_mwatts, data = pred_matrix)
# 
# # full dbRDA 
# dbrda_partial <- dbrda(bray_dist ~ mean_mwatts + t.leps + elev_mean , data = pred_matrix_std)

dbrda_full<- dbrda(comm_matrix_std ~ mean_mwatts + t.leps + elev_mean, pred_matrix_std,dist="bray")
# add species scores 

sppscores(dbrda_full)<-wisconsin(bray_dist)




#Assess model significance ---------------------------------------------------

# Test marginal (type III) effects of each predictor
anova_marginal <- anova.cca(dbrda_full, by = "margin", permutations = 999)
anova_marginal


plot(dbrda_full, display = c("sites", "bp"), scaling = 2)

library(ggvegan)

autoplot(dbrda_full, scaling = 2) +
  theme_minimal() +
  labs(title = "dbRDA of bat communities",
       subtitle = "Constrained by mean_mwatts, t.leps, and elev_mean")


# Get site, species, and environmental scores
site_scores <- scores(dbrda_full, display = "sites", scaling = 2)
site_scores_df <- as.data.frame(site_scores)
site_scores_df$site <- rownames(site_scores)

species_scores <- scores(dbrda_full, display = "species", scaling = 2)
species_scores_df <- as.data.frame(species_scores)
species_scores_df$species <- rownames(species_scores)

env_scores <- scores(dbrda_full, display = "bp", scaling = 2)
env_scores_df <- as.data.frame(env_scores)
env_scores_df$variable <- rownames(env_scores)

# Plot
p1<-ggplot() +
  # Sites (small, light gray)
  geom_point(data = site_scores_df, 
             aes(x = dbRDA1, y = dbRDA2), 
             color = "gray60", size = 2, alpha = 0.6, shape = 16) +
  geom_text_repel(data = site_scores_df, 
                  aes(x = dbRDA1, y = dbRDA2, label = site),
                  size = 3, color = "gray60",
                  max.overlaps = 30) +
  
  # Species (blue, bigger font)
  geom_text_repel(data = species_scores_df, 
                  aes(x = dbRDA1, y = dbRDA2, label = species),
                  size = 4, fontface = "bold", color = "black",
                  max.overlaps = 30) +
  
  # Environmental vectors (arrows, red)
  geom_segment(data = env_scores_df,
               aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               arrow = arrow(length = unit(0.25, "cm")), 
               color = "#D55E00", size = 1.2) +
  geom_text_repel(data = env_scores_df, 
                  aes(x = dbRDA1, y = dbRDA2, label = variable),
                  size = 4, color = "#D55E00",
                  max.overlaps = 30) +
  
  labs(title = "dbRDA of bat communities",
       x = "dbRDA1", y = "dbRDA2") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
p1

ggsave("figures/dbrd/dbRDA_v1.tiff", p1, width = 10, height = 8)
