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

# add light as a treatment 
litsites<-c("iron01","iron03","iron05","long01","long03")


mlight$treatmt<-ifelse(mlight$site %in% litsites , "lit", "dark")

head(mlight)
# elevation

elev<-read_csv('data_for_analysis/elev/elevation.csv', name_repair = "universal")
elev<-rename(elev, site=name)

# insects

insect<-read_csv('data_for_analysis/insect_wranglin/ins_bm.csv', name_repair = "universal")

# we need to summarize total lepidoptera by site 

insect<-insect %>% 
  group_by(site) %>% 
  summarise(t.leps= sum(t.lepidoptera))



# matrix build ------------------------------------------------------------



# We need a site × species matrix with abundances.

bat_comm <- bm %>%
  group_by(site, sp) %>%     # site = sampling unit, auto_id = species
  summarise(abundance = sum(n), .groups = "drop") %>%
  pivot_wider(names_from = sp, values_from = abundance, values_fill = 0) %>% 
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

head(bray_dist)

# predictor matrix

pred_matrix<-left_join(mlight,insect, by="site")
pred_matrix<-left_join(pred_matrix, elev, by="site")
pred_matrix<-pred_matrix %>% select(-c(time, buff_area))

# standardize predictors Gelman method

pred_matrix_std <- pred_matrix %>%
  mutate(across(
    c(mean_mwatts, t.leps, elev_mean),
    ~ (.-mean(. , na.rm = TRUE)) / (2*sd(. , na.rm = TRUE))
  ))


head(pred_matrix_std)

# check predictors for colinearity

cor(pred_matrix %>% select(mean_mwatts, t.leps, elev_mean), use = "pairwise.complete.obs")



# Run partial dbRDA---------------------------------------------------

# # just light predictor
 dbrda_full <- dbrda(bray_dist ~ treatmt, data = pred_matrix)

# 
# # full dbRDA 
# dbrda_partial <- dbrda(bray_dist ~ mean_mwatts + t.leps + elev_mean , data = pred_matrix_std)

dbrda_full<- dbrda(comm_matrix_std ~ mean_mwatts + t.leps + elev_mean, pred_matrix_std,dist="bray")

# run light as treatment 

dbrda_full<- dbrda(comm_matrix_std ~ treatmt + t.leps + elev_mean, pred_matrix_std,dist="bray")

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

library(ggrepel)

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

ggsave("figures/dbrd/dbRDA_v1.tiff",
       p1,
       dev = "tiff",
       dpi = 600,
       width = 10,
       height = 8,
       bg = "white",)



# dbrda second version  ---------------------------------------------------
# here we want to make the site-year the sampling unit and see if there's any signal for light as treatment not as watts affecting the community. 

# sampling unit year site 

row_unit<-c("site", "yr")

# effort: nights per row unit 

effrot_tbl<- bm %>%
  distinct(across(c(all_of(row_unit), noche))) %>% #unique night sampled
  count(across(all_of(row_unit)), name = "n_nights") # count nights per row unitross)

effort_keep<-effrot_tbl %>%  filter (n_nights > 5) # we keep those nights site years with more than 5 nights sampled.

# species abundance per row unit

comm_long<- bm %>% 
  semi_join(effort_keep, by = row_unit) %>% 
  group_by(across(c(all_of(row_unit), sp))) %>%
  summarise(n_tot= sum(n, na.rm= TRUE), .groups = "drop")

#remove noise and noid
comm_long<- comm_long %>% filter(!sp %in% c("Noise", "NoID")) # remove noise and noid

# standardize by effort

comm_std_long<- comm_long %>% 
  left_join(effort_keep, by = row_unit) %>% 
  mutate(activity = n_tot/n_nights)

#wide community matrix: one row per row unit, colums = species 

bat_com_std<-comm_std_long %>% 
  select(all_of(row_unit), sp, activity) %>% 
  pivot_wider(names_from = sp, values_from = activity, values_fill = 0) %>% 
  arrange(across(all_of(row_unit))) 

# Save the row id (for later joins) and set rownames for vegan
row_id <- bat_com_std %>% select(all_of(row_unit))
comm_matrix_std <- bat_com_std %>%
  select(-all_of(row_unit)) %>%
  as.data.frame()
rownames(comm_matrix_std) <- apply(row_id, 1, paste, collapse = "_")


# com matrix for mvGLM site-year x species count

bat_comm_counts <-comm_std_long %>% 
  select(all_of(row_unit), sp, n_tot) %>% 
  pivot_wider(names_from = sp, values_from = n_tot, values_fill = 0) %>% 
  arrange(across(all_of(row_unit)))

row_id <- bat_comm_counts %>% select(all_of(row_unit))

comm_matrix_counts <- bat_comm_counts %>%  # this one should be used in mvglm
  select(-all_of(row_unit)) %>%
  as.data.frame()
rownames(comm_matrix_counts) <- apply(row_id, 1, paste, collapse = "_")




# Bray–Curtis distance ----------------------------------------------
# we take the com matrix standardize to caculate the brady index matrix
bray_dist <- vegdist(comm_matrix_std, method = "bray")

# now we need the predictors at the same resolution as the comm matrix, 

# we read the bm2 data used in the glmmm_v4 analysis because it has all the predictors
bm2<- read_csv("data_for_analysis/dbrda/bm2.csv") %>% 
  clean_names() %>% 
  filter(!sp %in% c("Noise", "NoID")) # remove noise and noid

head(bm2)

# we need to add elevation 

bm2<- left_join(bm2, elev, by = "site")

# Example: mean light and Lepidoptera summed/averaged per row unit (site/year)

pred_matrix <- bm2 %>%
  group_by(across(all_of(row_unit))) %>%
  summarise(
    mean_mwatts = mean(mwatts, na.rm = TRUE),
    t_leps      = sum(t_lepidoptera, na.rm = TRUE),
    elev_mean   = mean(elev_mean, na.rm = TRUE), # lets try it without elevation for now
    .groups = "drop"
  ) %>%
  semi_join(effort_keep, by = row_unit)

# Keep only rows present in the community matrix
pred_matrix$RowKey <- apply(pred_matrix[row_unit], 1, paste, collapse = "_") # creates col Rowkey
pred_matrix <- pred_matrix %>% filter(RowKey %in% rownames(comm_matrix_std)) %>%
  column_to_rownames("RowKey")


#then we scale those variables

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


variables_to_scale <- c(
  "t.lepidoptera",
  "mean_mwatts",
  "elev_mean"
  )

pred_matrix_std <- pred_matrix %>%
  ungroup() %>%                                   # <- key
  scale_by_2sd_tidy(c("t_leps","mean_mwatts","elev_mean" ))

summary(pred_matrix_std)

cor(pred_matrix_std %>% select(mean_mwatts, t_leps, elev_mean), use = "pairwise.complete.obs")

# add treatment 

pred_matrix_std$treatmt<-ifelse(pred_matrix_std$mean_mwatts_s > 0, "lit", "dark")



# now we run the dbrda

dbrda_full <- dbrda(comm_matrix_std ~ treatmt + t_leps_s + elev_mean_s  ,
                    data = pred_matrix_std, distance = "bray")
dbrda_light<- dbrda(comm_matrix_std ~ treatmt, pred_matrix_std, dist= "bray")

# add species scores 
# dbRDA has no information on species you we have to include it manually.
sppscores(dbrda_full)<-wisconsin(bray_dist)

anova.cca(dbrda_full,dbrda_light)

# Test marginal (type III) effects of each predictor
anova_marginal <- anova.cca(dbrda_full, by = "margin", permutations = 999)
anova_marginal

anova.cca(dbrda_light)                # overall
anova.cca(dbrda_light, by = "margin") # marginal (type III) effects

# explained variation 

RsquareAdj(dbrda_full) # total variation explaned is r.quared 
RsquareAdj(dbrda_light)


# Variation partitioning
# 1. Apply Hellinger transformation to community data
comm_hell <- decostand(comm_matrix_std, method = "hellinger")

# 2. Run variation partitioning on transformed data
varpart_result <- varpart(comm_hell,
                          ~ mean_mwatts,
                          ~ t_leps,
                          ~ elev_mean,
                          data = pred_matrix)

# 3. Check the adjusted R² fractions
varpart_result
plot(varpart_result,
     bg = c("gold", "lightblue", "forestgreen"),
     Xnames = c("Light", "Lepidoptera", "Elevation"))


plot(dbrda_full, display = c("sites", "bp"), scaling = 2)
plot(dbrda_light)

library(ggvegan)

autoplot(dbrda_full, scaling = 2) +
  theme_minimal() +
  labs(title = "dbRDA of bat communities",
       subtitle = "Constrained by mean_mwatts, t.leps, and elev_mean")

autoplot(dbrda_light, scaling = 2) +
  theme_minimal() +
  labs(title = "dbRDA of bat communities",
       subtitle = "Constrained by mean_mwatts, t.leps, and elev_mean")


library(ggrepel)
# Get site, species, and environmental scores
site_scores <- scores(dbrda_full, display = "sites", scaling = 2)
site_scores_df <- as.data.frame(site_scores)
site_scores_df$site <- rownames(site_scores)

# compute species scores manually
species_scores <- wascores(site_scores, comm_matrix_std)
species_scores_df <- as.data.frame(species_scores)
species_scores_df$species <- rownames(species_scores_df)

env_scores <- scores(dbrda_full, display = "bp", scaling = 2)
env_scores_df <- as.data.frame(env_scores)
env_scores_df$variable <- rownames(env_scores)

# adding light
# make rows names to columns firt 
pred_matrix_std2 <- pred_matrix_std %>% rownames_to_column(var = "ID")
# do the same for stie_scores
site_scores_df<- site_scores_df %>% rownames_to_column(var = "ID")
# join them 
site_scores_df<- left_join(pred_matrix_std2, site_scores_df, by= "ID")

# # change names of the variables to be able to plot the actual names. 
# 
# site_scores_df<- site_scores_df %>% 
#   rename(light = mean_mwatts, elevation = elev_mean, leps = t_leps)

# Plot
p1<-ggplot() +
  # Sites (small, light gray)
  geom_point(data = site_scores_df, 
             aes(x = dbRDA1, y = dbRDA2, color = treatmt ), 
             size = 4, alpha = 0.6) +
  # geom_text_repel(data = site_scores_df,
  #                 aes(x = dbRDA1, y = dbRDA2, label = site.x, color = mean_mwatts),
  #                 size = 4, max.overlaps = 30) +

  # Species (blue, bigger font)
  geom_text_repel(data = species_scores_df,
                  aes(x = dbRDA1, y = dbRDA2, label = species),
                  size = 4, fontface = "bold", color = "black",
                  max.overlaps = 30, alpha=0.5) +

  # Environmental vectors (arrows, red)
  geom_segment(data = env_scores_df,
               aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               arrow = arrow(length = unit(0.25, "cm")), 
               color = "#D55E00", size = 1.2) +
  geom_text_repel(data = env_scores_df, 
                  aes(x = dbRDA1, y = dbRDA2, label = variable),
                  size = 6, color = "black",
                  max.overlaps = 30) +
  
  labs(title = "dbRDA of bat community",
       x = "dbRDA1", y = "dbRDA2") +
  scale_color_manual(values = c("dark" = "black", "lit" = "grey50"))+
  # scale_color_viridis_c(option = "inferno", end = 0.9) +  # nice continuous color scale
  theme_minimal(base_size = 14, base_family = "serif") +
  theme(legend.position = "right",
  )
p1

species_scores_df %>%
  arrange(desc(abs(dbRDA1))) %>%
  head(10)

# mvGLM -----------------------------------------------------------------

library(mvabund)
# Prepare data -------------------------------------------------------------
# Response = bat community abundance matrix (species × sites)
# Convert comm_matrix_std into an mvabund object
Y <- mvabund(comm_matrix)

# Predictors = environmental variables and offset oftotal effort. 
pred_matrix<- pred_matrix %>% left_join(total_effort) %>% 
  mutate(log_nights= log(total_effort))

X <- pred_matrix %>%
  column_to_rownames("site") %>%
  mutate(across(everything(), scale))   # center & scale predictors

#add treatment
X$treatmt<-ifelse(X$mean_mwatts > 0, "lit", "dark")

# Ensure matching row order
X <- X[rownames(comm_matrix), ]

# Build mvGLM --------------------------------------------------------------
# Negative binomial family is often best for overdispersed count data
fit <- manyglm(Y ~ mean_mwatts + t.leps + elev_mean, offset("log_nights"),
               data = X,
               family = "negative.binomial")

# Model diagnostics --------------------------------------------------------
# Residual checks (optional but recommended)
plot(fit)   # residual vs fitted, QQ plot etc.

# Hypothesis tests ---------------------------------------------------------
# Overall test: does the full predictor set explain variation in the community?
anova_full <- anova.manyglm(fit, resamp = "pit.trap", nBoot = 999)
anova_full

# Marginal tests: importance of each predictor individually (Type III)
anova_marginal <- anova.manyglm(fit, p.uni = "adjusted", resamp = "pit.trap", nBoot = 999)
anova_marginal

# Visualize effects --------------------------------------------------------
# Turn fitted values into long format for plotting
fitted_vals <- fitted(fit)
fitted_df <- as.data.frame(fitted_vals)
fitted_df$site <- rownames(comm_matrix)
fitted_long <- pivot_longer(fitted_df, -site, names_to = "species", values_to = "fit")

ggplot(fitted_long, aes(x = species, y = fit)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Fitted bat abundances across species",
       y = "Fitted abundance (calls)", x = "Species")




# mvGLM updated-------------------------------------------------------------------


# Prepare data -------------------------------------------------------------
# Response = bat community abundance matrix (species × site/year)
# Convert comm_matrix_std into an mvabund object
Y <- mvabund(comm_matrix_counts)

# Predictors = environmental variables and offset of total effort. 
pred_matrix_std <- pred_matrix_std %>%
  left_join(effort_keep, by = c("site", "yr")) %>%
  mutate(log_nights = log(n_nights))

X <- pred_matrix_std 
# mutate(across(everything(), scale))   # center & scale predictors

# # Ensure matching row order
# X <- X[rownames(comm_matrix), ]


# Build mvGLM --------------------------------------------------------------
# Negative binomial family is often best for overdispersed count data
fit <- manyglm(Y ~ treatmt + t_leps_s + elev_mean_s,
               data = X,
               family = "negative.binomial",
               offset = X$log_nights)
summary(fit)
# Model diagnostics --------------------------------------------------------
# Residual checks (optional but recommended)
plot(fit)   # residual vs fitted, QQ plot etc.

# Hypothesis tests ---------------------------------------------------------
# Overall test: does the full predictor set explain variation in the community?
anova_full <- anova.manyglm(fit, resamp = "pit.trap", nBoot = 999)
anova_full

# Marginal tests: importance of each predictor individually (Type III)
anova_marginal <- anova.manyglm(fit, p.uni = "adjusted", resamp = "pit.trap", nBoot = 999)
anova_marginal

anova_species <- anova.manyglm(
  fit,
  resamp = "pit.trap",
  nBoot = 999,
  p.uni = "adjusted"  # gives species-level tests
)



mv_tab_manual <- tibble::tibble(
  term = c("Treatment", "Lepidoptera", "Elevation"),
  Df = c(1, 1, 1),
  Dev = c(85.38, 18.82, 41.66),
  p_value = c(0.001, 0.251, 0.020)
)

flextable::flextable(mv_tab_manual)
# Extract univariate test table
uni_table <- as.data.frame(anova_species$uni.test)



# --- Clean and reshape univariate results for plotting ---

uni_table <- as.data.frame(anova_species$uni.test) %>%
  tibble::rownames_to_column("predictor")

uni_long <- uni_table %>%
  pivot_longer(
    cols = -predictor,
    names_to = "species",
    values_to = "p_value"
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )


# filter the intercept 

uni_long <- uni_long %>%
  filter(predictor != "intercept")


ggplot(uni_long, aes(x = predictor, y = reorder(species, p_value), color = p_value)) +
  geom_point(size = 4) +
  geom_text(aes(label = significance), nudge_x = 0.2, size = 4) +
  scale_color_viridis_c(option = "C", direction = 1, name = "p-value") +
  scale_x_discrete(
    labels = c(
      "mean_mwatts_s" = "Light (Watts)",
      "t_leps_s" = "Lepidoptera",
      "elev_mean_s" = "Elevation"
    )
  ) +
  labs(
    x = "Predictor",
    y = "Species",
    title = "Species-level significance (mvGLM)",
    subtitle = "P-values per predictor (PIT-trap resampling, 999 iterations)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(family = "Times", size = 12, face = "italic"),
    axis.text.x = element_text(family = "Times", size = 13, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

 # Visualize effects --------------------------------------------------------
# Turn fitted values into long format for plotting
fitted_vals <- fitted(fit)
fitted_df <- as.data.frame(fitted_vals)
fitted_df$site <- rownames(comm_matrix_counts)
fitted_long <- pivot_longer(fitted_df, -site, names_to = "species", values_to = "fit")

ggplot(fitted_long, aes(x = species, y = fit)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Fitted bird abundances across species",
       y = "Fitted abundance (calls)", x = "Species")

# rank abundance curves ---------------------------------------------------



# Calculate sampling effort per site and year
effort <- bm %>%
  group_by(site, yr) %>%
  summarise(nights = n_distinct(noche), .groups = "drop")


bm_sub <- bm 

# drop Noise and NoID

bm_sub<- bm_sub %>%
  filter(!auto_id %in% c("Noise", "NoID"))

bm_std <- bm_sub %>%
  group_by(site, yr, auto_id) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  left_join(
    bm_sub %>% group_by(site, yr) %>% summarise(nights = n_distinct(noche), .groups = "drop"),
    by = c("site", "yr")
  ) %>%
  mutate(activity = n / nights)

# Convert to long format for codyn
bm_long <- bm_std %>%
  select(site, yr, auto_id, activity)

# Calculate RAC change per site (compare reference year 2021 to others)
rac_change <- RAC_change(
  df = bm_long,
  time.var = "yr",
  species.var = "auto_id",
  abundance.var = "activity",
  replicate.var = "site",
  reference.time = 2021)

# Get treatment info per site
trmt_info <- bm_sub %>% select(site, treatmt) %>% distinct()
rac_change <- left_join(rac_change, trmt_info, by = "site")


# add light

rac_change<-left_join(rac_change, mlight, by= "site")

# add moths

rac_change<- left_join(rac_change, insect, by= "site")

# add elev

rac_change<- left_join(rac_change, elev, by= "site")


# Example model for evenness change
library(lme4)
m1 <- lmer(richness_change ~ treatmt+  yr2+ (1|site), data = rac_change)
summary(m1)


library(ggeffects)
pred <- ggpredict(m1, terms = c("treatmt"))
plot(pred)


m1 <- lmer(richness_change ~ treatmt+ yr2 + + (1|site), data = rac_change)
summary(m1)


library(ggeffects)

pred <- ggpredict(m1, terms = c("treatmt", "yr2"))

ggplot(pred, aes(x = x, y = predicted, color = group)) +
  geom_point(position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2, position = position_dodge(0.3)) +
  labs(x = "Treatment", y = "Change in richness",
       title = "Effect of lighting on bat community change") +
  theme_minimal()


m1 <- lmer(richness_change ~ mean_mwatts + t.leps + yr2 +elev_mean+ (1|site), data = rac_change)
summary(m1)




#updated rank abundance curves

library(codyn)

#sample unit

row_unit<-c("site", "yr")

# 1) Effort (nights) per row unit --------------------------------------
effort_tbl <- bm %>%
  distinct(across(c(all_of(row_unit), noche))) %>%   # unique nights sampled
  count(across(all_of(row_unit)), name = "n_nights")

effort_keep <- effort_tbl %>% filter(n_nights >= 10) # we filter nights site,years that have more than 10 nights 
range(effort_keep$n_nights) # check range of nights sampled per site-year)
# most sites have plenty of data the range of number of night sampled is 28-78


# we use community long table 

comm_long
# add effort
comm_long<-left_join(comm_long, effort_tbl, by= c("site","yr"))
summary(comm_long)

## standardize detection with the effort as detections per day*
comm_long<- comm_long %>% 
  mutate(dpd = n_tot/  n_nights)


# run RCA from cody package 
# 

rac_cambio<-RAC_change(df=comm_long, # this function seems to be the one Jennings 2025 used. 
                       time.var = "yr",
                       species.var = "sp",
                       abundance.var = "dpd", # detections per day
                       replicate.var = "site",
                       reference.time = 2021)

rac_dif<-RAC_difference(df=comm_long,
                        time.var = 'yr',
                        species.var = "sp",
                        abundance.var = "dpd",
                        replicate.var = 'site')


# 4. Join treatments

litsites<-c("iron01","iron03","iron05","long01","long03")

rac_cambio$treatmt<- ifelse(rac_cambio$site %in% litsites, 1, -1)
rac_cambio$tmt<- ifelse(rac_cambio$site %in% litsites, "lit", "dark")

rac_cambio <- rac_cambio %>%
  mutate(yr_s = case_when(
    yr2 == 2021 ~ -1,
    yr2 == 2022 ~ 0,
    TRUE ~ 1
  ))

# add light.
rac_cambio<- rac_cambio %>% left_join(light, by=c("site", "yr2"="yr"))

# add elevation

rac_cambio<- rac_cambio %>% left_join(elev, by=c("site"))


# standardize predictors

rac_cambio<- rac_cambio %>% 
  mutate(mwatts_s= (mwatts - mean(mwatts, na.rm=TRUE))/(2*sd(mwatts, na.rm=TRUE)),
         elev_mean_s= (elev_mean - mean(elev_mean, na.rm=TRUE))/(2*sd(elev_mean, na.rm=TRUE))
  )

# standardize predictors

rac_cambio<- rac_cambio %>% 
  mutate(mwatts_s= (mwatts - mean(mwatts, na.rm=TRUE))/(2*sd(mwatts, na.rm=TRUE)),
         elev_mean_s= (elev_mean - mean(elev_mean, na.rm=TRUE))/(2*sd(elev_mean, na.rm=TRUE))
  )


# 5. Model example (rank change)

library(lme4)
library(lmerTest)

# we ran the models using the treament as light/dark and years as factor (2022,2023) with 2021 as reference (see above).

m_rich <- lmer(richness_change  ~ tmt * as.factor(yr2) + (1|site), data = rac_cambio)
m_even <- lmer(evenness_change ~ tmt * as.factor(yr2) + (1|site), data = rac_cambio)
m_rank <- lmer(rank_change ~ tmt * as.factor(yr2) + (1|site), data = rac_cambio)
m_gain <- lmer(gains ~ tmt * as.factor(yr2) + (1|site), data = rac_cambio)
m_loss <- lmer(losses ~ tmt * as.factor(yr2) + (1|site), data = rac_cambio)

confidenceIntervals <- function(model) {
  confint(model, level = 0.95, method = "Wald")
}


confidenceIntervals(m_rich)


summary(m_rich)
summary(m_even)
summary(m_rank)
summary(m_gain)
summary(m_loss)

plot(m_rich)
plot(m_even)
plot(m_rank)
plot(m_gain)
plot(m_loss)

# 6. Plot marginal effects
library(ggeffects)

pred_rich <- ggpredict(m_rich, terms = c("yr2 [2022,2023]", "tmt")) %>% as.data.frame()
pred_even <- ggpredict(m_even, terms = c("yr2 [2022,2023]", "tmt")) %>% as.data.frame()
pred_rank <- ggpredict(m_rank, terms = c("yr2 [2022,2023]", "tmt")) %>% as.data.frame()
pred_gain <- ggpredict(m_gain, terms = c("yr2 [2022,2023]", "tmt")) %>% as.data.frame()
pred_loss <- ggpredict(m_loss, terms = c("yr2 [2022,2023]", "tmt")) %>% as.data.frame()

# rename 

pred_rich <- pred_rich %>% rename(yr=x,
                                  tmt = group)
pred_even <- pred_even %>% rename(yr=x,
                                  tmt = group)
pred_rank <- pred_rank %>%  rename(yr=x,
                                   tmt=group)
pred_gain <- pred_gain %>% rename(yr=x,
                                  tmt=group)

pred_loss <- pred_loss %>% rename(yr=x,
                                  tmt=group)



rich<-ggplot() +
  # raw site points (semi-transparent)
  geom_point(
    data = rac_cambio,
    aes(x = tmt, y = richness_change, color = as.factor(yr2)),
    position = position_jitter(width = 0.1, height = 0, seed = 1),
    alpha = 0.35, size = 2
  ) +
  # predicted mean and CI
  geom_errorbar(
    data = pred_rich,
    aes(x = tmt, ymin = conf.low, ymax = conf.high, color = as.factor(yr)),
    width = 0.08, linewidth = 1,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    data = pred_rich,
    aes(x = tmt, y = predicted, color = as.factor(yr)),
    size = 3,
    position = position_dodge(width = 0.4)
  ) +
  scale_color_manual(values = c("2022" = "gray50", "2023" = "#1b9e77"), name = "Year") +
  labs(
    x = "", y = "Change in Richness",
    title = "RAC: change in richness by light treatment 2021 as reference"
  ) +
  theme_minimal(base_size = 14)


rich


even<-ggplot() +
  # raw site points (semi-transparent)
  geom_point(
    data = rac_cambio,
    aes(x = tmt, y = evenness_change, color = as.factor(yr2)),
    position = position_jitter(width = 0.1, height = 0, seed = 1),
    alpha = 0.35, size = 2
  ) +
  # predicted mean and CI
  geom_errorbar(
    data = pred_even,
    aes(x = tmt, ymin = conf.low, ymax = conf.high, color = as.factor(yr)),
    width = 0.08, linewidth = 1,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    data = pred_even,
    aes(x = tmt, y = predicted, color = as.factor(yr)),
    size = 3,
    position = position_dodge(width = 0.4)
  ) +
  scale_color_manual(values = c("2022" = "gray50", "2023" = "#1b9e77"), name = "Year") +
  labs(
    x = "", y = "Change in Evenness",
    title = "RAC: change in Evenness by light treatment 2021 as reference"
  ) +
  theme_minimal(base_size = 14)


even




rank<-ggplot() +
  # raw site points (semi-transparent)
  geom_point(
    data = rac_cambio,
    aes(x = tmt, y = rank_change, color = as.factor(yr2)),
    position = position_jitter(width = 0.1, height = 0, seed = 1),
    alpha = 0.35, size = 2
  ) +
  # predicted mean and CI
  geom_errorbar(
    data = pred_rank,
    aes(x = tmt, ymin = conf.low, ymax = conf.high, color = as.factor(yr)),
    width = 0.08, linewidth = 1,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    data = pred_rank,
    aes(x = tmt, y = predicted, color = as.factor(yr)),
    size = 3,
    position = position_dodge(width = 0.4)
  ) +
  scale_color_manual(values = c("2022" = "gray50", "2023" = "#1b9e77"), name = "Year") +
  labs(
    x = "", y = "Change in Rank",
    title = "RAC: change in Rank by light treatment 2021 as reference"
  ) +
  theme_minimal(base_size = 14)


rank


gains<-ggplot() +
  # raw site points (semi-transparent)
  geom_point(
    data = rac_cambio,
    aes(x = tmt, y = gains, color = as.factor(yr2)),
    position = position_jitter(width = 0.1, height = 0, seed = 1),
    alpha = 0.35, size = 2
  ) +
  # predicted mean and CI
  geom_errorbar(
    data = pred_gain,
    aes(x = tmt, ymin = conf.low, ymax = conf.high, color = as.factor(yr)),
    width = 0.08, linewidth = 1,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    data = pred_gain,
    aes(x = tmt, y = predicted, color = as.factor(yr)),
    size = 3,
    position = position_dodge(width = 0.4)
  ) +
  scale_color_manual(values = c("2022" = "gray50", "2023" = "#1b9e77"), name = "Year") +
  labs(
    x = "", y = "Gains",
    title = "RAC: change in Gains by light treatment 2021 as reference"
  ) +
  theme_minimal(base_size = 14)


gains

losses<-ggplot() +
  # raw site points (semi-transparent)
  geom_point(
    data = rac_cambio,
    aes(x = tmt, y = losses, color = as.factor(yr2)),
    position = position_jitter(width = 0.1, height = 0, seed = 1),
    alpha = 0.35, size = 2
  ) +
  # predicted mean and CI
  geom_errorbar(
    data = pred_loss,
    aes(x = tmt, ymin = conf.low, ymax = conf.high, color = as.factor(yr)),
    width = 0.08, linewidth = 1,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    data = pred_loss,
    aes(x = tmt, y = predicted, color = as.factor(yr)),
    size = 3,
    position = position_dodge(width = 0.4)
  ) +
  scale_color_manual(values = c("2022" = "gray50", "2023" = "#1b9e77"), name = "Year") +
  labs(
    x = "", y = "Losses",
    title = "RAC: change in Losses by light treatment 2021 as reference"
  ) +
  theme_minimal(base_size = 14)

losses

RAC<-rich + gains + losses 

# Now I want to make the graph into a three panels with richness, gains, and losses. 
# below there is the code to polish the graph above according to chatgpt

# --- 1) Make sure "Year" is a factor in both raw + predicted data ----------
rac_cambio <- rac_cambio %>%
  mutate(yr2_f = factor(yr2))

pred_rich <- pred_rich %>% mutate(yr_f = factor(yr))
pred_gain <- pred_gain %>% mutate(yr_f = factor(yr))
pred_loss <- pred_loss %>% mutate(yr_f = factor(yr))

# --- 2) Shared palette + shared theme --------------------------------------
year_cols <- c("2022" = "gray50", "2023" = "#1b9e77")

theme_serif <- theme_minimal(base_size = 14, base_family = "serif") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    plot.title   = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

# Helper: same scales in every panel (required for guide collection)
scale_year <- list(
  scale_color_manual(values = year_cols, name = "Year"),
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
)

# --- 3) Build plots (note: keep legends ON; patchwork will collect them) ----
rich <- ggplot() +
  geom_point(
    data = rac_cambio,
    aes(x = tmt, y = richness_change, color = yr2_f),
    position = position_jitter(width = 0.1, height = 0, seed = 1),
    alpha = 0.35, size = 2
  ) +
  geom_errorbar(
    data = pred_rich,
    aes(x = tmt, ymin = conf.low, ymax = conf.high, color = yr_f),
    width = 0.08, linewidth = 1,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    data = pred_rich,
    aes(x = tmt, y = predicted, color = yr_f),
    size = 3,
    position = position_dodge(width = 0.4)
  ) +
  scale_year +
  labs(x = NULL, y = "Change in Richness", title = "RAC: Richness") +
  theme_serif

gains <- ggplot() +
  geom_point(
    data = rac_cambio,
    aes(x = tmt, y = gains, color = yr2_f),
    position = position_jitter(width = 0.1, height = 0, seed = 1),
    alpha = 0.35, size = 2
  ) +
  geom_errorbar(
    data = pred_gain,
    aes(x = tmt, ymin = conf.low, ymax = conf.high, color = yr_f),
    width = 0.08, linewidth = 1,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    data = pred_gain,
    aes(x = tmt, y = predicted, color = yr_f),
    size = 3,
    position = position_dodge(width = 0.4)
  ) +
  scale_year +
  labs(x = NULL, y = "Gains", title = "RAC: Gains") +
  theme_serif

losses <- ggplot() +
  geom_point(
    data = rac_cambio,
    aes(x = tmt, y = losses, color = yr2_f),
    position = position_jitter(width = 0.1, height = 0, seed = 1),
    alpha = 0.35, size = 2
  ) +
  geom_errorbar(
    data = pred_loss,
    aes(x = tmt, ymin = conf.low, ymax = conf.high, color = yr_f),
    width = 0.08, linewidth = 1,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    data = pred_loss,
    aes(x = tmt, y = predicted, color = yr_f),
    size = 3,
    position = position_dodge(width = 0.4)
  ) +
  scale_year +
  labs(x = NULL, y = "Losses", title = "RAC: Losses") +
  theme_serif

# --- 4) Combine with ONE legend --------------------------------------------
RAC <- (rich + gains + losses) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom",
        panel.grid = element_blank()   # removes all grid lines
  )  # applies after collecting

RAC

# trash -------------------------------------------------------------------


# Example for 2021 and 2022
plot_data <- bm_long %>%
  left_join(trmt_info, by = "site") 

ggplot(plot_data, aes(x = reorder(auto_id, -activity), y = activity, fill = treatmt)) +
  geom_bar(stat = "identity") +
  facet_grid(yr ~ treatmt) +
  labs(x = "Species", y = "Standardized Activity per Night") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

