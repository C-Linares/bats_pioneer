# ---------------------------
##
## Script name:  bat_diversity.R
##
## Purpose of script: asses bat diversity
##
## Author: Carlos Linares
##
## Date Created: 05/21/2024
##
## Email: carlosgarcialina@u.boisestate.edu
##
## inputs
##
##
##
## ---------------------------
##
## Notes: started in 2022 and reworked in 2024 
##   
##
## ---------------------------
## 

# load libraries 
library(vegan)
library(ggplot2)
library(iNEXT)
library(tidyverse)
library(lubridate)
library(viridis)


# load data ---------------------------------------------------------------

# load kpro data. 

bm<-read_csv('data_for_analysis/prep_for_glmm/bm.csv')


# bat2021$datetime<-ymd_hms(bat2021$datetime) # makes dates as dates 
# bat2021$noche<-ymd(bat2021$noche)# makes dates as dates

# unique(bat2021$site)# tell us what sites we have



# Matrix building ------------------

# Summarize data to get total counts per species per site
abundance_data <- bm %>%
  group_by(site, AUTO.ID.) %>%
  summarize(count = sum(n), .groups = 'drop') %>%
  spread(key = site, value = count, fill = 0)


# Convert to matrix and set species names as row names
abundance_matrix <- as.matrix(abundance_data[,-1])
rownames(abundance_matrix) <- abundance_data$AUTO.ID.

# Check the matrix
head(abundance_matrix)
summary(abundance_data) # there is no NAs


#Sum the abundances of each species across all sites
total_abundance <- rowSums(abundance_matrix)

species_counts <- data.frame(
  Species = names(total_abundance),
  Count = as.numeric(total_abundance)
)
species_counts<- species_counts %>% filter(!Species %in% c("Noise", "NoID"))

sum(species_counts$Count)

summary(species_counts)
# Create and style the table
kable(species_counts, format = "html", caption = "Species Count Table") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

p1<-ggplot(species_counts, aes(x = reorder(Species, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Species Call Counts 2021-2023", x = "Species", y = "Count") +
  annotate("text", x = Inf, y = Inf, label = "Total = 560,645", hjust = 2, vjust = 1, size = 5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_blank(base_size = 12, base_family = "")
p1


p1 <- ggplot(species_counts, aes(x = reorder(Species, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "white") +  # Change bar fill color to black
  labs(title = "Species Call Counts 2021-2023", x = "Species", y = "Count") +
  annotate("text", x = Inf, y = Inf, label = "Total = 560,645", hjust = 2, vjust = 1, size = 5, color = "white") +  # Set text color to white
  theme_minimal(base_size = 16) + 
  theme(
    plot.background = element_rect(fill = "black"),  # Set plot background to black
    panel.background = element_rect(fill = "black"),  # Set panel background to black
    axis.text = element_text(color = "white"),  # Set axis text color to white
    axis.title = element_text(color = "white"),  # Set axis title color to white
    plot.title = element_text(color = "white"),  # Set plot title color to white
    axis.ticks = element_line(color = "white"),  # Set axis ticks color to white
    panel.grid.major = element_line(color = "black"),  # Optional: set grid lines color (adjust as needed)
    panel.grid.minor = element_line(color = "black"),# Optional: set minor grid lines color (adjust as needed)
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text by 45 degrees
  )

p1

ggsave("speciescounts_v2.tiff", plot = p1, device = "tiff",path = 'figures/bat_diversity', width =10 , height =5,units = "in" )

# Identify the species with the highest total abundance
most_abundant_species <- names(which.max(total_abundance))

# Print the most abundant species
print(paste("The most abundant species is:", most_abundant_species))

# Perform iNEXT analysis
iNEXT_result <- iNEXT(abundance_matrix, q=c(0,1,2),datatype = "abundance")

# Plot results
ggiNEXT(iNEXT_result)

# Manually specify shapes and colors
shapes <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)  # Up to 15 different shapes
colors <- rainbow(15)  # Up to 15 different colors

# Customize the ggiNEXT plot
ggiNEXT(iNEXT_result,type = 1) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  theme_minimal()

ggiNEXT(iNEXT_result, type = 3) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  theme_minimal()


# heat map
abundance_data<-filter(abundance_data, AUTO.ID.!=c("Noise","NoID")) #

long_data <- abundance_data %>%
  pivot_longer(cols = -AUTO.ID., names_to = "Site", values_to = "Abundance") %>%
  mutate(Site = factor(Site, levels = unique(Site))) %>%
  mutate(Species = factor(AUTO.ID.))

# Plot the heatmap
ggplot(long_data, aes(x = Site, y = Species, fill = Abundance)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "Site number", y = "Species", fill = "No. detections") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#Group by site and treatment, then summarize the total counts for each species

bm<-filter(bm, !AUTO.ID. %in% c("Noise","NoID")) #
species_abundance <- bm %>%
  group_by(site, treatmt, AUTO.ID.) %>%
  summarize(total_count = sum(n), .groups = 'drop') %>%
  arrange(site, treatmt, desc(total_count))

# Get the most abundant species for each site and treatment
most_abundant_species <- species_abundance %>%
  group_by(site, treatmt) %>%
  slice_max(total_count, n = 1, with_ties = FALSE)

print(most_abundant_species)


 





####################################################### plots #########################


# Alpha Diversity -------------------------------------------------------

bat.sp.m<- bm %>% 
  arrange(site) %>% 
  group_by(site) %>% 
  count(AUTO.ID.) %>%  #counts the occurrences not the number of calls
  spread(AUTO.ID., n, fill=0)

bat.sp.m<-bat.sp.m %>% select(-c(NoID,Noise))

tab<-t( bat.sp.m[,-1])
totals2021<- rowSums(tab) # bat passes totals for 2021
totals2021
write.table(totals2021)
# Alpha

sps<-specnumber(bat.sp.m)
shanon<- diversity(bat.sp.m[-1])
ens<- exp(shanon)
totbats <-
  rowSums(bat.sp.m != 0)
sites<-unique(bat.sp.m$site)
alpha.div <-cbind.data.frame(sites = bat.sp.m[1] , sps, totbats, shanon, ens) 
alpha.div

litsites<-c("iron01","iron03","iron05","long01","long03")


alpha.div$treatmt<-ifelse(alpha.div$site %in% litsites , "lit", "dark")


#----- plot alpha ------
  
p1.1 <- ggplot(data=alpha.div, aes(x=site, y=ens, color=treatmt)) +
  geom_point(size=4) +
  ylim(10,14)+
  ylab("effective number of species")+
  theme_classic() +
  scale_color_viridis(discrete = T, option = "A",begin = 0, end = .5)
p1.1

ggsave("spnum.tiff", plot = p1.1, device = "tiff",path = 'figures/bat_diversity')



bat.sp.m<-column_to_rownames(bat.sp.m, var = "site") #makes sites row names 

h<-diversity(bat.sp.m)
simp<-diversity(bat.sp.m, "simpson")
invsimp<-diversity(bat.sp.m, "inv")
unb.simp<- simpson.unb(bat.sp.m)
# fisher<-fisher.alpha(bat.sp.m)# not working 

pairs(cbind(h, simp, invsimp, unb.simp), pch="+", col="blue") # don't know how to interpret...

# lets make a species curve. 

curve.sp<-specaccum(bat.sp.m,method = "random")

plot(curve.sp, ci.type = c("bar"), ylab = curve.sp$method)
boxplot(curve.sp)



# simple  plot of sites by bat call abundance. (not here I am using all the data we need to filter out the UV bucket days)

p1<- ggplot(bat.m, aes(x=site,y=n, fill=treatmt))+
  geom_col()+
  theme_minimal()

p1


bat.site.m<-column_to_rownames(bat.site.m, var = "SppAccp")
# b <- as.data.frame(bat.site.m[, -1]) previous line does the same and might work better

out1<-iNEXT(bat.site.m, q=c(0,1,2), datatype = "abundance", endpoint = 100000)

q0<- fortify(out1, type = 1)
q0.point<- q0[which(q0$Method=="Observed"),]

ggplot(q0, aes(x=x, y=y, fill=Assemblage)) +
  geom_point(aes(shape=Assemblage), size=5, data=q0.point)

ggiNEXT(out1, type = 1, facet.var = "None", color.var = "Assemblage")

p1<-ggiNEXT(out1, type = 1)+
  facet_wrap(, scales = "free")+
  theme_classic()
p1

p1<- ggiNEXT(out1,type = 1,)+
  theme_classic()
p1

plot(out1, type = 1)
plot(out1, type = 2)
plot(out1, type = 3)

ggiNEXT(out1, se=T,facet.var = "None", type = 1, color.var = "Assamblage", grey = F)


beta.bat<-betadiver(bat.site.m, "w")
range(beta.bat - vegdist(bat.site.m, binary = T))

# -----------------Beta -----------------------

metadata <- data.frame(
  site = bat.sp.m$site,
  Group = c("lit","d","lit","d","lit","d","lit","d","lit","d","d","d","d","d","d")) # Example grouping


# Calculate Bray-Curtis dissimilarity
beta_div <- vegdist(bat.sp.m[,-1], method = "bray")

# Ensure the site order in metadata matches the site order in bat.sp.m
metadata <- metadata[match(bat.sp.m$site, metadata$site), ]

# Perform PERMANOVA
adonis_result <- adonis2(beta_div ~ Group, data = metadata)

# -----------------------betapart-----------------

library(betapart)
library(dendextend)
library(ggdendro)
binary_bat_sp_m <- bat.sp.m %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))

beta.div<-beta.pair(binary_bat_sp_m[,-1], index.family = "sorensen")
print(beta.div)

# Perform hierarchical clustering
hc <- hclust(beta.div$beta.sor, method = "average")
# Convert to dendrogram object
dend <- as.dendrogram(hc)
#Convert to data frame for ggplot
dend_data <- ggdendro::dendro_data(dend)

# Plot using ggplot2
ggplot() +
  geom_segment(data = dend_data$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dend_data$labels, aes(x = x, y = y, label = label), hjust = 1, size = 3) +
  labs(title = "Dendrogram of Beta Diversity 2021") +
  ylab("Sorensen disimilarity")+
  theme(plot.title = element_text(hjust = 0.5))

plot(hclust(beta.div$beta.sim, method="average"), hang=-1)
plot(hclust(beta.div$beta.sor, method="average"))
plot(density(beta.div$beta.sor), xlim=c(0,0.8), ylim=c(0, 19), xlab='Beta diversity', main="", lwd=3)

# adespatial

# install.packages("adespatial")
library(adespatial)

# BAT package

# install.packages("BAT")
library(BAT)

beta.div<-beta(bat.sp.m[,-1], abund = T) # jaccard beta diversity. 
print(beta.div)

# Convert beta_div to a distance matrix
beta_dist <- as.dist(beta.div$Brich)

# Perform PCoA
pcoa_result <- cmdscale(beta_dist, eig = TRUE, k = 2)  # k is the number of dimensions
pcoa_result$eig
# Extract coordinates
coordinates <- as.data.frame(pcoa_result$points)

# Plot the PCoA results
plot(coordinates$V1, coordinates$V2, xlab = "PCoA1", ylab = "PCoA2", main = "PCoA of Beta Diversity")
text(coordinates$V1, coordinates$V2, labels = rownames(coordinates), pos = 3)


coordinates_df <- as.data.frame(pcoa_result$points)
coordinates_df$Site <- rownames(coordinates_df)

# Create a ggplot
ggplot(coordinates_df, aes(x = V1, y = V2, label = Site)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 0.5) +
  labs(x = "PCoA1", y = "PCoA2", title = "PCoA of Beta Diversity") +
  theme_minimal()

ggplot(coordinates_df, aes(x = V1, y = V2)) +
  geom_point(size = 3) +
  geom_text(aes(label = c("iron01", "iron02", "iron03", "iron04", "iron05", 
                          "iron06", "long01", "long02", "long03", "long04", 
                          "long05", "viz01", "viz02", "viz03", "viz04")),
            vjust = -1, hjust = 0.5) +  # Adjust vjust and hjust for label position
  labs(title = "PCoA of Beta Diversity", x = "PCoA1", y = "PCoA2") +
  theme_minimal()

# beta1<-beta.multi(bat.sp.m[,-1], func = "jaccard", abund = TRUE, raref = 2, runs = 100) not working
# con<-contribution(bat.sp.m[,-1], abund = TRUE) Not wortking 
#---------------------------------


bat.sp.m<-column_to_rownames(bat.sp.m, var = "site")
bat.sp.m<-bat.sp.m[-12,]
bat.sp.m<-bat.sp.m[,-4] # reomove sp EDU
# lets try without vizc01

bats.bc.nmds <- metaMDS(bat.sp.m, k=7, trymax=500,trace = T) 
ordiplot(bats.bc.nmds, type = "t",display = "sites", cex = 0.6)
points(bats.bc.nmds, "sites", pch=21, col="blue", bg="yellow")
scores(bats.bc.nmds)

p6<-ordiplot(bats.bc.nmds, type = "none")
points(p6, "sites", pch=21, col="red", bg="yellow")
text(p6, "sites", col="blue", cex=0.9)









# ------------------------------Figure bat response by site lit vs dark ----------------------

p2<- ggplot(bat.m, aes(x=hr,y=n))+
  geom_col()+
  facet_grid(.~ treatmt)+
  theme_classic()

p2

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p2<- ggplot(data = bat.m, aes(factor(hr, levels = c(20,21,22,23,0,1,2,3,4,5,6)), n, fill=treatmt))+
         geom_col()+
  facet_grid(.~treatmt)+
  theme_classic()+
  xlab("hours")+
  ylab("Calls counts 2021")+
  scale_fill_manual(values=cbbPalette) # inclued blind friendly 
  # guides(fill = "none") # removes the legend

p2
ggsave(p2,filename = "p2.png", path = "data_analysis/figures/")
       
p3<- ggplot(data = bat.m, aes(factor(hr, levels = c(20,21,22,23,0,1,2,3,4,5,6)), n, fill=SppAccp))+
  geom_col()+
  facet_grid(.~treatmt)+
  theme_classic()+
  xlab("hours")+
  ylab("Calls counts 2021")+
  scale_fill_manual(values=(hcl.colors(n = 12, palette = "Temps"))) # inclued blind friendly 
# guides(fill = "none") # removes the legend

p3

ggsave(filename = "p2.png", path = "data_analysis/figures/")
ggsave(filename = "p3.png", path = "data_analysis/figures/")




# save --------------------------------------------------------------------

save.image(file = "working_env/bat_divesity_v2.RData")

load(file = 'working_env/bat_divesity_v2.RData')



# junk --------------------------------------------------------------------
# species matrix with species as columns 

bat.site.m<- bat2021 %>% 
  arrange(site) %>% 
  group_by(site) %>% 
  count(SppAccp) %>% 
  spread(site, n, fill=0) #site matrix with sites as columns

bat.m<- bat2021 %>% 
  group_by(site) %>% 
  count(SppAccp,noche,hr) %>% 
  mutate(jday=yday(noche)) %>% 
  ungroup()# for running models 

bat.pred<- select(bat2021, "noche","site", "fraction", 'treatmt','trmt_bin') %>% 
  distinct()


bat.m<- left_join(bat.m, bat.pred)


bat.sp.m<- bat2021 %>% 
  arrange(site) %>% 
  group_by(site) %>% 
  count(SppAccp) %>% 
  spread(SppAccp, n, fill=0
#---------------Activity Index --------------------#
# with these two matrix we try to calculate the miller activity index.


bat.m2<-bat2021 %>% # counts every minute but I make them presence absence later
  group_by(site,SppAccp,noche,hr,min) %>% 
  count(SppAccp) %>% 
  mutate(pre.index= ifelse(n>=1,1,0)) %>% # this makes any counts >1 into 1 
  spread(SppAccp,pre.index,fill = 0)

# bat.m2long<- 

bat.m3<- bat2021 %>% 
  group_by(site,SppAccp,noche,hr,min) %>% 
  count(SppAccp) %>% spread(SppAccp, n, fill=0)

bat.m4<-bat.m3 %>% pivot_longer(-c(site,noche,hr,min,), names_to = "spp", values_to = "ai_miller")# long version of bat.m3 I don't know if I need a long format. I guess that for running the analysis I do. 


recdates<-bat2021 %>% # this gives us the number of days each recorder was on. 
  group_by(site) %>% 
  summarise(startdate= min(datetime), enddate=max(datetime)) %>% 
  mutate(n.days= time_length(enddate - startdate, unit = "days"))



