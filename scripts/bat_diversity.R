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
library(janitor)
library(faux)


# load data ---------------------------------------------------------------

# load kpro data. 

t<-read.csv('data_for_analysis/data_glmm/bat_counts.csv', header = T) # bat counts by jday  
bm<- read_csv('data_for_analysis/prep_for_glmm/bm.csv') %>% clean_names()
  

# bat2021<-read.csv('data_analysis/bat2021_v2.csv',check.names = T)
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

# Print the total abundances
print(total_abundance)

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




#Group by site and treatment, then summarize the total counts for each species

bm<-filter(bm, AUTO.ID.!=c("Noise","NoID"))
species_abundance <- bm %>%
  group_by(site, treatmt, AUTO.ID.) %>%
  summarize(total_count = sum(n), .groups = 'drop') %>%
  arrange(site, treatmt, desc(total_count))

# Get the most abundant species for each site and treatment
most_abundant_species <- species_abundance %>%
  group_by(site, treatmt) %>%
  slice_max(total_count, n = 1, with_ties = FALSE)

print(most_abundant_species)


bat.sp.m<- bat2021 %>% 
  arrange(site) %>% 
  group_by(site) %>% 
  count(SppAccp) %>% 
  spread(SppAccp, n, fill=0) # species matrix with species as columns 

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







####################################################### plots #########################


# Alpha Diversity -------------------------------------------------------

tab<-t( bat.sp.m[,-1])
totals<- rowSums(tab) # bat passes totals for 2021
totals
write.table(totals)
# Alpha

sps<-specnumber(bat.sp.m)
shan<- diversity(bat.sp.m[-1])
ens<- exp(shan)
totbats <-
  rowSums(bat.sp.m != 0)

alpha <-cbind.data.frame(sites = bat.sp.m[1] , sps, totbats, shan, ens) 
alpha

#----- plot alpha ------
  
p1.1 <- ggplot(data=alpha, aes(x=site, y=sps,size=5)) +
  geom_point(colour="brown") +
  # ylim(0,10) +
  ylab("bat species richness")+
  theme_classic() 
p1.1


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
# mod1<- fitspecaccum(curve.sp, "lomolino")# error
# what doe this mean? there is great variation between sites, but as we add more sites we seemed to have capture all the specie present in the area. 



# simple  plot of sites by bat call abundance. (not here I am using all the data we need to filter out the UV bucket days)

p1<- ggplot(bat.m, aes(x=site,y=n, fill=treatmt))+
  geom_col()+
  theme_minimal()

p1

# need to relearn the I next package. 


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


bat.sp.m<-column_to_rownames(bat.sp.m, var = "site")
bat.sp.m<-bat.sp.m[-11,]

bats.bc.nmds <- metaMDS(bat.sp.m, k=2, trymax=100,trace = T) 
ordiplot(bats.bc.nmds, type = "t",display = "sites", cex = 0.6)
points(p5, "sites", pch=21, col="blue", bg="yellow")


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



# -----------------------modeling ---------------------
library(lme4)
library(sjPlot)
library(brms)
library(glmmTMB)
library(car)
library(ggeffects)
library(emmeans)

bat.m$site<-(bat.m$site)
bat.m$jday<-yday(bat.m$noche)
hist(bat.m$jday)# I don't expect jday to have much of an effect
hist(bat.m$n )

m1 <- glmer( n ~ jday+ hr + treatmt + (1|site), data = bat.m, family = poisson(link = "log"))
summary(m1)
plot_model(m1)
lattice::dotplot(ranef(m1, postVar = TRUE))

m1.zi<- glmmTMB(n ~ jday+hr+treatmt+ (1|site), data = bat.m, ziformula = ~1)
summary(m1.zi)

# is the data zero inflated?
ggplot(bat.m, aes(n)) + geom_histogram() + scale_x_log10() # I don't think the data is zero inflated

m2<- glmer(n~ jday +hr +treatmt + phase + (1|site), data = bat.m, family = poisson(link = "log"))

summary(m2)
plot_model(m2)
plot_model(m2,type = "re")
plot_model(m2, type = "pred")


AIC(m1,m2)
anova(m1,m2)

m3<- glmer(n ~ jday+hr+treatmt+phase+fraction+(1|site), data = bat.m, family = poisson(link = "log"))

summary(m3)
plot_model(m3)
plot_model(m3,type = "re")
plot_model(m3, type = "pred")


AIC(m1,m2,m3)
anova(m1,m2,m3)

# with scaled variables phase, fraction, and jday.
m4 <-
  glmer(
    n ~ scale(jday) + scale(hr) + trmt_bin + scale(phase) + scale(fraction) + scale(dlt.sunset)+
      (1 | site),
    data = bat.m,
    family = poisson(link = "log")
  )
summary(m4)
vif(m4)
# which variables should be scaled? 
# apparelty we should scale prdictor varialbles that have different ranges. 

plot_model(m4)
plot_model(m4, type = "std")
plot_model(m4, type = "pred", terms = "treatmt")

m1<- 
  glmer(n~ treatmt + (1|site),
        data = bat.m,
        family = poisson)
plot_model(m1)



m1.2<- glmer(n ~ scale(jday) + scale(fraction) + SppAccp + treatmt+ (1|site), data = bat.m, family = poisson)
summary(m1.2)
vif(m1.2)
plot_model(m1.2,transform = "exp", sort.est = T, title = "")

m1.3<- glmer(n ~ scale(jday) + scale(fraction)  + trmt_bin + (1|site), data = bat.m, family = poisson)
summary(m1.3)
vif(m1.3)
plot_model(m1.3,transform = "exp", sort.est = T, title = "")
plot_model(m1.3, "re")

e<-emmeans(m1.3, "trmt_bin", type= "response")
plot(e)

pr<- ggpredict(m1.3, "trmt_bin", type = "count")

m1.3<- glmmTMB(n ~ treatmt + (1|site), data = bat.m, family = poisson(link = "log"))
pr <- ggpredict(m1.3, terms = "treatmt")
plot(m1.3)
# modeling individual bat response ----------------------------------------

spss<- unique(bat.m$SppAccp)
epfu<- bat.m[bat.m$SppAccp == "Epfu", ]
laci<- bat.m[bat.m$SppAccp == "Laci", ]
lano<- bat.m[bat.m$SppAccp== "Lano",]
myca<- bat.m[bat.m$SppAccp== "Myca",]
myci<- bat.m[bat.m$SppAccp== "Myci",]
mylu<- bat.m[bat.m$SppAccp== "Mylu",]
myvo<- bat.m[bat.m$SppAccp== "Myvo",]
myev<- bat.m[bat.m$SppAccp== "Myev",]

m.epfu<- glmer(n ~ scale(jday) + scale(fraction) + trmt_bin + (1|site), data = epfu, family = poisson)
summary(m.epfu)
vif(m.epfu)

plot_model(m.epfu,transform = "exp", sort.est = T, title = "Eptesicus fuscus")
plot_model(m.epfu, type = "re")

m.laci<- glmer(n~ scale(jday)+ trmt_bin + scale(fraction) + (1|site), data = laci, family = poisson)
summary(m.laci)
plot_model(m.laci, title = "Lasiurus cinereus")

m.lano<- glmer(n ~ scale(jday) + scale(fraction) + trmt_bin + (1|site), data = lano, family = poisson)
summary(m.lano)
plot_model(m.lano, title= "Lasionycteris noctivagans")

m.myca<- glmer(n ~ scale(jday) + scale(fraction) + trmt_bin + (1|site), data = myca, family = poisson)
summary(m.myca)
plot_model(m.myca, title = "Myotis californicus")

m.myci<- glmer(n ~ scale(jday) + scale(fraction) + trmt_bin + (1|site), data = myci, family = poisson)
summary (m.myci)
plot_model(m.myca, title = "Myotis ciliolabrum ")

m.myvo<- glmer(n ~ scale(jday) + scale(fraction) + trmt_bin + (1|site), data = myvo, family = poisson)
summary (m.myvo)
plot_model(m.myvo, title = "Myotis volans")

m.myev<- glmer(n ~ scale(jday) + scale(fraction) + trmt_bin + (1|site), data = myev, family = poisson)
summary (m.myev)
plot_model(m.myev, title = "Myotis evotis")





# junk --------------------------------------------------------------------


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



