# ---------------------------
##
## Script name:  phylo_reg.R
##
## Purpose
# this script will is to answer the question: Do bats respond to light and is there a phylogenetic signal?
# 
## Author: Carlos Linares, 
## Date Created: 9/11/2025
##
## Email: carlosgarcialina@u.boisestate.edu
##
## ---------------------------
##
## Notes: 
##   
##
## ---------------------------
## # inputs ------------------------------------------------------------------
#' - bat_combined, file = 'data_for_analysis/prep_for_glmm/bat_combined.csv'
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
  "phyr",
  "ape"
)

library(phangorn)  # RF.dist

# data 

# one of the first steps for this analysis is to get a pruned phylogenetic tree with just the species of interest from the following source (https://vertlife.org/phylosubsets/)



# Load data
bats <- read_csv("data_for_analysis/Species_bats.csv") %>%
  clean_names()

bm_species <- read_csv("data_for_analysis/prep_for_glmm/bm.csv") %>%
  clean_names() %>%
  rename(sp = auto_id) %>%
  filter(!sp %in% c("Noise", "NoID")) %>%   # remove unwanted rows
  left_join(bats, by = c("sp" = "six_letter_species_code")) %>% 
  distinct(sp, scientific_name)             # keep one row per species


# export just the scientific names as a txt
# 
write_lines(
  bm_species$scientific_name,
  "data_for_analysis/bat_species.txt"
)

# load vertilife taxonomy because it did not work. let's try to trim it. 

vligfe<- read_csv('data_for_analysis/vertlife_taxonomies.csv') %>% 
  clean_names()


# read tree

arbol<- read.nexus('data_for_analysis/tree/output.nex') # tree from vertilife 
class(arbol)
length(arbol)


# Extract the first tree
single_tree <- arbol[[1]]



plot.phylo(single_tree, show.tip.label = TRUE)

# code below should help us evaluate if the trees are stable. 

tipsets_same <- all(sapply(arbol, function(tr) 
  identical(sort(tr$tip.label), sort(arbol[[1]]$tip.label))
))
if(!tipsets_same) stop("Tip labels differ across trees: you must standardize/prune tip sets first")

# 2) Build a majority-rule consensus tree (splits present in >= 50% of trees)

consensus_tree <- consensus(arbol, p = 0.5)   # p=0.5 is majority-rule
plot(consensus_tree, show.tip.label = TRUE, cex = 0.7)
title("Majority-rule consensus (p = 0.5)")

# 3) Calculate clade (node) support: fraction of trees that contain each split in the consensus
# prop.clades returns counts of trees that contain each clade in 'consensus_tree'
clade_counts <- prop.clades(consensus_tree, arbol)
clade_support <- clade_counts / length(arbol)  # proportion 0..1

# attach support to consensus tree node labels (internal nodes)
consensus_tree$node.label <- round(clade_support, 3)
plot(consensus_tree, show.node.label = TRUE, cex = 0.7)
title("Consensus tree with node support (proportion of trees)")

# Interpretation helper: print nodes with low support
low_support_nodes <- which(clade_support < 0.5)
if(length(low_support_nodes)) {
  message("Nodes with support < 0.5 (low): ", paste(low_support_nodes, collapse = ", "))
} else {
  message("All internal nodes have >= 50% support")
}


# load 