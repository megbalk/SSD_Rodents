#1) Calculate the SSD index for each species 
#2) Look at Rensch’s Rule across the group – phylogenetic reduced major axis (pRMA) regression
#3) Rates of SSD (index) evolution across the phylogeny
#4) PGLS to examine SSD with life history traits

library(dplyr)
library(ape)
library(geiger)
library(phytools)
library(diversitree)
library(nlme)
library(qpcR)
library(reshape)
library(calibrate)

setwd("/Users/Maggie/Dropbox/Mammal_SSD/Rodents/")

#Body mass
rodent_mass <- read.csv("BM_rodent_data_SSD.csv", header = TRUE, stringsAsFactors = FALSE)
str(rodent_mass)
rodent_mass$scientificName <- as.factor(rodent_mass$scientificName)
plyr::count(rodent_mass$scientificName)

######################################################
#SSD Index

BM.SDI <- rodent_mass %>%
  group_by(scientificName) %>%
  dplyr::summarize(SSD = (-(mean_BM_male/mean_BM_female)+1)*100) 

rodent_mass <- right_join(rodent_mass, BM.SDI, by="scientificName")

#####################################################################
# Read in tree
mam.tree <- read.tree("mammal2.tre")
str(mam.tree)
table(mam.tree$tip.label)
plyr::count(mam.tree$tip.label)

#Body mass phylogeny 
BM_species <- read.csv("BM_species.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
#checking match
matching.files <- name.check(mam.tree, BM_species) 
matching.files$data_not_tree #all match

#re-upload
BM_species <- read.csv("BM_species.csv")
plyr::count(BM_species$scientificName)

pruned.tree.BM <-drop.tip(mam.tree, mam.tree$tip.label[-na.omit(match(BM_species[,1],mam.tree$tip.label))])
plyr::count(pruned.tree.BM$tip.label)
str(pruned.tree.BM)
pruned.tree.BM$tip.label
plot(pruned.tree.BM)
plotTree(pruned.tree.BM,ftype="i",lwd=1,fsize=0.3,type="fan",part=1)

#####################################################################

SSD_frame_BM <- data.frame(rodent_mass[,2:10]) 
rownames(SSD_frame_BM) <- rodent_mass[,1]
SSD_frame_BM

#####################################################################
# function for matching tips and sorting 
treeresults_BM <- treedata(pruned.tree.BM, SSD_frame_BM, sort=TRUE, warnings=TRUE)
treeresults_BM

# analysis tree is the pruned phy
analysis_tree_BM <- treeresults_BM$phy

# grab matched data (sorted)
matched_data_BM <- treeresults_BM$data
matched_data_BM

#####################################################################
#stand-alone variables for regression

log_mass_m <- as.numeric(matched_data_BM[,7])
names(log_mass_m) <- rownames(matched_data_BM)

log_mass_f <- as.numeric(matched_data_BM[,8])
names(log_mass_f) <- rownames(matched_data_BM)

ssd_bm <- as.numeric(matched_data_BM[,9])
names(ssd_bm) <- rownames(matched_data_BM)

summary(log_mass_m)
summary(log_mass_f)
summary(ssd_bm)
#####################################################################

# Rensch - body mass

RMA_BM <- phyl.RMA(log_mass_f, log_mass_m, analysis_tree_BM, method="BM", lambda=NULL, fixed=FALSE)
RMA_BM$test
RMA_BM$RMA.beta

plot(y=log_mass_m, x=log_mass_f, ylab="log mass males", xlab="log mass females", pch=21, bg="grey", cex=2.5)
abline(RMA_BM$RMA.beta, lwd=3)
# plot the isometry line
abline(0,1, lwd=3, lty=3)




