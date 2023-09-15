#birds phyloPCA
rm(list = ls(all=TRUE))
library(ape)
library(phangorn)
library(phytools)
library(geiger)
#Reading and plotting tree
birds_tree<-read.tree(file="anc_kg.treefile")
birds_tree
plotTree(birds_tree,ftype="i", type = 'fan', fsize=0.4,lwd=1)
#Preparing data
df_mtdna = read.csv('Birds_dataset_paper.csv')
df_eco = read.csv('flying_birds.csv')
df_names = as.data.frame(unique(df_mtdna$species_name))
names(df_names) = c('species_name')
df_eco = df_eco[,c(2,3,4)]
names(df_eco) = c('species_name', 'Flight','Diving')
df_eco = merge(df_names, df_eco)
df_eco$species_name = gsub(' ', '_', df_eco$species_name)
name.check(birds_tree,df_eco$species_name)
