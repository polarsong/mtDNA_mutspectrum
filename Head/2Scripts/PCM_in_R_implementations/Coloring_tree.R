rm(list = ls(all=TRUE))
library(phytools)
library(nlme)
library(geiger)
library(ggtree)
peng_data_tree = read.tree('peng_and_fly_all_data.tre')
birds_tree = read.tree('anc_kg.treefile')
ggtree(birds_tree)

df_flight = read.csv('flight_and_gene.csv')
df_flight = df_flight[,c(2,3,4,5,6)]
df_flight$flight_num = 0
df_flight[df_flight$ability_to_fly != 'Flying',]$flight_num = 1
rownames(df_flight) = df_flight$species_name
flying_tree = read.tree('flying_birds_tree.tre')
name.check(flying_tree, df_flight)
df_flight$ability_to_fly = as.character(df_flight$ability_to_fly)
df_flight$GhAhSkew = as.numeric(as.character(df_flight$GhAhSkew))
df_flight$fAn = as.numeric(as.character(df_flight$fAn))
df_flight$fGn = as.numeric(as.character(df_flight$fGn))
df_flight$flight_num = as.numeric(as.character(df_flight$flight_num))

df_mutspec = read.table('C:/Users/User/Desktop/mutspec12.tsv', header = TRUE, fill = TRUE)
df_ff = df_mutspec[df_mutspec$Label == 'ff',]
df_syn = df_mutspec[df_mutspec$Label == 'syn',]
df_ff = df_ff[!grepl('Node', df_ff$AltNode),]
df_ff = df_ff[,c(1,5,7)]
df_ff1 = reshape(data = df_ff, idvar = 'AltNode',
                 timevar = 'Mut',
                 direction = 'wide')
names(df_ff1) = c('species_name', 'Mutation_TG', 'Mutation_TC', 'Mutation_TA', 'Mutation_GT', 'Mutation_GC', 
                  'Mutation_GA', 'Mutation_CT', 'Mutation_CG', 'Mutation_CA', 'Mutation_AT', 'Mutation_AG', 'Mutation_AC')
df_syn = df_syn[!grepl('Node', df_syn$AltNode),]
df_syn = df_syn[,c(1,5,7)]
df_syn1 = reshape(data = df_syn, idvar = 'AltNode',
                  timevar = 'Mut',
                  direction = 'wide')
names(df_syn1) = c('species_name', 'Mutation_TG_syn', 'Mutation_TC_syn', 'Mutation_TA_syn', 'Mutation_GT_syn', 'Mutation_GC_syn', 
                   'Mutation_GA_syn', 'Mutation_CT_syn', 'Mutation_CG_syn', 'Mutation_CA_syn', 'Mutation_AT_syn', 'Mutation_AG_syn', 'Mutation_AC_syn')
df_flight1 = merge(df_flight, df_ff1)
df_flight1 = merge(df_flight1, df_syn1)

need_species = setdiff(df_flight$species_name, df_flight1$species_name)
correct_need_species = setdiff(df_flight$species_name, need_species)
birds_ms_tree<-keep.tip(flying_tree,correct_need_species)


df_flight1$Mutation_AG = as.numeric(df_flight1$Mutation_AG)
df_flight1$Mutation_AG_syn = as.numeric(df_flight1$Mutation_AG_syn)
df_flight1$Mutation_CT = as.numeric(df_flight1$Mutation_CT)
df_flight1$Mutation_CT_syn = as.numeric(df_flight1$Mutation_CT_syn)



#Starting coloring tree
p <- ggtree(birds_ms_tree, color = df_for_color$Mutation_AG_syn) %<+% df_flight1 + 
  geom_tippoint(aes(color=ability_to_fly)) 
p

df_syn = df_mutspec[df_mutspec$Label == 'syn',]
df_syn = df_syn[,c(1,5,7)]
df_syn = df_syn[df_syn$Mut == 'T>C',]
df_syn = df_syn[,c(2,3)]
df_syn_tips = df_syn[!grepl('Node', df_syn$AltNode),]
df_syn_nodes = df_syn[grepl('Node', df_syn$AltNode),]
#install.packages("phylobase")
library(phylobase)
ab = phylo4d(birds_ms_tree)

names(df_syn) = c('Mut', 'MutSpec', 'species_name')
df_for_color = merge(df_flight1, df_syn)




tip <- get.tree(birds_ms_tree)$tip.label

beast_tree <- groupOTU(birds_ms_tree, tip[grep("Anas", tip)], 
                       group_name = "host")
