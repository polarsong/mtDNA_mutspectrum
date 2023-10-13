rm(list = ls(all=TRUE))
library(phytools)
library(nlme)
library(geiger)
library(ggtree)
peng_data_tree = read.tree('peng_and_fly_all_data.tre')
birds_tree = read.tree('anc_kg.treefile')
df_mtdna = read.csv('Birds_dataset_paper.csv')

df_need = data.frame()
for (i in unique(df_mtdna$species_name))
{
  a = df_mtdna[df_mtdna$species_name == i,]
  b = sum(a$chthSkew)/12
  ab = c(i,b)
  df_need = rbind(df_need, ab)
}

names(df_need) = c('species_name', 'ChThSkew')
df_need$species_name = gsub(' ', '_', df_need$species_name)

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

df_flight1$GhAhSkew = as.numeric(df_flight1$GhAhSkew)
df_flight1$Mutation_AG = as.numeric(df_flight1$Mutation_AG)
df_flight1$Mutation_AG_syn = as.numeric(df_flight1$Mutation_AG_syn)
df_flight1$Mutation_CT = as.numeric(df_flight1$Mutation_CT)
df_flight1$Mutation_CT_syn = as.numeric(df_flight1$Mutation_CT_syn)
rownames(df_flight1)= df_flight1$species_name
df_flight1$double_fly = 0
df_flight1[df_flight1$ability_to_fly == 'Flying',]$double_fly = 1
df_temp_fly = read.csv('birds_metrics.csv')
df_temp_fly = df_temp_fly[,c(2:17)]
df_flight2 = merge(df_flight1, df_temp_fly)
need_species = setdiff(df_flight1$species_name, df_flight2$species_name)
correct_need_species = setdiff(df_flight1$species_name, need_species)
birds_ms_and_temp_tree<-keep.tip(birds_ms_tree,correct_need_species)
row.names(df_flight2) = df_flight2$species_name
df_flight2$Latitude = gsub(',', '.', df_flight2$Latitude)
df_flight2$TempRange = gsub(',', '.', df_flight2$TempRange)
df_flight2$AnnualPrecip = gsub(',', '.', df_flight2$AnnualPrecip)
df_flight2$PrecipRange= gsub(',', '.', df_flight2$PrecipRange)
df_flight2$AnnualTemp= gsub(',', '.', df_flight2$AnnualTemp)

df_flight2$AnnualTemp = as.numeric(as.character(df_flight2$AnnualTemp))
df_flight2$Beak_length_Culmen = as.numeric(as.character(df_flight2$Beak_length_Culmen))
df_flight2$Beak_length_Nares = as.numeric(as.character(df_flight2$Beak_length_Nares))
df_flight2$Beak_width = as.numeric(as.character(df_flight2$Beak_width))
df_flight2$Beak_depth = as.numeric(as.character(df_flight2$Beak_depth))
df_flight2$Tarsus_length = as.numeric(as.character(df_flight2$Tarsus_length))
df_flight2$Wing_length = as.numeric(as.character(df_flight2$Wing_length))
df_flight2$Kipps_distance = as.numeric(as.character(df_flight2$Kipps_distance))
df_flight2$Hand_wing_index = as.numeric(as.character(df_flight2$Hand_wing_index))
df_flight2$Tail_length = as.numeric(as.character(df_flight2$Tail_length))
df_flight2$Mass = as.numeric(as.character(df_flight2$Mass))
df_flight2$Latitude = as.numeric(as.character(df_flight2$Latitude))
df_flight2$AnnualTemp = as.numeric(as.character(df_flight2$AnnualTemp))
df_flight2$TempRange = as.numeric(as.character(df_flight2$TempRange))
df_flight2$AnnualPrecip = as.numeric(as.character(df_flight2$AnnualPrecip))
df_flight2$PrecipRange = as.numeric(as.character(df_flight2$PrecipRange))

df_flight2$Mass = log10(df_flight2$Mass)
df_flight2$Beak_length_Culmen = log10(df_flight2$Beak_length_Culmen)
df_flight2$Beak_length_Nares = log10(df_flight2$Beak_length_Nares)
df_flight2$Beak_width = log10(df_flight2$Beak_width)
df_flight2$Beak_depth = log10(df_flight2$Beak_depth)
df_flight2$Tarsus_length = log10(df_flight2$Tarsus_length)
df_flight2$Wing_length = log10(df_flight2$Wing_length)
df_flight2$Kipps_distance = log10(df_flight2$Kipps_distance)
df_flight2$Hand_wing_index = log10(df_flight2$Hand_wing_index)
df_flight2$Tail_length = log10(df_flight2$Tail_length)


df_flight2$TempRange = log10(df_flight2$TempRange)
df_flight2$AnnualPrecip = log10(df_flight2$AnnualPrecip)
df_flight2$PrecipRange = log10(df_flight2$PrecipRange)
df_flight2$Mass = log10(df_flight2$Mass)
df_flight2 = merge(df_flight2, df_need)
df_flight2$ChThSkew = as.numeric(as.character(df_flight2$ChThSkew))
row.names(df_flight2) = df_flight2$species_name
#Starting coloring tree
## extract total body length and log-transform
lnTL<-setNames(df_flight2$Tail_length,rownames(df_flight2))
head(lnTL)
## estimate ancestral states using fastAnc
fit.lnTL<-fastAnc(birds_ms_and_temp_tree,lnTL,vars=TRUE,CI=TRUE)
print(fit.lnTL,printlen=10)
## compute "contMap" object
birds_contMap<-contMap(birds_ms_and_temp_tree,lnTL,
                     plot=FALSE)
## plot "contMap" object
plot(birds_contMap,sig=2,fsize=c(0.4,0.9),
     lwd=c(2,3))

## identify the tips descended from node 102
#write.tree(birds_ms_tree, 'flying_mt_data.tre')
tips<-extract.clade(birds_ms_tree,'Node690')$tip.label #699 - peng, 690 peng + ant
tips
## prune "contMap" object to retain only these tips
pruned.contMap<-keep.tip.contMap(birds_contMap,tips)
## plot object
plot(pruned.contMap)
## add error bars
