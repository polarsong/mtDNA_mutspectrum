#new PGLS for birds
rm(list = ls(all=TRUE))
library(phytools)
library(nlme)
library(geiger)
#flying birds
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
df_fly_peng = df_flight[df_flight$ability_to_fly == 'Flying' | df_flight$ability_to_fly =='Sphenisciformes',]
df_fly_not = df_flight[df_flight$ability_to_fly !='Sphenisciformes',]

peng_fly_tree = read.tree('flying_and_peng.tre')
no_peng_fly_tree = read.tree('flying_and_no_peng.tre')
spp = rownames(df_fly_peng)
corBM = corBrownian(phy=peng_fly_tree,form=~spp)

pgls_flying = gls(GhAhSkew~ability_to_fly,
                  data=df_fly_peng,correlation=corBM)
summary(pgls_flying)
anova(pgls_flying)

df_fly_not[df_fly_not$ability_to_fly != 'Flying',]$ability_to_fly = 'Non-flying'
spp = rownames(df_fly_not)
corBM = corBrownian(phy=no_peng_fly_tree,form=~spp)

pgls_flying = gls(GhAhSkew~ability_to_fly,
                  data=df_fly_not,correlation=corBM)
summary(pgls_flying)
anova(pgls_flying)





#mutspec PGLS
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
df_fly_peng1 = merge(df_fly_peng, df_ff1)
df_fly_peng1 = merge(df_fly_peng1, df_syn1)
need_species = setdiff(df_fly_peng$species_name, df_fly_peng1$species_name)
correct_need_species = setdiff(df_fly_peng$species_name, need_species)
peng_fly_tree1<-keep.tip(peng_fly_tree,correct_need_species)
df_fly_peng1$Mutation_AG = as.numeric(df_fly_peng1$Mutation_AG)
df_fly_peng1$Mutation_AG_syn = as.numeric(df_fly_peng1$Mutation_AG_syn)
df_fly_peng1$Mutation_CT = as.numeric(df_fly_peng1$Mutation_CT)
df_fly_peng1$Mutation_CT_syn = as.numeric(df_fly_peng1$Mutation_CT_syn)
df_fly_peng1$ability_to_fly_binar = 0
df_fly_peng1[df_fly_peng1$ability_to_fly == 'Sphenisciformes',]$ability_to_fly_binar = 1
df_fly_peng1$ability_to_fly_binar_extra = 1
df_fly_peng1[df_fly_peng1$ability_to_fly == 'Sphenisciformes',]$ability_to_fly_binar_extra = 0
row.names(df_fly_peng1) = df_fly_peng1$species_name
spp1 = rownames(df_fly_peng1)
corBM1 = corBrownian(phy=peng_fly_tree1,form=~spp1)


pgls_mt_peng = gls(GhAhSkew~ability_to_fly_binar_extra,
                  data=df_fly_peng1,correlation=corBM1)
summary(pgls_mt_peng)


df_fly_not1 = merge(df_fly_not, df_ff1)
df_fly_not1 = merge(df_fly_not1, df_syn1)
need_species = setdiff(df_fly_not$species_name, df_fly_not1$species_name)
correct_need_species = setdiff(df_fly_not$species_name, need_species)
not_fly_tree = read.tree('flying_and_no_peng.tre')
not_fly_tree1<-keep.tip(not_fly_tree,correct_need_species)
df_fly_not1$Mutation_AG = as.numeric(df_fly_not1$Mutation_AG)
df_fly_not1$Mutation_AG_syn = as.numeric(df_fly_not1$Mutation_AG_syn)
df_fly_not1$Mutation_CT = as.numeric(df_fly_not1$Mutation_CT)
df_fly_not1$Mutation_CT_syn = as.numeric(df_fly_not1$Mutation_CT_syn)
df_fly_not2 = df_fly_not1
df_fly_not2[df_fly_not2$ability_to_fly != 'Flying',]$ability_to_fly = 'Non-flying'
df_fly_not1$ability_to_fly_binar = 0
df_fly_not1[df_fly_not1$ability_to_fly != 'Flying',]$ability_to_fly_binar = 1
df_fly_not1$ability_to_fly_binar_extra = 1
df_fly_not1[df_fly_not1$ability_to_fly != 'Flying',]$ability_to_fly_binar_extra = 0
row.names(df_fly_not1) = df_fly_not1$species_name
spp2 = rownames(df_fly_not1)
corBM2 = corBrownian(phy=not_fly_tree1,form=~spp2)
pgls_mt_not_fly = gls(GhAhSkew~ability_to_fly_binar_extra,
                      data=df_fly_not1,correlation=corBM2)
summary(pgls_mt_not_fly)

row.names(df_fly_not2) = df_fly_not2$species_name
spp3 = rownames(df_fly_not2)
corBM3 = corBrownian(phy=not_fly_tree1,form=~spp3)
pgls_mt_not_fly = gls(GhAhSkew~ability_to_fly,
                   data=df_fly_not2,correlation=corBM3)
summary(pgls_mt_not_fly)
#write.tree(peng_fly_tree1,'peng_and_fly_all_data.tre')
#write.tree(not_fly_tree1,'no_peng_and_fly_all_data.tre')


#dive pgls
df_dive = read.csv('dive_and_gene.csv')
df_dive = df_dive[,c(2:6)]
dive_tree = read.tree('diving_birds_tree.tre')
row.names(df_dive) = df_dive$species_name
df_dive1 = merge(df_dive, df_syn1)
unique(df_dive1$ability_to_dive)
df_dive2 = df_dive1
df_dive1[df_dive1$ability_to_dive == 'Anseriformes',]$ability_to_dive = 'Diving_all_species'
df_dive1[df_dive1$ability_to_dive == "Sphenisciformes",]$ability_to_dive = 'Diving_all_species'
df_dive1[df_dive1$ability_to_dive == "Podicipediformes",]$ability_to_dive = 'Diving_all_species'
df_dive1[df_dive1$ability_to_dive == "Gaviiformes",]$ability_to_dive = 'Diving_all_species'
df_dive1[df_dive1$ability_to_dive == "Suliformes",]$ability_to_dive = 'Diving_all_species'

df_dive1[df_dive1$ability_to_dive == "Charadriiformes",]$ability_to_dive = 'Diving_some_species'
df_dive1[df_dive1$ability_to_dive == "Coraciiformes",]$ability_to_dive = 'Diving_some_species'
df_dive1[df_dive1$ability_to_dive == "Passeriformes",]$ability_to_dive = 'Diving_some_species'
df_dive1[df_dive1$ability_to_dive == "Procellariiformes",]$ability_to_dive = 'Diving_some_species'
df_dive1[df_dive1$ability_to_dive == "Gruiformes",]$ability_to_dive = 'Diving_some_species'
df_dive1$GhAhSkew = as.numeric(as.character(df_dive1$GhAhSkew))
df_dive1$Mutation_AG_syn = as.numeric(as.character(df_dive1$Mutation_AG_syn))
row.names(df_dive1) = df_dive1$species_name
need_species = setdiff(df_dive$species_name, df_dive1$species_name)
correct_need_species = setdiff(df_dive$species_name, need_species)
dive_tree1 = keep.tip(dive_tree,correct_need_species)
df_dive1$ability_to_dive_binar = 0
df_dive1[df_dive1$ability_to_dive != 'Non-diving',]$ability_to_dive_binar = 1

spp4 = rownames(df_dive1)
corBM4 = corBrownian(phy=dive_tree1,form=~spp4)
pgls_dive = gls(GhAhSkew~ability_to_dive_binar,
                      data=df_dive1,correlation=corBM4)
summary(pgls_dive)



df_dive11 = df_dive1[df_dive1$ability_to_dive != 'Diving_some_species',]
row.names(df_dive11) = df_dive11$species_name
need_species = setdiff(df_dive1$species_name, df_dive11$species_name)
correct_need_species = setdiff(df_dive1$species_name, need_species)
dive_tree11 = keep.tip(dive_tree1,correct_need_species)

spp5 = rownames(df_dive11)
corBM5 = corBrownian(phy=dive_tree11,form=~spp5)
pgls_dive1 = gls(GhAhSkew~ability_to_dive_binar,
                data=df_dive11,correlation=corBM5)
summary(pgls_dive1)

df_dive12 = df_dive1[df_dive1$ability_to_dive != 'Diving_all_species',]
row.names(df_dive12) = df_dive12$species_name
need_species = setdiff(df_dive1$species_name, df_dive12$species_name)
correct_need_species = setdiff(df_dive1$species_name, need_species)
dive_tree12 = keep.tip(dive_tree1,correct_need_species)

spp6 = rownames(df_dive12)
corBM6 = corBrownian(phy=dive_tree12,form=~spp6)
pgls_dive2 = gls(GhAhSkew~ability_to_dive_binar,
                 data=df_dive12,correlation=corBM6)
summary(pgls_dive2)

#TR_TS
df_tr_ts = read.csv('TR_TS_dataset.csv')
df_tr_ts = df_tr_ts[,c(2:19)]
df_tr_ts = df_tr_ts[,c(1,16:18)]
