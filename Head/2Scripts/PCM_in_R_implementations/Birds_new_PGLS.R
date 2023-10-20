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
row.names(df_fly_peng1) = df_fly_peng1$species_name
spp1 = rownames(df_fly_peng1)
corBM1 = corBrownian(phy=peng_fly_tree1,form=~spp1)


pgls_mt_peng = gls(Mutation_AG_syn~ability_to_fly,
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
row.names(df_fly_not1) = df_fly_not1$species_name
spp2 = rownames(df_fly_not1)
corBM2 = corBrownian(phy=not_fly_tree1,form=~spp2)

row.names(df_fly_not2) = df_fly_not2$species_name
spp3 = rownames(df_fly_not2)
corBM3 = corBrownian(phy=not_fly_tree1,form=~spp3)
pgls_mt_not_fly = gls(Mutation_AG_syn~ability_to_fly,
                   data=df_fly_not2,correlation=corBM2)
summary(pgls_mt_not_fly)
#write.tree(peng_fly_tree1,'peng_and_fly_all_data.tre')
#write.tree(not_fly_tree1,'no_peng_and_fly_all_data.tre')


#dive pgls


