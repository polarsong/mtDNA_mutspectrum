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
names(df_ff1) = c('species_name', 'T>G', 'T>C', 'T>A', 'G>T', 'G>C', 'G>A', 'C>T', 'C>G', 'C>A', 'A>T', 'A>G', 'A>C')
df_fly_peng = merge(df_fly_peng, df_ff1)
#do pgls later