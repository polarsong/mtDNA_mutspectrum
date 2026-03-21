rm(list = ls(all=TRUE))
library(ggplot2)
library(ape); library(phytools);  library(geiger)
library(nlme)
#gene data
df_mtdna = read.csv('../Work_with_Andrey/Birds_dataset_paper.csv', header = TRUE, sep = ';')
df_mtdna$Mass = gsub(',', '.', df_mtdna$Mass)
df_mtdna$ghahSkew = gsub(',', '.', df_mtdna$ghahSkew)
df_mtdna$chthSkew = gsub(',', '.', df_mtdna$chthSkew)
df_mtdna$Mass = as.numeric(as.character(df_mtdna$Mass))
df_mtdna$ghahSkew = as.numeric(as.character(df_mtdna$ghahSkew))
df_mtdna$chthSkew = as.numeric(as.character(df_mtdna$chthSkew))
names_v = unique(df_mtdna$Species)
df_short = data.frame()
for (i in names_v)
{
  df1 = df_mtdna[df_mtdna$Species == i,]
  a = sum(df1$ghahSkew)/12
  b = sum(df1$chthSkew)/12
  v = sum(df1$Mass)/12
  ab = c(i, a, b, v)
  df_short = rbind(df_short, ab)
}
names(df_short) = c('Species', 'GhAhSkew', 'ThChSkew', 'Mass')
df_short$Species = gsub(' ', '_', df_short$Species)

#fly data
df_fly = read.csv('../flying_birds.csv')
df_fly = df_fly[,c(2,3,4)]
names(df_fly) = c('species_name', 'flightless', 'diving')
df_fly_clean1 = df_fly[df_fly$flightless =='Flightless',]
df_fly_clean= df_fly[df_fly$flightless == 'Almost_flightless',]
df_fly_clean = na.omit(df_fly_clean)
df_fly_clean1 = na.omit(df_fly_clean1)
df_fly = df_fly[df_fly$flightless != 'Flightless',]
df_fly = df_fly[df_fly$flightless != 'Almost_flightless',]
df_fly_clean$flightless = 'Tinamiformes'
df_fly_clean1$flightless = 'Casuariiformes'
df_fly_big = rbind(df_fly, df_fly_clean, df_fly_clean1)
names(df_fly_big) = c("Species", 'flightless', 'diving')
df_fly_big$Species = gsub(' ', '_', df_fly_big$Species)
df_fly_final = merge(df_fly_big, df_short)
df_fly_final = df_fly_final[df_fly_final$flightless != 'Galliformes',]
df_fly_final[df_fly_final$flightless == '0',]$flightless = 'Flying birds'
df_fly_final$flightless1 = factor(df_fly_final$flightless, levels = c('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes"))
df_fly_final$GhAhSkew = as.numeric(as.character(df_fly_final$GhAhSkew))
ggplot(df_fly_final, aes(x = flightless, y = GhAhSkew, color = flightless1))+
  geom_point(position = position_jitter(width = 0.2))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Birds groups')+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
t.test(df_fly_final[df_fly_final$flightless == 'Flying birds',]$GhAhSkew, df_fly_final[df_fly_final$flightless == "Tinamiformes" | df_fly_final$flightless == "Apterygiformes" | df_fly_final$flightless == "Casuariiformes" | df_fly_final$flightless == "Struthioniformes" | df_fly_final$flightless == "Rheiformes",]$GhAhSkew)
t.test(df_fly_final[df_fly_final$flightless == 'Flying birds',]$GhAhSkew, df_fly_final[df_fly_final$flightless == "Psittaciformes" | df_fly_final$flightless == "Columbiformes" | df_fly_final$flightless == "Eurypygiformes" | df_fly_final$flightless == "Gruiformes" | df_fly_final$flightless == "Sphenisciformes",]$GhAhSkew)
t.test(df_fly_final[df_fly_final$flightless == 'Flying birds',]$GhAhSkew, df_fly_final[df_fly_final$flightless == "Psittaciformes" | df_fly_final$flightless == "Columbiformes" | df_fly_final$flightless == "Eurypygiformes" | df_fly_final$flightless == "Gruiformes",]$GhAhSkew)
t.test(df_fly_final[df_fly_final$flightless == 'Flying birds',]$GhAhSkew, df_fly_final[df_fly_final$flightless == "Sphenisciformes",]$GhAhSkew)

#phylogenetics
df_fly_final$abtf = 1
df_fly_final[df_fly_final$flightless1 != 'Flying birds',]$abtf = 0
feathertree <- read.nexus("../Work_with_Andrey/Ultrametric_feathertree.nex")
feathertree$node.label <- NULL
row.names(df_fly_final) = df_fly_final$Species
name.check(feathertree, df_fly_final)
df_fly_final[df_fly_final$Species == "Agapornis_pullarius" | df_fly_final$Species == "Mergus_squamatus" | df_fly_final$Species == "Vestiaria_coccinea",] = NA
df_fly_final = na.omit(df_fly_final)
listSkew_fly = df_fly_final$Species
listTree_fly <- feathertree$tip.label
SpeciesToDrop <- setdiff(listTree_fly, listSkew_fly)
drop.tip(feathertree, SpeciesToDrop) -> fly_tree
name.check(feathertree, df_short)
#with penguins
spp_1 = rownames(df_fly_final)
corLambda_1 = corPagel(value = 1, phy = fly_tree, form=~spp_1)
pgls_1 = gls(GhAhSkew~abtf,
                data=df_fly_final, correlation=corLambda_1)
summary(pgls_1)
#drop penguins
df_paleo = df_fly_final[df_fly_final$flightless != 'Sphenisciformes',]
name.check(fly_tree, df_paleo)
listSkew_paleo = df_paleo$Species
listTree_paleo <- fly_tree$tip.label
SpeciesToDrop <- setdiff(listTree_paleo, listSkew_paleo)
drop.tip(fly_tree, SpeciesToDrop) -> paleo_tree
spp_2 = rownames(df_paleo)
corLambda_2 = corPagel(value = 1, phy = paleo_tree, form=~spp_2)
pgls_2 = gls(GhAhSkew~abtf,
             data=df_paleo, correlation=corLambda_2)
summary(pgls_2)

paleo = setNames(df_paleo[,"abtf"], rownames(df_paleo))
Gh_paleo = setNames(df_paleo[,"GhAhSkew"], rownames(df_paleo))
paleopic = pic(paleo, paleo_tree)
Gh_paleoPIC = pic(Gh_paleo, paleo_tree)
fit_pic_paleo = lm(paleopic~Gh_paleoPIC+0)
fit_pic_paleo
summary(fit_pic_paleo)
plot(Gh_paleoPIC~paleopic)

