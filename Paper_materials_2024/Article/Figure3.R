rm(list = ls(all=TRUE))
library(ape); library(phytools);  library(geiger)
df_nd6 = read.csv('../Birds_mtDNA_data.csv')
df_nd6$GhAhSkew = (df_nd6$neutral_c- df_nd6$neutral_T)/(df_nd6$neutral_c + df_nd6$neutral_T)
df_nd6$ThChSkew = (df_nd6$neutral_A - df_nd6$neutral_g)/(df_nd6$neutral_A + df_nd6$neutral_g)

df_need = data.frame()
for (i in unique(df_nd6$species_name))
{
  a = df_nd6[df_nd6$species_name == i,]
  b = sum(a$GhAhSkew)/12
  ab = c(i,b)
  df_need = rbind(df_need, ab)
}
names(df_need) = c('species_name', 'GhAhSkew')
df_fly = read.csv('../flying_birds.csv')
df_fly1 = df_fly[df_fly$Flightless != 'Almost_flightless',]
df_fly1 = df_fly1[df_fly1$Flightless != 'Galliformes',]
df_fly1 = df_fly1[df_fly1$Flightless != 'Flightless',]
df_fly1 = df_fly1[,c(2,3,4)]
names(df_fly1) = c('species_name', 'ability to fly', 'ability to dive')
df_tree = merge(df_fly1, df_need)
#merge data, grab tree and go
#old KG tree
tree = read.tree('../../Paper_materials_2024/anc_kg.treefile')
df_tree$species_name = gsub(' ', '_', df_tree$species_name)
row.names(df_tree) = df_tree$species_name
name.check(tree, df_tree)
df_tree[df_tree$species_name == "Agapornis_pullarius" | df_tree$species_name == "Mergus_squamatus",] = NA
df_tree = na.omit(df_tree)
listSkew = df_tree$species_name
listTree <- tree$tip.label
SpeciesToDrop <- setdiff(listTree, listSkew)
drop.tip(tree, SpeciesToDrop) -> nonf_tree
df_tree$GhAhSkew = as.numeric(as.character(df_tree$GhAhSkew))
lnTL<-setNames(df_tree$GhAhSkew,rownames(df_tree))
fit.lnTL<-fastAnc(nonf_tree,lnTL,vars=TRUE,CI=TRUE)
print(fit.lnTL,printlen=10)
nonf_birds_contMap<-contMap(nonf_tree,lnTL,
                       plot=FALSE)
plot(nonf_birds_contMap,sig=2,fsize=c(0.45,0.9),
     lwd=c(2,3))
tips<-extract.clade(nonf_tree,'Node690')$tip.label #699 - peng, 690 peng + ant 582 non-flying 496 ducks
tips
pruned.contMap<-keep.tip.contMap(nonf_birds_contMap,tips)
plot(pruned.contMap)

#Andrey tree
feathertree <- read.nexus("../Work_with_Andrey/Ultrametric_feathertree.nex")
df_need$species_name = gsub(' ', '_', df_need$species_name)
rownames(df_need) = df_need$species_name
name.check(feathertree, df_need)
df_need[df_need$species_name == "Agapornis_pullarius" | df_need$species_name == "Mergus_squamatus" | df_need$species_name == "Coturnix_chinensis" | df_need$species_name == "Serinus_albogularis" | df_need$species_name == "Vestiaria_coccinea",] = NA
df_need = na.omit(df_need)
listSkew = df_need$species_name
listTree <- feathertree$tip.label
SpeciesToDrop <- setdiff(listTree, listSkew)
drop.tip(feathertree, SpeciesToDrop) -> big_tree

df_need$GhAhSkew = as.numeric(as.character(df_need$GhAhSkew))
lnTL_big<-setNames(df_need$GhAhSkew,rownames(df_need))
fit.lnTL_big<-fastAnc(big_tree,lnTL_big,vars=TRUE,CI=TRUE)
print(fit.lnTL_big,printlen=10)
big_birds_contMap<-contMap(big_tree,lnTL_big,
                            plot=FALSE)
plot(big_birds_contMap,sig=2,fsize=c(0.45,0.9),
     lwd=c(2,3))
tips<-extract.clade(big_tree,'I1497')$tip.label #699 - peng, 690 peng + ant 582 non-flying 496 ducks
tips
pruned.contMap<-keep.tip.contMap(nonf_birds_contMap,tips)
plot(pruned.contMap)