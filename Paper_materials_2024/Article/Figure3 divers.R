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
  pro = sum(df1$CCT) + sum(df1$CCA) + sum(df1$CCG) + sum(df1$CCC) 
  pheleu = sum(df1$TTT) + sum(df1$TTC) + sum(df1$TTG) + sum(df1$TTA) 
  ab = c(i, a, b, v, pro, pheleu)
  df_short = rbind(df_short, ab)
}
names(df_short) = c('Species', 'GhAhSkew', 'ThChSkew', 'Mass', 'Pro', 'PheLeu')
df_short$Species = gsub(' ', '_', df_short$Species)
#dive data
df_fly = read.csv('../flying_birds.csv')
df_fly = df_fly[,c(2,3,4)]
names(df_fly) = c('Species', 'flightless', 'diving')
df_fly$Species = gsub(' ', '_', df_fly$Species)
df_dive_final = merge(df_fly, df_short, by = 'Species')
df_dive_final = df_dive_final[df_dive_final$diving != 'waterbird',]
df_dive_final[df_dive_final$diving == '0',]$diving = 'Non-diving birds'
df_dive_final$GhAhSkew = as.numeric(as.character(df_dive_final$GhAhSkew))
df_dive_final$diving1 = factor(df_dive_final$diving, levels = c('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes"))
ggplot(df_dive_final, aes(x = diving, y = GhAhSkew, colour = diving1))+
  geom_point(position = position_jitter(width = 0.2))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Birds groups')+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
t.test(df_dive_final[df_dive_final$diving == "Non-diving birds",]$GhAhSkew, df_dive_final[df_dive_final$flightless == "Anseriformes" | df_dive_final$flightless == "Sphenisciformes" | df_dive_final$flightless == "Podicipediformes" | df_dive_final$flightless == "Gaviiformes" | df_dive_final$flightless == "Suliformes" | df_dive_final$flightless == "Coraciiformes" | df_dive_final$flightless == "Passeriformes" | df_dive_final$flightless == "Gruiformes" | df_dive_final$flightless == "Charadriiformes" | df_dive_final$flightless == "Procellariiformes",]$GhAhSkew)
t.test(df_dive_final[df_dive_final$diving == "Non-diving birds",]$GhAhSkew, df_dive_final[df_dive_final$flightless == "Anseriformes" | df_dive_final$flightless == "Sphenisciformes" | df_dive_final$flightless == "Podicipediformes" | df_dive_final$flightless == "Gaviiformes" | df_dive_final$flightless == "Suliformes",]$GhAhSkew)
t.test(df_dive_final[df_dive_final$diving == "Non-diving birds",]$GhAhSkew, df_dive_final[df_dive_final$flightless == "Coraciiformes" | df_dive_final$flightless == "Passeriformes" | df_dive_final$flightless == "Gruiformes" | df_dive_final$flightless == "Charadriiformes" | df_dive_final$flightless == "Procellariiformes",]$GhAhSkew)

#phylogenetics
feathertree <- read.nexus("../Work_with_Andrey/Ultrametric_feathertree.nex")
feathertree$node.label <- NULL
row.names(df_dive_final) = df_dive_final$Species
name.check(feathertree, df_dive_final)
df_dive_final[df_dive_final$Species == "Agapornis_pullarius" | df_dive_final$Species == "Mergus_squamatus" | df_dive_final$Species == "Vestiaria_coccinea",] = NA
df_dive_final = na.omit(df_dive_final)
listSkew_dive = df_dive_final$Species
listTree_dive <- feathertree$tip.label
SpeciesToDrop <- setdiff(listTree_dive, listSkew_dive)
drop.tip(feathertree, SpeciesToDrop) -> dive_tree
name.check(dive_tree, df_dive_final)
df_dive_final$abtd = 0
df_dive_final[df_dive_final$diving != 'Non-diving birds',]$abtd = 1
spp_1 = rownames(df_dive_final)
corLambda_1 = corPagel(value = 1, phy = dive_tree, form=~spp_1)
pgls_1 = gls(GhAhSkew~abtd,
             data=df_dive_final, correlation=corLambda_1)
summary(pgls_1)

#cut some divers
df_dive_cut = df_dive_final[df_dive_final$diving != "Coraciiformes" & df_dive_final$diving != "Passeriformes" & df_dive_final$diving != "Gruiformes" & df_dive_final$diving != "Charadriiformes" & df_dive_final$diving != "Procellariiformes",]
name.check(dive_tree, df_dive_cut)
listSkew_dive1 = df_dive_cut$Species
listTree_dive1 <- dive_tree$tip.label
SpeciesToDrop <- setdiff(listTree_dive1, listSkew_dive1)
drop.tip(dive_tree, SpeciesToDrop) -> dive_tree_cut
name.check(dive_tree_cut, df_dive_cut)
spp_2 = rownames(df_dive_cut)
corLambda_2 = corPagel(value = 1, phy = dive_tree_cut, form=~spp_2)
pgls_2 = gls(GhAhSkew~abtd,
             data=df_dive_cut, correlation=corLambda_2)
summary(pgls_2)

#AA shift
#AA shift
df_dive_cut$Pro = as.numeric(as.character(df_dive_cut$Pro))
df_dive_cut$PheLeu = as.numeric(as.character(df_dive_cut$PheLeu))
df_dive_cut$propheleu = df_dive_cut$Pro/df_dive_cut$PheLeu
ggplot(df_dive_cut, aes(x = diving, y = propheleu))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes")
