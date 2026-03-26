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
df_dive_final$exp_div = as.factor(df_dive_final$abtd)
spp_1 = rownames(df_dive_final)
corLambda_1 = corPagel(value = 1, phy = dive_tree, form=~spp_1)
pgls_1 = gls(GhAhSkew~abtd,
             data=df_dive_final, correlation=corLambda_1)
summary(pgls_1)
pgls_1_1 = gls(GhAhSkew~exp_div,
             data=df_dive_final, correlation=corLambda_1)
summary(pgls_1_1)

df_dive_final$loggh = log10(df_dive_final$GhAhSkew + 0.3)
pgls_1_2 = gls(loggh~abtd,
             data=df_dive_final, correlation=corLambda_1)
summary(pgls_1_2)

pgls_1_2_1 = gls(loggh~exp_div,
               data=df_dive_final, correlation=corLambda_1)
summary(pgls_1_2_1)

df_dive_final$Pro = as.numeric(as.character(df_dive_final$Pro))
df_dive_final$PheLeu = as.numeric(as.character(df_dive_final$PheLeu))
df_dive_final$propheleu = df_dive_final$Pro/df_dive_final$PheLeu
pgls_1_3 = gls(propheleu~abtd,
               data=df_dive_final, correlation=corLambda_1)
summary(pgls_1_3)

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
pgls_2_1 = gls(loggh~abtd,
             data=df_dive_cut, correlation=corLambda_2)
summary(pgls_2_1)
pgls_2_2 = gls(propheleu~abtd,
               data=df_dive_cut, correlation=corLambda_2)
summary(pgls_2_2)
df_dive_cut$logaa = log10(df_dive_cut$propheleu)
pgls_2_3 = gls(logaa~abtd,
               data=df_dive_cut, correlation=corLambda_2)
summary(pgls_2_3)
pgls_2_4 = gls(GhAhSkew~exp_div,
               data=df_dive_cut, correlation=corLambda_2)
summary(pgls_2_4)
pgls_2_5 = gls(loggh~exp_div,
               data=df_dive_cut, correlation=corLambda_2)
summary(pgls_2_5)
#AA shift
#AA shift
df_dive_cut$Pro = as.numeric(as.character(df_dive_cut$Pro))
df_dive_cut$PheLeu = as.numeric(as.character(df_dive_cut$PheLeu))
df_dive_cut$propheleu = df_dive_cut$Pro/df_dive_cut$PheLeu
ggplot(df_dive_cut, aes(x = diving, y = propheleu))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes")

#big pgls
df_dive_clean = df_dive_final[df_dive_final$flightless =='Flightless',]
df_dive_clean1= df_dive_final[df_dive_final$flightless == 'Almost_flightless',]
df_dive_clean = na.omit(df_dive_clean)
df_dive_clean1 = na.omit(df_dive_clean1)
df_divegls = df_dive_final[df_dive_final$flightless != 'Flightless',]
df_divegls = df_divegls[df_divegls$flightless != 'Almost_flightless',]
df_dive_clean$flightless = 'Tinamiformes'
df_dive_clean1$flightless = 'Casuariiformes'
df_divegls = rbind(df_divegls, df_dive_clean, df_dive_clean1)
df_divegls = df_divegls[df_divegls$flightless != 'Galliformes',]
df_divegls[df_divegls$flightless == '0',]$flightless = 'Flying birds'
df_divegls$abtf = 1
df_divegls[df_divegls$flightless != 'Flying birds',]$abtf = 0
df_divegls$exp_flight = as.factor(df_divegls$abtf)
df_divegls$Pro = as.numeric(as.character(df_divegls$Pro))
df_divegls$PheLeu = as.numeric(as.character(df_divegls$PheLeu))
df_divegls$propheleu = df_divegls$Pro/df_divegls$PheLeu
df_divegls$Mass = as.numeric(as.character(df_divegls$Mass))
name.check(df_divegls, dive_tree)
df_divegls$loggh = log10(df_divegls$GhAhSkew + 0.3)
df_divegls$logmass = log10(df_divegls$Mass)
df_divegls$logppl = log10(df_divegls$propheleu)
#listSkew_dive = df_dive_final$Species
#listTree_dive <- feathertree$tip.label
#SpeciesToDrop <- setdiff(listTree_dive, listSkew_dive)
#drop.tip(feathertree, SpeciesToDrop) -> dive_tree
spp_3 = rownames(df_divegls)
corLambda_3 = corPagel(value = 1, phy = dive_tree, form=~spp_3)
corBM<-corBrownian(phy=dive_tree,form=~spp_3)
pgls_3 = gls(GhAhSkew~abtd+abtf+Mass,
             data=df_divegls, correlation=corBM)
summary(pgls_3)
anova(pgls_3)
pgls_3_1 = gls(GhAhSkew~abtd+abtf+Mass,
             data=df_divegls, correlation=corLambda_3)
summary(pgls_3_1)
pgls_3_1_1 = gls(GhAhSkew~exp_div+exp_flight+Mass,
               data=df_divegls, correlation=corLambda_3)
summary(pgls_3_1_1)
anova(pgls_3_1)


pgls_3_2 = gls(loggh~abtd+abtf+logmass,
               data=df_divegls, correlation=corLambda_3)
summary(pgls_3_2)
pgls_3_2_1 = gls(loggh~exp_div+exp_flight+logmass,
               data=df_divegls, correlation=corLambda_3)
summary(pgls_3_2_1)
anova(pgls_3_2)
pgls_3_3 = gls(loggh~abtd+abtf+logmass,
               data=df_divegls, correlation=corBM)
summary(pgls_3_3)
anova(pgls_3_3)

pgls_3_4 = gls(GhAhSkew~abtd*abtf+Mass,
               data=df_divegls, correlation=corLambda_3)
summary(pgls_3_4)
anova(pgls_3_4)

pgls_3_5 = gls(loggh~abtd*abtf+logmass,
               data=df_divegls, correlation=corLambda_3)
summary(pgls_3_5)
anova(pgls_3_5)

pgls_3_5_1 = gls(loggh~exp_div*exp_flight+logmass,
               data=df_divegls, correlation=corLambda_3)
summary(pgls_3_5_1)



pgls_3_6 = gls(GhAhSkew~abtd*abtf,
               data=df_divegls, correlation=corLambda_3)
summary(pgls_3_6)
pgls_3_6_1 = gls(GhAhSkew~exp_div*exp_flight,
               data=df_divegls, correlation=corLambda_3)
summary(pgls_3_6_1)


pgls_3_7 = gls(loggh~abtd*abtf,
               data=df_divegls, correlation=corLambda_3)
summary(pgls_3_7)
anova(pgls_3_7)

pgls_3_7_1 = gls(loggh~exp_div*exp_flight,
               data=df_divegls, correlation=corLambda_3)
summary(pgls_3_7_1)

#cut again
dfglscut = df_divegls[df_divegls$diving != "Coraciiformes" & df_divegls$diving != "Passeriformes" & df_divegls$diving != "Gruiformes" & df_divegls$diving != "Charadriiformes" & df_divegls$diving != "Procellariiformes",]
name.check(dfglscut, dive_tree)
spp_4 = rownames(dfglscut)
corLambda_4 = corPagel(value = 1, phy = dive_tree, form=~spp_4)
pgls_4 = gls(GhAhSkew~abtd*abtf,
             data=dfglscut, correlation=corLambda_4)
summary(pgls_4)
pgls_4_1 = gls(GhAhSkew~exp_div*exp_flight,
             data=dfglscut, correlation=corLambda_4)
summary(pgls_4_1)
#Andrey graph
library(sjPlot)
plot_model(pgls_3_7, type = "int", terms = c("abtd", "abtf"))
plot_model(pgls_3_7, type = "int", terms = c("abtf", "abtd"))

df_divegls$exp_div = factor(df_divegls$abtd, levels = c(1, 0))
df_divegls$exp_flight = factor(df_divegls$abtf, levels = c(1, 0))
pgls_3_8 = gls(GhAhSkew~exp_div*exp_flight,
               data=df_divegls, correlation=corLambda_3)
summary(pgls_3_8)
plot_model(pgls_3_8, type = "int", terms = c("exp_div", "exp_flight"))
plot_model(pgls_3_8, type = "int", terms = c("exp_flight", "exp_div"))
plot(pgls_3_8)

pgls_3_9 = gls(GhAhSkew~exp_div*exp_flight+logmass,
               data=df_divegls, correlation=corLambda_3)
summary(pgls_3_9)
anova(pgls_3_9)

#FD shenanigans
df_divegls$FD = 1


#cut sphe - doesn't work
df_divegls_exp = df_divegls[df_divegls$diving != 'Sphenisciformes',]
name.check(df_divegls_exp, dive_tree)
spp_exp = rownames(df_divegls_exp)
corLambda_exp = corPagel(value = 1, phy = dive_tree, form=~spp_exp)
pgls_exp = gls(GhAhSkew~abtf*abtd+Mass,
             data=df_divegls_exp, correlation=corLambda_exp)
summary(pgls_exp)
anova(pgls_exp)

#propheleu 
pgls_4 = gls(propheleu~abtf*abtd+Mass,
             data = df_divegls, correlation = corLambda_3)
summary(pgls_4)
anova(pgls_4)
pgls_4_1 = gls(logppl~abtf*abtd+logmass,
             data = df_divegls, correlation = corLambda_3)
summary(pgls_4_1)
anova(pgls_4_1)

pgls_4_2 = gls(logppl~abtf*abtd,
               data = df_divegls, correlation = corLambda_3)
summary(pgls_4_2)
anova(pgls_4_2)
