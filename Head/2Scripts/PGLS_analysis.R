#changing data and PGLS analysis
#1 - antarctic 0 - afrotropic
rm(list=ls(all=T))

library(dplyr)
library(ggplot2)
library(ape)
library(geiger)
library(nlme)
library(phytools)

df = read.csv('../../Body/3Results/For_Bogdan.csv') #reading file
df_aa = df[df$realm == 'Antarctic',]
df_aa$realm_num = 1
df_af = df[df$realm == 'Afrotropic',]
df_af$realm_num = 0
df_all = rbind(df_aa, df_af)
rownames(df_all) = 1:nrow(df_all)
df_tree = read.tree('../../Body/3Results/phylo.treefile')
plot(df_tree) 
name.check(df_all,df_tree)
pglsModel <- gls(realm_num ~ neutral_c, correlation = corBrownian(phy = df_tree),
                 data = df_all, method = "ML")
summary(pglsModel)
#trying smth else
df_na = df[df$realm != 'Antarctic',]
df_na$realm_num = 0
df_final = rbind(df_na, df_aa)
pglsModel <- gls(realm_num ~ neutral_c, correlation = corBrownian(phy = df_tree),
                 data = df_final, method = "ML")
