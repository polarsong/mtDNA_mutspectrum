rm(list = ls(all=TRUE))
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggbiplot)
library(phytools)
library(nlme)
library(geiger)
library(ggtree)
library(stringr)
library(dplyr)
df_mtdna = read.csv('../../Paper_materials_2024/Birds_dataset_paper.csv')
names_v = unique(df_mtdna$species_name)
df_mtdna$Mass = as.numeric(as.character(df_mtdna$Mass))
df_short = data.frame()
for (i in names_v)
{
  df1 = df_mtdna[df_mtdna$species_name == i,]
  a = sum(df1$ghahSkew)/12
  b = sum(df1$chthSkew)/12
  v = sum(df1$Mass)/12
  ab = c(i, a, b, v)
  df_short = rbind(df_short, ab)
}
names(df_short) = c('Species', 'GhAhSkew', 'ThChSkew', 'Mass')
df_short$GhAhSkew = as.numeric(df_short$GhAhSkew)
df_temp = read.csv('temp_new_data.csv')
df_big = merge(df_short, df_temp, by = 'Species')

lm_temp = lm(GhAhSkew~TNZ, data=df_big)
modelSummary = summary(lm_temp)

ggplot(df_big, aes(x = TNZ, y = GhAhSkew))+
  geom_point()+
  geom_smooth(method = 'lm', formula = y~x)+                            
  annotate('text', x = (min(df_big$TNZ)+1), y = min(df_big$GhAhSkew), label = round(modelSummary$coefficients[,4][2], digits = 4))
plot(lm_temp)
shapiro.test(rstandard(lm_temp))
