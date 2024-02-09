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
library(RSQLite)
df_long = read.csv('AVES_longevity.csv')
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
df_long$scinam = firstup(df_long$scinam)
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
names(df_long) = c('Species', 'Longevity', 'Origin', 'Data')
df_long_short = merge(df_long, df_short)
df_long_short$Mass = as.numeric(df_long_short$Mass)
df_long_short$logmass = log10(df_long_short$Mass)
df_long_short$GhAhSkew = as.numeric(df_long_short$GhAhSkew)
lm_long = lm(GhAhSkew~Longevity, data=df_long_short)
modelSummary = summary(lm_long)
modelSummary
ggplot(df_long_short, aes(x = Longevity, y = GhAhSkew))+
  geom_point()+
  geom_smooth(method = 'lm', formula = y~x)+                            
  annotate('text', x = (min(df_long_short$Longevity)+1), y = min(df_long_short$GhAhSkew), label = round(modelSummary$coefficients[,4][2], digits = 4))
plot(lm_long)
shapiro.test(rstandard(lm_long))

lm_long = lm(GhAhSkew~Longevity+logmass, data=df_long_short)
modelSummary = summary(lm_long)
modelSummary
ggplot(df_long_short, aes(x = Longevity, y = GhAhSkew))+
  geom_point()+
  geom_smooth(method = 'lm', formula = y~x)+                            
  annotate('text', x = (min(df_long_short$Longevity)+1), y = min(df_long_short$GhAhSkew), label = round(modelSummary$coefficients[,4][2], digits = 4))
plot(lm_long)
shapiro.test(rstandard(lm_long))

df_long_short$loglong = log10(df_long_short$Longevity)
lm_long = lm(GhAhSkew~loglong+logmass, data=df_long_short)
modelSummary = summary(lm_long)
modelSummary
ggplot(df_long_short, aes(x = Longevity, y = GhAhSkew))+
  geom_point()+
  geom_smooth(method = 'lm', formula = y~x)+                            
  annotate('text', x = (min(df_long_short$Longevity)+1), y = min(df_long_short$GhAhSkew), label = round(modelSummary$coefficients[,4][2], digits = 4))
plot(lm_long)
shapiro.test(rstandard(lm_long))
