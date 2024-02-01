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
df11 = read.csv ('Species_life-histories.csv')
df12 = read.csv('GlobalBMRbase.csv')
df_mtdna = read.csv('../../Paper_materials_2024/Birds_dataset_paper.csv')
df_mtdna_orders = unique(df_mtdna[,c(2,5)])
df_mtdna_orders$taxonomy = sub('formes.*', '', df_mtdna_orders$taxonomy)
df_mtdna_orders$taxonomy = sub('.*\\ ', '', df_mtdna_orders$taxonomy)
df_mtdna_orders$taxonomy = paste(df_mtdna_orders$taxonomy, 'formes', sep = '')
names(df_mtdna_orders) = c('Species', 'Orders')
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

df12_bmr = df12[df12$Trait == 'BMR',]
names_v = unique(df12_bmr$Original_Name)
df_short1 = data.frame()
for (i in names_v)
{
  df1 = df12_bmr[df12_bmr$Original_Name == i,]
  a1 = sum(df1$TraitValue)
  a2 = nrow(df1)
  a = a1/a2
  b = 'BMR'
  ab = c(i, b, a)
  df_short1 = rbind(df_short1, ab)
}
names(df_short1) = c('Species', 'Trait', 'BMR_value')
df_short$GhAhSkew = as.numeric(df_short$GhAhSkew)

df_short2 = merge(df_short, df_short1)
df_short2$BMR_value = as.numeric(df_short2$BMR_value)
df_short2 = merge(df_short2, df_mtdna_orders)
lm1 = lm(log10(BMR_value)~GhAhSkew, data = df_short2)
abc = summary(lm1)
abc1 = abc$coefficients[,4][2]
ggplot(df_short2, aes(x = log10(BMR_value), y = GhAhSkew))+
  geom_point()+
  geom_smooth(method = 'lm', formula = y~x)+                            
  annotate('text', x = (min(log10(df_short2$BMR_value))+1), y = min(df_short2$GhAhSkew), label = round(abc1, digits = 4))
df_short2$log10BMR = log10(df_short2$BMR_value)
df_short2$Mass = as.numeric(df_short2$Mass)
df_short2$log10Mass = log10(df_short2$Mass)
counter = 0
for (i in unique(df_short2$Orders))
{
  a = df_short2[df_short2$Orders == i,]
  linearMod <- lm(log10BMR~GhAhSkew, data=a)
  modelSummary = summary(linearMod)
  c = modelSummary$coefficients[,4][2]
  b = ggplot(a, aes(x = log10BMR, y = GhAhSkew))+
    geom_point()+
    ggtitle(i)+
    geom_smooth(method = 'lm', formula = y~x)+                            
    annotate('text', x = (min(a$log10BMR)+1), y = min(a$GhAhSkew), label = round(c, digits = 4))
  assign(paste0("Order_", i), b)
  counter = counter + 1
}
pdf('Birds_BMR.pdf',         
    width = 30,height = 100)
ggarrange(Order_Accipitriformes, Order_Anseriformes, Order_Apodiformes, Order_Apterygiformes, Order_Casuariiformes, Order_Charadriiformes, Order_Coliiformes, Order_Columbiformes, Order_Coraciiformes, Order_Cuculiformes, Order_Falconiformes, Order_Galliformes, Order_Gruiformes, Order_Passeriformes, Order_Pelecaniformes, Order_Piciformes, Order_Procellariiformes, Order_Psittaciformes, Order_Sphenisciformes, Order_Strigiformes, Order_Struthioniformes, Order_Tinamiformes, Order_Upupiformes, 
          ncol = 3, nrow = 11)
dev.off()

lmMassBMR = lm(GhAhSkew~log10Mass+log10BMR , data=df_short2)
summary(lmMassBMR)
res = resid(lmMassBMR)
plot(lmMassBMR)
plot(fitted(lmMassBMR), res)
abline(0,0)
qqnorm(res)
qqline(res)
shapiro.test(rstandard(lmMassBMR))

df_short_pos = df_short2[df_short2$GhAhSkew >= 0,]
df_short_pos$log10GhAhSkew = log10(df_short_pos$GhAhSkew)
lmMassBMR1 = lm(log10GhAhSkew~log(BMR_value) + log(Mass), data=df_short_pos)
lmMassBMR12 = lm(GhAhSkew~log(BMR_value) + log(Mass), data=df_short_pos)
summary(lmMassBMR1)
summary(lmMassBMR12)
AIC(lmMassBMR12, lmMassBMR1)



for (i in unique(df_short2$Orders))
{
  a = df_short2[df_short2$Orders == i,]
  linearMod <- lm(GhAhSkew~log10BMR + log10Mass, data=a)
  modelSummary = summary(linearMod)
  c = modelSummary$coefficients[,4][2]
  b = ggplot(a, aes(x = log10BMR, y = GhAhSkew))+
    geom_point()+
    ggtitle(i)+
    geom_smooth(method = 'lm', formula = y~x)+                            
    annotate('text', x = (min(a$log10BMR)+1), y = min(a$GhAhSkew), label = round(c, digits = 4))
  assign(paste0("Order_", i), b)
  counter = counter + 1
}
pdf('Birds_BMR.pdf',         
    width = 30,height = 100)
ggarrange(Order_Accipitriformes, Order_Anseriformes, Order_Apodiformes, Order_Apterygiformes, Order_Casuariiformes, Order_Charadriiformes, Order_Coliiformes, Order_Columbiformes, Order_Coraciiformes, Order_Cuculiformes, Order_Falconiformes, Order_Galliformes, Order_Gruiformes, Order_Passeriformes, Order_Pelecaniformes, Order_Piciformes, Order_Procellariiformes, Order_Psittaciformes, Order_Sphenisciformes, Order_Strigiformes, Order_Struthioniformes, Order_Tinamiformes, Order_Upupiformes, 
          ncol = 3, nrow = 11)
dev.off()




df_mtdna11 = merge(df_short, df11, by = 'Species')
ggplot(df_mtdna11, aes(x = Migration, y = GhAhSkew))+
  geom_boxplot()
ggplot(df_mtdna11, aes(x = Habitat, y = GhAhSkew))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(df_mtdna11, aes(x = Foraging_environment, y = GhAhSkew))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(df_mtdna11, aes(x = Daily_activity, y = GhAhSkew))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(df_mtdna11, aes(x = Diet, y = GhAhSkew))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(df_mtdna11, aes(x = Clutch, y = GhAhSkew))+
  geom_point()
ggplot(df_mtdna11, aes(x = NestingPeriod, y = GhAhSkew))+
  geom_point()
