rm(list = ls(all.names = TRUE))
gc() 
install.packages('ggplot2')
library('ggplot2')
sup = read.table("../../Body/3Results/Birds supplementary materials - DatabaseS1.csv", header = TRUE, sep = ',') #supplements materials
gold = read.table("../../Body/3Results/Golden birds dataset.csv", header = TRUE, sep = ',') #reading golden dataset
clsup = data.frame(sup$Binomial, sup$Realm, sup$TrophicLevel, sup$TrophicNiche, sup$ForagingNiche) #get rid of PC
names(clsup) = c('Species.name','Realm', 'TrophicLevel', 'TrophicNiche', 'ForagingNiche')
clsup$Species.name = gsub("_", " ", clsup$Species.name)
names(brds) = c('Species.name')
brds_clsup = merge(clsup, gold) #merge golden dataset and supplements materials

#comparing neutral nucleotides in genes
#neutralA graph
grafA = ggplot(data = brds_clsup, aes(x = Gene.name, y = neutralA))+
  geom_boxplot(aes(fill = TrophicLevel), notch = TRUE)
grafA + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]"))
#neutralG graph
grafG = ggplot(data = brds_clsup, aes(x = Gene.name, y = neutralG))+
  geom_boxplot(aes(fill = TrophicLevel), notch = TRUE)
grafG + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]"))
#neutralC graph
grafC = ggplot(data = brds_clsup, aes(x = Gene.name, y = neutralC))+
  geom_boxplot(aes(fill = TrophicLevel),notch = TRUE)
grafC + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]"))
#neutralT graph
grafT = ggplot(data = brds_clsup, aes(x = Gene.name, y = neutralT))+
  geom_boxplot(aes(fill = TrophicLevel),notch = TRUE)
grafT + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]"))


brds_clsup$GhAhSkew = ((brds_clsup$neutralC - brds_clsup$neutralT)/(brds_clsup$neutralC + brds_clsup$neutralT))
GhAhgraf = ggplot(data = brds_clsup, aes(x = Gene.name, y = GhAhSkew))+
  geom_boxplot(aes(fill = TrophicLevel), notch = TRUE)
GhAhgraf + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]"))
