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
ggplot(data = brds_clsup, aes(x = Gene.name, y = neutralA))+
  geom_boxplot(notch = TRUE)
ggplot(data = brds_clsup, aes(x = Gene.name, y = neutralG))+
  geom_boxplot(notch = TRUE)
ggplot(data = brds_clsup, aes(x = Gene.name, y = neutralC))+
  geom_boxplot(notch = TRUE)
ggplot(data = brds_clsup, aes(x = Gene.name, y = neutralT))+
  geom_boxplot(notch = TRUE)

ggplot(data = brds_clsup, aes(x = Gene.name , y = neutralA))+
  geom_boxplot(aes(fill = TrophicLevel), notch = TRUE)

brds_clsup$Gene.name
