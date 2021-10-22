rm(list = ls(all.names = TRUE))
gc() 
sup = read.table("../../Body/3Results/Birds supplementary materials - DatabaseS1.csv", header = TRUE, sep = ',') #supplements materials
gold = read.table("../../Body/3Results/Golden birds dataset.csv", header = TRUE, sep = ',') #reading golden dataset
clsup = data.frame(sup$Binomial, sup$Realm, sup$TrophicLevel, sup$TrophicNiche, sup$ForagingNiche) #get rid of PC
names(clsup) = c('Species.name','Realm', 'TrophicLevel', 'TrophicNiche', 'ForagingNiche')
brds = data.frame(gold$Species.name)
names(brds) = c('Species.name')
brds_clsup = merge(clsup, brds)
