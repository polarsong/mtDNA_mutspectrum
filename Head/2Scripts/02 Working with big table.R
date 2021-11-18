rm(list = ls(all.names = TRUE))
gc() 
install.packages('ggplot2')
install.packages('ggpubr')
library('ggplot2')
library('ggpubr')
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

#GhAhSkew for birds with different trophic levels.
brds_clsup$GhAhSkew = ((brds_clsup$neutralC - brds_clsup$neutralT)/(brds_clsup$neutralC + brds_clsup$neutralT))
GhAhgraf = ggplot(data = brds_clsup, aes(x = Gene.name, y = GhAhSkew))+
  geom_boxplot(aes(fill = TrophicLevel), notch = TRUE)
GhAhgraf + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]","[ND1]","[ND2]"))
#GhAhSkew for birds with different trophic niche
GhAhgrafTN = ggplot(data = brds_clsup, aes(x = Gene.name, y = GhAhSkew))+
  geom_boxplot(aes(fill = TrophicNiche), notch = TRUE)
GhAhgrafTN + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]","[ND1]","[ND2]"))
#GhAhSkew for birds with different foraging niche 
GhAhgrafFN = ggplot(data = brds_clsup, aes(x = Gene.name, y = GhAhSkew))+
  geom_boxplot(aes(fill = ForagingNiche), notch = TRUE)
GhAhgrafFN + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]","[ND1]","[ND2]"))
#GhAhSkew for birds with different realm
GhAhgrafR = ggplot(data = brds_clsup, aes(x = Gene.name, y = GhAhSkew))+
  geom_boxplot(aes(fill = Realm), notch = TRUE)
GhAhgrafR + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]","[ND1]","[ND2]"))


#Stg-Sac
brds_clsup$frneuA = brds_clsup$neutralA/brds_clsup$Neutral.count
brds_clsup$frneuG = brds_clsup$neutralG/brds_clsup$Neutral.count
brds_clsup$frneuC = brds_clsup$neutralC/brds_clsup$Neutral.count
brds_clsup$frneuT = brds_clsup$neutralT/brds_clsup$Neutral.count
brds_clsup$Stg = ((brds_clsup$frneuA + brds_clsup$frneuC)-(brds_clsup$frneuG+brds_clsup$frneuT))
#trophiclevel
Stlgraph = ggplot(data = brds_clsup, aes(x = Gene.name, y = Stg))+
  geom_boxplot(aes(fill = TrophicLevel), notch = TRUE)
Stlgraph + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]","[ND1]","[ND2]"))
#trophicniche
Stngraph = ggplot(data = brds_clsup, aes(x = Gene.name, y = Stg))+
  geom_boxplot(aes(fill = TrophicNiche), notch = TRUE)
Stngraph + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]","[ND1]","[ND2]"))
#foragingniche
Sfngraph = ggplot(data = brds_clsup, aes(x = Gene.name, y = Stg))+
  geom_boxplot(aes(fill = ForagingNiche), notch = TRUE)
Sfngraph + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]","[ND1]","[ND2]"))
#realm
Srgraph = ggplot(data = brds_clsup, aes(x = Gene.name, y = Stg))+
  geom_boxplot(aes(fill = Realm), notch = TRUE)
Srgraph + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]","[ND1]","[ND2]"))

#lines
brds_clsup$neutralF = brds_clsup$Neutral.count/brds_clsup$mtDNA.length

ggplot(data = brds_clsup)+
  geom_point(mapping = aes(x = Gene.name, y = neutralF))

#lines 
Atl = ggplot(data = brds_clsup, aes(x = Gene.name, y = frneuA))+ 
  geom_boxplot()
Atl = Atl + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]","[ND1]","[ND2]"))
Gtl = ggplot(data = brds_clsup, aes(x = Gene.name, y = frneuG))+ 
  geom_boxplot()
Gtl = Gtl + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]","[ND1]","[ND2]"))
Ctl = ggplot(data = brds_clsup, aes(x = Gene.name, y = frneuC))+ 
  geom_boxplot()
Ctl = Ctl + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]","[ND1]","[ND2]"))
Ttl = ggplot(data = brds_clsup, aes(x = Gene.name, y = frneuT))+ 
  geom_boxplot()
Ttl = Ttl + xlim(c("[COX1]","[COX2]","[ATP8]","[ATP6]","[COX3]", "[ND3]", "[ND4L]","[ND4]","[ND5]","[CYTB]","[ND6]","[ND1]","[ND2]"))

ggarrange(Atl, Gtl, Ctl, Ttl, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)
