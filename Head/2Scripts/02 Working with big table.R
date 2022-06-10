rm(list = ls(all.names = TRUE))
gc() 
install.packages('ggplot2')
install.packages('ggpubr')
library('ggplot2')
library('ggpubr')
sup = read.table("../../Body/3Results/Birds supplementary materials - DatabaseS1.csv", header = TRUE, sep = ',') #supplements materials
gold = read.table("../../Body/3Results/birds_list.csv", header = TRUE, sep = ',') #reading golden dataset
clsup = data.frame(sup$Binomial, sup$Realm, sup$TrophicLevel, sup$TrophicNiche, sup$ForagingNiche) #get rid of PC
names(clsup) = c('Species.name','Realm', 'TrophicLevel', 'TrophicNiche', 'ForagingNiche')
clsup$Species.name = gsub("_", " ", clsup$Species.name)
names(brds) = c('Species.name')
brds_clsup = merge(clsup, gold) #merge golden dataset and supplements materials
brds_clsup$temp = 1
df2 = aggregate(brds_clsup$temp,by = list(brds_clsup$Species.name), FUN = sum)
df3 = brds_clsup[brds_clsup$Species.name == 'Sarothrura ayresi',]

brds_clsup = read.csv("../../Body/3Results/For_Bogdan.csv")
#comparing neutral nucleotides in genes
#neutralA graph
grafA = ggplot(data = brds_clsup, aes(x = gene_name, y = neutral_A))+
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
brds_clsup = brds_clsup[!(brds_clsup$Species.name == 'Paradoxornis heudei'),]  

write.csv(brds_clsup, "final_birds_list_with_no_mistakes.csv")



tree = read.csv("../../Body/3Results/skeleton.csv")
names(tree) = c('Species.name')
tree$Species.name = gsub("_", " ", tree$Species.name)
brd_tree = merge(brds_clsup, tree)


chdf = read.csv("../../Body/3Results/For_Bogdan.csv")
chdf$GhAhSkew = ((chdf$neutral_c - chdf$neutral_T)/(chdf$neutral_c + chdf$neutral_T))
GhAhgraf = ggplot(data = chdf, aes(x = gene_name, y = GhAhSkew))+
  geom_boxplot(aes(fill = trophic_level), notch = TRUE)
GhAhgraf + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
GhAhgraf1 = ggplot(data = chdf, aes(x = gene_name, y = GhAhSkew))+
  geom_boxplot(aes(fill = realm), notch = TRUE)
GhAhgraf1 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))


chdf$frneuA = chdf$neutral_A/chdf$neutral_amount
chdf$frneuG = chdf$neutral_g/chdf$neutral_amount
chdf$frneuC = chdf$neutral_c/chdf$neutral_amount
chdf$frneuT = chdf$neutral_T/chdf$neutral_amount
chdf$Stg = ((chdf$frneuA + chdf$frneuC)-(chdf$frneuG+chdf$frneuT))
Stggraph = ggplot(data = chdf, aes(x = gene_name, y = Stg))+
  geom_boxplot(aes(fill = trophic_level), notch = TRUE)
Stggraph + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
Stggraph1 = ggplot(data = chdf, aes(x = gene_name, y = Stg))+
  geom_boxplot(aes(fill = realm), notch = TRUE)
Stggraph1 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
Ab = ggplot(data = chdf, aes(x = gene_name, y = frneuA))+
  geom_boxplot()
Ab = Ab + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
Gb = ggplot(data = chdf, aes(x = gene_name, y = frneuG))+
  geom_boxplot()
Gb = Gb + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
Cb = ggplot(data = chdf, aes(x = gene_name, y = frneuC))+
  geom_boxplot()
Cb = Cb + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
Tb = ggplot(data = chdf, aes(x = gene_name, y = frneuT))+
  geom_boxplot()
Tb = Tb + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
all = ggarrange(Ab, Gb, Cb, Tb, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)
all



#course work grafs
extradf = brds_clsup
extradf$nA_fraction = extradf$neutral_A/extradf$neutral_amount
extradf$nG_fraction = extradf$neutral_g/extradf$neutral_amount
extradf$nC_fraction = extradf$neutral_c/extradf$neutral_amount
extradf$nT_fraction = extradf$neutral_T/extradf$neutral_amount
g1 = ggplot(data = extradf, aes(x = gene_name, y = nA_fraction))+
  geom_boxplot()+
  xlab('Митохондриальные гены')+
  ylab('Частота Аденина')
g1 = g1 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
g1
g2 = ggplot(data = extradf, aes(x = gene_name, y = nG_fraction))+
  geom_boxplot()+
  xlab('Митохондриальные гены')+
  ylab('Частота Гуанина')
g2 = g2 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
g2
g3 = ggplot(data = extradf, aes(x = gene_name, y = nC_fraction))+
  geom_boxplot()+
  xlab('Митохондриальные гены')+
  ylab('Частота Цитозина')
g3 = g3 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
g3
g4 = ggplot(data = extradf, aes(x = gene_name, y = nT_fraction))+
  geom_boxplot()+
  xlab('Митохондриальные гены')+
  ylab('Частота Тимина')
g4 = g4 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
g4
g5 = ggarrange(g1, g2, g3, g4, 
                labels = c("A", "B", "C","D"),
                ncol = 2, nrow = 2)
g5
ecology = sup = read.table("../../Body/3Results/Birds supplementary materials - DatabaseS1.csv", header = TRUE, sep = ',') 
ecology = ecology[c(1,15,16)]
names(ecology) = c('species_name', 'Realm', "Trophic_level" )
extradf$species_name = gsub(' ','_', extradf$species_name)
extradf = merge(extradf, ecology, by = 'species_name')
extradf$GhAhSkew = (extradf$neutral_c - extradf$neutral_T)/(extradf$neutral_c + extradf$neutral_T)
extradf$Stg_ac = (extradf$nA_fraction + extradf$nC_fraction) - (extradf$nT_fraction + extradf$nG_fraction)
g6 = ggplot(data = extradf, aes(x = gene_name, y = GhAhSkew))+
  geom_boxplot(aes(fill = Realm))+
  xlab('Митохондриальные гены')+
  ylab('Скос аденина и гуанина')
g6 = g6 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
g6
g7 = ggplot(data = extradf, aes(x = gene_name, y = GhAhSkew))+
  geom_boxplot(aes(fill = Trophic_level))+
  xlab('Митохондриальные гены')+
  ylab('Скос аденина и гуанина')
g7 = g7 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
g7
g8 = ggplot(data = extradf, aes(x = gene_name, y = Stg_ac))+
  geom_boxplot(aes(fill = Realm))+
  xlab('Митохондриальные гены')+
  ylab('Stg - Sac')
g8 = g8 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
g8
g9 = ggplot(data = extradf, aes(x = gene_name, y = Stg_ac))+
  geom_boxplot(aes(fill = Trophic_level))+
  xlab('Митохондриальные гены')+
  ylab('Stg - Sac')
g9 = g9 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
g9
