###################################
###### 
###################################

rm(list=ls(all=TRUE))

############ list of hibernating animals
Hib = read.table("../../Body/1Raw/HibernatingMammals.txt", sep = '\t')
ListOfHibSpecies = gsub(' ','_',Hib$V1); length(ListOfHibSpecies)

############ Syn mut
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE)

### make ND6 complementary:
NotND6 = SynNuc[SynNuc$Gene != 'ND6',]
NotND6$FrA = NotND6$NeutralA / (NotND6$NeutralA + NotND6$NeutralT + NotND6$NeutralG + NotND6$NeutralC)
NotND6$FrT = NotND6$NeutralT / (NotND6$NeutralA + NotND6$NeutralT + NotND6$NeutralG + NotND6$NeutralC) 
NotND6$FrG = NotND6$NeutralG / (NotND6$NeutralA + NotND6$NeutralT + NotND6$NeutralG + NotND6$NeutralC) 
NotND6$FrC = NotND6$NeutralC / (NotND6$NeutralA + NotND6$NeutralT + NotND6$NeutralG + NotND6$NeutralC) 

ND6 = SynNuc[SynNuc$Gene == 'ND6',]
ND6$FrA = ND6$NeutralT / (ND6$NeutralA + ND6$NeutralT + ND6$NeutralG + ND6$NeutralC)
ND6$FrT = ND6$NeutralA / (ND6$NeutralA + ND6$NeutralT + ND6$NeutralG + ND6$NeutralC) 
ND6$FrG = ND6$NeutralC / (ND6$NeutralA + ND6$NeutralT + ND6$NeutralG + ND6$NeutralC) 
ND6$FrC = ND6$NeutralG / (ND6$NeutralA + ND6$NeutralT + ND6$NeutralG + ND6$NeutralC) 

SynNuc = rbind(NotND6,ND6)

VecOfTaxa = unique(SynNuc$Class)

############# Hibernation in different classes
#### who is in our list? only mammals!
for (i in 1:length(VecOfTaxa))
{ # i = 1
  Species = unique(SynNuc[SynNuc$Class == VecOfTaxa[i],]$Species); length(Species)
  L = length(intersect(Species,ListOfHibSpecies))
  print(c(as.character(VecOfTaxa[i]), L))
}

########## compare nucleotide frequencies between hibernating and other mammals
SynNuc = SynNuc[SynNuc$Class == 'Mammalia',]; length(unique(SynNuc$Species))
AGG = aggregate(list(SynNuc$FrA,SynNuc$FrT,SynNuc$FrG,SynNuc$FrC), by = list(SynNuc$Species), FUN = mean)
names(AGG)=c('Species','FrA','FrT','FrG','FrC')
summary(AGG[AGG$Species %in% ListOfHibSpecies,]$FrA)
summary(AGG[!AGG$Species %in% ListOfHibSpecies,]$FrA)

pdf("../../Body/4Figures/WholeGenomeAnalyses.HibernationWithinMammals.NoOverlap.R.01.pdf")

par(mfrow=c(2,2))
boxplot(AGG[AGG$Species %in% ListOfHibSpecies,]$FrA,AGG[!AGG$Species %in% ListOfHibSpecies,]$FrA, notch = TRUE, main = 'A', names = c('Hib','Others'))
boxplot(AGG[AGG$Species %in% ListOfHibSpecies,]$FrT,AGG[!AGG$Species %in% ListOfHibSpecies,]$FrT, notch = TRUE, main = 'T', names = c('Hib','Others'))
boxplot(AGG[AGG$Species %in% ListOfHibSpecies,]$FrG,AGG[!AGG$Species %in% ListOfHibSpecies,]$FrG, notch = TRUE, main = 'G', names = c('Hib','Others'))
boxplot(AGG[AGG$Species %in% ListOfHibSpecies,]$FrC,AGG[!AGG$Species %in% ListOfHibSpecies,]$FrC, notch = TRUE, main = 'C', names = c('Hib','Others'))

dev.off()

### TO DO
# hibernation vs torpor
# effect of generation time (hib are short or long lived...)
# temp of envir of cold-blooded






