###################################
###### 
###################################

rm(list=ls(all=TRUE))

############ lists of hibernating and torpor animals
cold = read.table("../../Body/1Raw/cold_water_fishes.txt", sep = '\t', header=FALSE)
trop = read.table("../../Body/1Raw/tropical_water_fishes.txt", sep = '\t', header=FALSE)
ListOfColdSpecies = gsub(' ','_', as.character(cold$V1)); length(ListOfColdSpecies)
ListOfTropicalSpecies = gsub(' ','_', as.character(trop$V1)); length(ListOfTropicalSpecies)

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep='\t')
file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")

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



########## compare nucleotide frequencies between hibernating and DT mammals
SynNuc = SynNuc[SynNuc$Class == 'Actinopterygii',]; length(unique(SynNuc$Species))
AGG = aggregate(list(SynNuc$FrA,SynNuc$FrT,SynNuc$FrG,SynNuc$FrC), by = list(SynNuc$Species), FUN = mean)
names(AGG)=c('Species','FrA','FrT','FrG','FrC')
summary(AGG[AGG$Species %in% ListOfColdSpecies,]$FrA)
summary(AGG[AGG$Species %in% ListOfTropicalSpecies,]$FrA)

pdf("../../Body/4Figures/WholeGenomeAnalyses.coldFishVStropicalFish_alya.NoOverlap.R.01.pdf")

par(mfrow=c(2,2))
boxplot(AGG[AGG$Species %in% ListOfColdSpecies,]$FrA,AGG[AGG$Species %in% ListOfTropicalSpecies,]$FrA, notch = TRUE, main = 'A', names = c('Cold','Tropical'))
boxplot(AGG[AGG$Species %in% ListOfColdSpecies,]$FrT,AGG[AGG$Species %in% ListOfTropicalSpecies,]$FrT, notch = TRUE, main = 'T', names = c('Cold','Tropical'))
boxplot(AGG[AGG$Species %in% ListOfColdSpecies,]$FrG,AGG[AGG$Species %in% ListOfTropicalSpecies,]$FrG, notch = TRUE, main = 'G', names = c('Cold','Tropical'))
boxplot(AGG[AGG$Species %in% ListOfColdSpecies,]$FrG,AGG[AGG$Species %in% ListOfTropicalSpecies,]$FrC, notch = TRUE, main = 'C', names = c('Cold','Tropical'))
dev.off()

### TO DO
# hibernation vs torpor
# effect of generation time (hib are short or long lived...)
# temp of envir of cold-blooded






