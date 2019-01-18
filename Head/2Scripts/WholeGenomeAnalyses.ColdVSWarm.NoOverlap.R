###################################
###### 
###################################

rm(list=ls(all=TRUE))

############ lists of cold and tropical fishes
all = read.table("../../Body/1Raw/TemperatureAllChordataDataSet.txt", sep = '\t', header=TRUE)

amph = all[all$Tax == "amphibia",]; amph$Species = gsub(' ','_', amph$Species); amph$Tax=gsub('a','A', amph$Tax)
rept = all[all$Tax == "reptilia",]; rept$Species = gsub(' ','_', rept$Species); rept$Tax=gsub('r','R', rept$Tax)
fish = all[all$Tax == "fishes",]; fish$Species = gsub(' ','_', fish$Species); fish$Tax= gsub('fishes','Actinopterygii', fish$Tax)

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, fill = TRUE)
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


########## compare nucleotide frequencies 
SynNucAMPH = SynNuc[SynNuc$Class == 'Amphibia',]; length(unique(SynNucAMPH$Species))
SynNucREPT = SynNuc[SynNuc$Class == 'Reptilia',]; length(unique(SynNucREPT$Species))
SynNucACTI = SynNuc[SynNuc$Class == 'Actinopterygii',]; length(unique(SynNucACTI$Species))
AGG = aggregate(list(SynNuc$FrA,SynNuc$FrT,SynNuc$FrG,SynNuc$FrC), by = list(SynNuc$Species), FUN = mean)
names(AGG)=c('Species','FrA','FrT','FrG','FrC')


###########Cold+Warm fishes
ColdFishes = fish[fish$T..oC. < quantile(fish$T..oC., 0.400),]; ColdFishesSpecies=ColdFishes$Species
WarmFishes = fish[fish$T..oC. > quantile(fish$T..oC., 0.600),]; WarmFishesSpecies=WarmFishes$Species

pdf("../../Body/4Figures/WholeGenomeAnalyses.coldFishVSwarmFish_alya.NoOverlap.R.01.pdf")

par(mfrow=c(2,2))
boxplot(AGG[AGG$Species %in% ColdFishesSpecies,]$FrA,AGG[AGG$Species %in% WarmFishesSpecies,]$FrA, notch = TRUE, main = 'A', names = c('cold','warm'))
boxplot(AGG[AGG$Species %in% ColdFishesSpecies,]$FrT,AGG[AGG$Species %in% WarmFishesSpecies,]$FrT, notch = TRUE, main = 'T', names = c('cold','warm'))
boxplot(AGG[AGG$Species %in% ColdFishesSpecies,]$FrG,AGG[AGG$Species %in% WarmFishesSpecies,]$FrG, notch = TRUE, main = 'G', names = c('cold','warm'))
boxplot(AGG[AGG$Species %in% ColdFishesSpecies,]$FrC,AGG[AGG$Species %in% WarmFishesSpecies,]$FrC, notch = TRUE, main = 'C', names = c('cold','warm'))

dev.off()

###########Cold+Warm reptilia
ColdReptilia = rept[rept$T..oC. < quantile(rept$T..oC., 0.400),]; ColdReptiliaSpecies=ColdReptilia$Species
WarmReptilia = rept[rept$T..oC. > quantile(rept$T..oC., 0.600),]; WarmReptiliaSpecies=WarmReptilia$Species

pdf("../../Body/4Figures/WholeGenomeAnalyses.coldReptiliaVSwarmReptilia_alya.NoOverlap.R.01.pdf")

par(mfrow=c(2,2))
boxplot(AGG[AGG$Species %in% ColdReptiliaSpecies,]$FrA,AGG[AGG$Species %in% WarmReptiliaSpecies,]$FrA, notch = TRUE, main = 'A', names = c('cold','warm'))
boxplot(AGG[AGG$Species %in% ColdReptiliaSpecies,]$FrT,AGG[AGG$Species %in% WarmReptiliaSpecies,]$FrT, notch = TRUE, main = 'T', names = c('cold','warm'))
boxplot(AGG[AGG$Species %in% ColdReptiliaSpecies,]$FrG,AGG[AGG$Species %in% WarmReptiliaSpecies,]$FrG, notch = TRUE, main = 'G', names = c('cold','warm'))
boxplot(AGG[AGG$Species %in% ColdReptiliaSpecies,]$FrC,AGG[AGG$Species %in% WarmReptiliaSpecies,]$FrC, notch = TRUE, main = 'C', names = c('cold','warm'))

dev.off()

###########Cold+Warm amphibia
ColdAmphibia = amph[amph$T..oC. < quantile(amph$T..oC., 0.400),]; ColdAmphibiaSpecies=ColdAmphibia$Species
WarmAmphibia = amph[amph$T..oC. > quantile(amph$T..oC., 0.600),]; WarmAmphibiaSpecies=WarmAmphibia$Species

pdf("../../Body/4Figures/WholeGenomeAnalyses.coldAmphibiaVSwarmAmphibia_alya.NoOverlap.R.01.pdf")

par(mfrow=c(2,2))
boxplot(AGG[AGG$Species %in% ColdReptiliaSpecies,]$FrA,AGG[AGG$Species %in% WarmReptiliaSpecies,]$FrA, notch = TRUE, main = 'A', names = c('cold','warm'))
boxplot(AGG[AGG$Species %in% ColdReptiliaSpecies,]$FrT,AGG[AGG$Species %in% WarmReptiliaSpecies,]$FrT, notch = TRUE, main = 'T', names = c('cold','warm'))
boxplot(AGG[AGG$Species %in% ColdReptiliaSpecies,]$FrG,AGG[AGG$Species %in% WarmReptiliaSpecies,]$FrG, notch = TRUE, main = 'G', names = c('cold','warm'))
boxplot(AGG[AGG$Species %in% ColdReptiliaSpecies,]$FrC,AGG[AGG$Species %in% WarmReptiliaSpecies,]$FrC, notch = TRUE, main = 'C', names = c('cold','warm'))

dev.off()





### TO DO
# hibernation vs torpor
# effect of generation time (hib are short or long lived...)
# temp of envir of cold-blooded






