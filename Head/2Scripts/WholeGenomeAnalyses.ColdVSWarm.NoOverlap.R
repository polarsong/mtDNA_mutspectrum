###################################
###### 
###################################

rm(list=ls(all=TRUE))
pdf("../../Body/4Figures/WholeGenomeAnalyses.ColdVSWarm.NoOverlap.R.01.pdf")
par(mfrow=c(2,2))

############ AnAge (for alternative info)
AnAge = read.table("../../Body/1Raw/anage_data.txt", sep = '\t', header=TRUE)
AnAge$Species = paste(AnAge$Genus, AnAge$Species, sep = '_')

table(AnAge$Class)
AnAge$TempInCelcius = AnAge$Temperature..K. - 273.15 # T(°C) = T(K) - 273.15
MammalinBodyTemp = AnAge[AnAge$Class == 'Mammalia',]$TempInCelcius

############ lists of cold and tropical vertebrates
all = read.table("../../Body/1Raw/TemperatureAllChordataDataSet.txt", sep = '\t', header=TRUE)
table(all$Tax)

amph = all[all$Tax == "amphibia",]; amph$Species = gsub(' ','_', amph$Species); amph$Tax= gsub('amphibia','Amphibia', amph$Tax)
rept = all[all$Tax == "reptilia",]; rept$Species = gsub(' ','_', rept$Species); rept$Tax=gsub('r','R', rept$Tax)
fish = all[all$Tax == "fishes",]; fish$Species = gsub(' ','_', fish$Species); fish$Tax= gsub('fishes','Actinopterygii', fish$Tax)
birds = all[all$Tax == "birds",]; fish$Species = gsub(' ','_', fish$Species); fish$Tax= gsub('birds','Aves', fish$Tax)  # 
summary(birds$T..oC.)  # ?????

boxplot(fish$T..oC.,amph$T..oC.,rept$T..oC.,MammalinBodyTemp,birds$T..oC., notch = TRUE, names = c('Actinopterygii','Amphibia','Reptilia','Mammalia','Aves'), outline = FALSE)  # dev.off()

## look for some eco correlations - BMR and body mass are strongly positive, others - nothing strong 
plot(log2(fish$Body.mass..gram.),log2(fish$T..oC.)); cor.test(fish$Body.mass..gram.,fish$T..oC.) # nothing!!
plot(log2(fish$Basal.metabolic.rate..watts.),log2(fish$T..oC.)); cor.test(fish$Basal.metabolic.rate..watts.,fish$T..oC.) # nothing or positive trend!
plot(log2(fish$Basal.metabolic.rate..watts.),log2(fish$Body.mass..gram.)); cor.test(log2(fish$Basal.metabolic.rate..watts.),log2(fish$Body.mass..gram.)) # super positive

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", sep = '\t', header = TRUE) # , fill = TRUE
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

####### compare nucleotide frequencies 
SynNucAMPH = SynNuc[SynNuc$Class == 'Amphibia',]; length(unique(SynNucAMPH$Species))
SynNucREPT = SynNuc[SynNuc$Class == 'Reptilia',]; length(unique(SynNucREPT$Species))
SynNucACTI = SynNuc[SynNuc$Class == 'Actinopterygii',]; length(unique(SynNucACTI$Species))
AGG = aggregate(list(SynNuc$FrA,SynNuc$FrT,SynNuc$FrG,SynNuc$FrC), by = list(SynNuc$Species), FUN = mean)
names(AGG)=c('Species','FrA','FrT','FrG','FrC')

##### merge amph, rept, fish with AGG (neutral syn  info)
amph = merge(amph,AGG, by = 'Species'); # just 10
rept = merge(rept,AGG, by = 'Species'); # just 17 
fish = merge(fish,AGG, by = 'Species'); # 86

######### try cor.test between G and temperature: amph - nothing; rept - nothing; fishes - T and G
cor.test(amph$T..oC.,amph$FrA, method = 'spearman')
cor.test(amph$T..oC.,amph$FrT, method = 'spearman')
cor.test(amph$T..oC.,amph$FrG, method = 'spearman')
cor.test(amph$T..oC.,amph$FrC, method = 'spearman')

cor.test(rept$T..oC.,rept$FrA, method = 'spearman')
cor.test(rept$T..oC.,rept$FrT, method = 'spearman')
cor.test(rept$T..oC.,rept$FrG, method = 'spearman')
cor.test(rept$T..oC.,rept$FrC, method = 'spearman')

cor.test(fish$T..oC.,fish$FrA, method = 'spearman') # a bit positive:                                   # 0.2095367, 0.05283
cor.test(fish$T..oC.,fish$FrT, method = 'spearman') # a bit negative. Warm fishes live longer!?         # -0.2669443; 0.01297
cor.test(fish$T..oC.,fish$FrG, method = 'spearman') # a bit negative. Effect of the temperature alone!? # -0.2170128; 0.04475
cor.test(fish$T..oC.,fish$FrC, method = 'spearman') # nothing

cor.test(fish$Basal.metabolic.rate..watts.,fish$FrA, method = 'spearman') # nothing
cor.test(fish$Basal.metabolic.rate..watts.,fish$FrT, method = 'spearman') # nothing
cor.test(fish$Basal.metabolic.rate..watts.,fish$FrG, method = 'spearman') # nothing
cor.test(fish$Basal.metabolic.rate..watts.,fish$FrC, method = 'spearman') # nothing

c<-lm(fish$T..oC. ~ scale(fish$FrA) + scale(fish$FrT) + scale(fish$FrG)); summary(c); # similar negative coefficients -> either T is also partially driven by temperature, or it is indurectly linked through ecology.
b<-lm(fish$T..oC. ~ scale(fish$FrT) + scale(fish$FrG)); summary(b); # similar negative coefficients -> either T is also partially driven by temperature, or it is indurectly linked through ecology.

#### ALINA, PICs!!!!!!!!!!!!!!
cor.test(fish$T..oC.,fish$FrT, method = 'spearman') # a bit negative. Warm fishes live longer!?         # -0.2669443; 0.01297
cor.test(fish$T..oC.,fish$FrG, method = 'spearman') # a bit negative. Effect of the temperature alone!? # -0.2170128; 0.04475

#### low T in warm fishes is a result of high lifespan (longevity) which is correlated with body size and toC or warm body directly? 
#### when I merge with AnAge (MaxLifespan) I loose half of species (-> 47) and with the rest we see only effect of T - probably because it is more variable than G?
fish = merge(fish,AnAge, by = 'Species')  # 47
cor.test(fish$T..oC.,fish$FrT, method = 'spearman'); plot(fish$T..oC.,fish$FrT) # very negative
cor.test(fish$T..oC.,fish$FrG, method = 'spearman'); plot(fish$T..oC.,fish$FrG)  # nothing
a<-lm(fish$T..oC. ~ fish$FrT + fish$FrG); summary(a) # only T is negatively associated with T..oC.

cor.test(fish$Maximum.longevity..yrs.,fish$FrT, method = 'spearman') # nothing
cor.test(fish$Maximum.longevity..yrs.,fish$FrG, method = 'spearman') # a bit negative
b<-lm(fish$Maximum.longevity..yrs. ~ fish$FrT + fish$FrG); summary(b) # only T is negatively correlated with MaxLifespan
plot(fish$FrG,fish$FrT)

dev.off()

########### Cold+Warm fishes
ColdFishes = fish[fish$T..oC. < quantile(fish$T..oC., 0.400),]; ColdFishesSpecies=ColdFishes$Species
WarmFishes = fish[fish$T..oC. > quantile(fish$T..oC., 0.600),]; WarmFishesSpecies=WarmFishes$Species

pdf("../../Body/4Figures/WholeGenomeAnalyses.coldFishVSwarmFish_alya.NoOverlap.R.01.pdf")

par(mfrow=c(2,2))
boxplot(AGG[AGG$Species %in% ColdFishesSpecies,]$FrA,AGG[AGG$Species %in% WarmFishesSpecies,]$FrA, notch = TRUE, main = 'A', names = c('cold','warm'))
boxplot(AGG[AGG$Species %in% ColdFishesSpecies,]$FrT,AGG[AGG$Species %in% WarmFishesSpecies,]$FrT, notch = TRUE, main = 'T', names = c('cold','warm'))
boxplot(AGG[AGG$Species %in% ColdFishesSpecies,]$FrG,AGG[AGG$Species %in% WarmFishesSpecies,]$FrG, notch = TRUE, main = 'G', names = c('cold','warm'))
boxplot(AGG[AGG$Species %in% ColdFishesSpecies,]$FrC,AGG[AGG$Species %in% WarmFishesSpecies,]$FrC, notch = TRUE, main = 'C', names = c('cold','warm'))

dev.off()

########### Cold+Warm reptilia
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






