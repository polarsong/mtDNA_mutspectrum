rm(list=ls(all=TRUE))
library(ggpubr)
library(caper)
library(geiger)

### reading whole genomes database
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

####### obtaining neutral nucleotide fractions in whole genomes
SynNuc = SynNuc[SynNuc$Gene != 'ND6',]
SynNuc = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species), FUN = sum)
names(SynNuc) = c('Species','NeutralA','NeutralT','NeutralG','NeutralC')
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 

######ecology table from us + Kuptsov 
kuptsovtable = read.table("../../Body/2Derived/EcologyMammalianTable01_KuptsovA_ver2_Full.txt", sep='\t', header=TRUE)

##### merging whole genomes of mammalia with ecology table
kuptsovtable$FrT= NULL
allparameters = merge(kuptsovtable, SynNuc, by="Species")
allparameters$Temper = as.numeric(gsub(",", ".", allparameters$Temperature.C._White2003.2006.other.close.species))
allparameters$GenerationLength_d = as.numeric(gsub(",", ".", allparameters$GenerationLength_d))
allparameters$TG = allparameters$FrT+allparameters$FrG
allparameters$AC = allparameters$FrA+allparameters$FrC
allparameters$TG_ACSkew = (allparameters$TG-allparameters$AC)/(allparameters$TG+allparameters$AC); summary(allparameters$TG_ACSkew)
allparameters$TtoCSkew = (allparameters$FrT-allparameters$FrC)/(allparameters$FrT+allparameters$FrC); summary(allparameters$TtoCSkew)
allparameters$AC_TGSkew = allparameters$TG_ACSkew *-1 
summary(allparameters$Temper)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#30.70   35.90   37.00   36.80   38.08   40.10     425 
nrow(allparameters[!is.na(allparameters$Temper),]) #224
allparametersm=allparameters

##### merging whole genomes of fishes with ecology table
SynNuc$CtoTSkew = (SynNuc$FrC-SynNuc$FrT)/(SynNuc$FrC+SynNuc$FrT); summary(SynNuc$CtoTSkew) 
SynNuc$GtoASkew = (SynNuc$FrG-SynNuc$FrA)/(SynNuc$FrG+SynNuc$FrA); summary(SynNuc$GtoASkew)
SynNuc$AtoGSkew = (SynNuc$FrA-SynNuc$FrG)/(SynNuc$FrA+SynNuc$FrG); summary(SynNuc$AtoGSkew)
SynNuc$TtoCSkew = (SynNuc$FrT-SynNuc$FrC)/(SynNuc$FrT+SynNuc$FrC); summary(SynNuc$TtoCSkew)
SynNuc$TG = SynNuc$FrT+SynNuc$FrG
SynNuc$AC = SynNuc$FrA+SynNuc$FrC
SynNuc$TG_ACSkew = (SynNuc$TG-SynNuc$AC)/(SynNuc$TG+SynNuc$AC); summary(SynNuc$TG_ACSkew)
SynNuc$AC_TGSkew = -(SynNuc$TG-SynNuc$AC)/(SynNuc$TG+SynNuc$AC); summary(SynNuc$AC_TGSkew)
TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)
TEMPE = aggregate(Temperature ~ ., median, data = TEMPE); summary(TEMPE$Temperature)
allparametersf = merge(TEMPE,SynNuc, by = 'Species', all = T); summary(SynNuc$Temperature)
MATUTM = read.table('../../Body/1Raw/FishBaseMaturity_Tm.txt',  header = TRUE)
MATUTM = aggregate(Tm ~ ., median, data = MATUTM) 
allparametersf = merge(allparametersf, MATUTM, by = 'Species', all = T); nrow(SynNuc)
allparametersf = allparametersf[!is.na(allparametersf$AC_TGSkew),]
str(allparametersf$Temperature)
str(allparametersf$Tm)
allparametersf$Longevity = round(allparametersf$Tm * 365.25)


######One integral table
allparametersf$LongevityD = 0
allparametersf = allparametersf[!is.na(allparametersf$Longevity),]
allparametersf[allparametersf$Longevity > median(allparametersf$Longevity),]$LongevityD = 1
onedatasetF = data.frame(allparametersf$Species, allparametersf$Temperature, allparametersf$Longevity, allparametersf$LongevityD, allparametersf$AC_TGSkew)
names(onedatasetF)=c("Species", "Temperature", "Longevity", "LongevityDummy", "AC_TGSkew")
onedatasetF$Ectothermy = 1

allparametersm$LongevityD = 0
str(allparametersm$GenerationLength_d)
allparametersm[allparametersm$GenerationLength_d > median(allparametersm$GenerationLength_d),]$LongevityD = 1
onedatasetM = data.frame(allparametersm$Species, allparametersm$Temper, allparametersm$GenerationLength_d, allparametersm$LongevityD, allparametersm$AC_TGSkew)
names(onedatasetM)=c("Species", "Temperature", "Longevity", "LongevityDummy", "AC_TGSkew")
onedatasetM$Ectothermy = 0

onedata = rbind(onedatasetF, onedatasetM)
str(onedata)
summary(lm(formula = AC_TGSkew ~ Temperature + Longevity + Ectothermy, data = onedata))
summary(lm(formula = AC_TGSkew ~ Temperature + LongevityDummy + Ectothermy, data = onedata))
nrow(onedata[!is.na(onedata$Temperature) & !is.na(onedata$AC_TGSkew),])

summary(onedata$Temperature)
summary(onedata$Longevity)
onedata[onedata$Longevity == 18980,]

#write.table(onedata, file="../../Body/2Derived/ForKGAllClasses.txt", row.names = F)
