###################################

rm(list=ls(all=TRUE))

setwd("../../Body/3Results")
############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

# SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE)
names(SynNuc)

### make ND6 complementary:
NotND6 = SynNuc[SynNuc$Gene != 'ND6',]
ND6 = SynNuc[SynNuc$Gene == 'ND6',]
A = ND6$NeutralT
T = ND6$NeutralA
G = ND6$NeutralC
C = ND6$NeutralG
ND6$NeutralA = A
ND6$NeutralT = T
ND6$NeutralG = G
ND6$NeutralC = C
SynNuc = rbind(NotND6,ND6)

### count fraction of nucleotides
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 

SynNuc$TAXON = SynNuc$Class
SynNucAll = SynNuc
VecOfTaxa = unique(SynNuc$TAXON)
VecOfTaxaShort = c('Actinopterygii','Reptilia','Aves','Mammalia','Amphibia')

i = 1

  VecOfTaxaShort[i]
  SynNuc = SynNucAll
  SynNuc = SynNuc[SynNuc$TAXON == VecOfTaxaShort[i],]

  
###merge with temperature
TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)
class(TEMPE$Temperature)
class(TEMPE$Species)
TEMPE = aggregate(Temperature ~ ., median, data = TEMPE)
SynNucTEMPE = merge(TEMPE,SynNuc)
class(SynNucTEMPE$Temperature)

AGG = aggregate(list(SynNucTEMPE$FrA,SynNucTEMPE$FrT,SynNucTEMPE$FrG,SynNucTEMPE$FrC), by = list(SynNucTEMPE$Species,SynNucTEMPE$Temperature), FUN = mean)
names(AGG) = c('Species','FemaleMaturityDays','FrA','FrT','FrG','FrC')
nrow(AGG) # 302
AGG = AGG[AGG$FemaleMaturityDays > 0,]

############## T - negative, G - negative
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrA)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrT)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrG)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrC)

######merge with fishbase maturation
MATULM = read.table('../../Body/1Raw/FishBaseMaturity_Lm.txt',  header = TRUE, stringsAsFactors=FALSE)
MATUTM = read.table('../../Body/1Raw/FishBaseMaturity_Tm.txt',  header = TRUE)
class(MATULM$Lm)
MATULM$Lm = as.numeric(MATULM$Lm)
MATUTM = aggregate(Tm ~ ., median, data = MATUTM) 
MATULM = aggregate(Lm ~ ., median, data = MATULM)
SynNucTM = merge(MATUTM,SynNuc)
SynNucLM = merge(MATULM,SynNuc)

AGG = aggregate(list(SynNucTM$FrA,SynNucTM$FrT,SynNucTM$FrG,SynNucTM$FrC), by = list(SynNucTM$Species,SynNucTM$Tm), FUN = mean)
names(AGG) = c('Species','FemaleMaturityDays','FrA','FrT','FrG','FrC')
nrow(AGG) # 188

############## T - negative
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrA)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrT)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrG)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrC)

AGG = aggregate(list(SynNucLM$FrA,SynNucLM$FrT,SynNucLM$FrG,SynNucLM$FrC), by = list(SynNucLM$Species,SynNucLM$Lm), FUN = mean)
names(AGG) = c('Species','FemaleMaturityDays','FrA','FrT','FrG','FrC')
nrow(AGG) # 192

############## T - negative, A - positive
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrA)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrT)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrG)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrC)
