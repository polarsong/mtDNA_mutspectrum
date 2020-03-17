###################################

rm(list=ls(all=TRUE))

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip", exdir = "../../Body/3Results/")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")

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

############# Longevity 
AA = read.table("../../Body/1Raw/anage_data.txt", header = TRUE, sep = '\t')
AA$Species = paste(AA$Genus,AA$Species,sep = '_')

pdf("../../Body/4Figures/WholeGenomeAnalyses.EcologyAndMutSpecChordataNoOverlap.ActinopterygiiOnly.R.01.pdf", width = 25, height = 25)
###########
# for (i in 1:length(VecOfTaxaShort))
i = 1
# { # i = 3  = aves; i = 2 -> reptilia; i = 5 -> amphibia; i = 1 => actionopterygii
  VecOfTaxaShort[i]
  SynNuc = SynNucAll
  SynNuc = SynNuc[SynNuc$TAXON == VecOfTaxaShort[i],]

  ############ merge with ecology
  SynNucAA = merge(AA,SynNuc)
  length(unique(SynNucAA$Species))  # 208 species

  ########### question 1: which nucleotides better correlate with MaximalLongevity of fish: T
  AGG = aggregate(list(SynNucAA$FrA,SynNucAA$FrT,SynNucAA$FrG,SynNucAA$FrC), by = list(SynNucAA$Species,SynNucAA$Maximum.longevity..yrs.), FUN = mean)
  names(AGG) = c('Species','FemaleMaturityDays','FrA','FrT','FrG','FrC')
  nrow(AGG) # 206
  
  ############## T - negative, C - positive, A - positive
  cor.test(log2(AGG$FemaleMaturityDays),AGG$FrA)
  cor.test(log2(AGG$FemaleMaturityDays),AGG$FrT)
  cor.test(log2(AGG$FemaleMaturityDays),AGG$FrG)
  cor.test(log2(AGG$FemaleMaturityDays),AGG$FrC)
  
  # don't use G from the very beginning, because - no rank cor: use three other nucleotides:
  A = lm(log2(AGG$FemaleMaturityDays) ~  scale(AGG$FrA) + scale(AGG$FrT) + scale(AGG$FrC));   summary(A)
  A = lm(log2(AGG$FemaleMaturityDays) ~  scale(AGG$FrT) + scale(AGG$FrC));   summary(A)
  A = lm(log2(AGG$FemaleMaturityDays) ~  scale(AGG$FrT)); summary(A)

########### question 2: which nucleotides better correlate with GenerationLength: T
AGG = aggregate(list(SynNucAA$FrA,SynNucAA$FrT,SynNucAA$FrG,SynNucAA$FrC), by = list(SynNucAA$Species,SynNucAA$Female.maturity..days.), FUN = mean)
names(AGG) = c('Species','FemaleMaturityDays','FrA','FrT','FrG','FrC')
nrow(AGG) # 91
 
############## T - negative, C - positive, 
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrA)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrT)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrG)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrC)

# don't use G and A from the very beginning, because - no rank cor: use three two nucleotides:
A = lm(log2(AGG$FemaleMaturityDays) ~  scale(AGG$FrT) + scale(AGG$FrC));   summary(A)
A = lm(log2(AGG$FemaleMaturityDays) ~  scale(AGG$FrT)); summary(A)

########### question 3: which nucleotides better correlate with BodyMass: T and C
AGG = aggregate(list(SynNucAA$FrA,SynNucAA$FrT,SynNucAA$FrG,SynNucAA$FrC), by = list(SynNucAA$Species,SynNucAA$Adult.weight..g.), FUN = mean)
names(AGG) = c('Species','FemaleMaturityDays','FrA','FrT','FrG','FrC')
nrow(AGG) # 126

############## T - negative, C - positive, 
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrA)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrT)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrG)
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrC)

# don't use G and A from the very beginning, because - no rank cor: use three two nucleotides:
A = lm(log2(AGG$FemaleMaturityDays) ~  scale(AGG$FrT) + scale(AGG$FrC));   summary(A)

#}  
dev.off()

