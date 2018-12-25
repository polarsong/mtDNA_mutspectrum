###################################
###### 
###################################

rm(list=ls(all=TRUE))

############ lists of cold and tropical fishes
cold = read.table("../../Body/1Raw/cold_water_fishes.txt", sep = '\t',)
tropical = read.table("../../Body/1Raw/tropical_water_fishes.txt", sep = '\t')
cold = gsub(' ','_', cold$V1); length(cold)
tropical = gsub(' ','_', tropical$V1); length(tropical)

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE)
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


########## compare nucleotide frequencies between cold and tropical fishes
SynNuc = SynNuc[SynNuc$Class == 'Actinopterygii',]; length(unique(SynNuc$Species))
AGG = aggregate(list(SynNuc$FrA,SynNuc$FrT,SynNuc$FrG,SynNuc$FrC), by = list(SynNuc$Species), FUN = mean)
names(AGG)=c('Species','FrA','FrT','FrG','FrC')
summary(AGG[AGG$Species %in% cold,]$FrA)
summary(AGG[AGG$Species %in% tropical,]$FrA)

pdf("../../Body/4Figures/WholeGenomeAnalyses.coldVStropical.NoOverlap.R.01.pdf")

par(mfrow=c(2,2))
boxplot(AGG[AGG$Species %in% cold,]$FrA,AGG[AGG$Species %in% tropical,]$FrA, notch = TRUE, main = 'A', names = c('cold','tropical'))
boxplot(AGG[AGG$Species %in% cold,]$FrT,AGG[AGG$Species %in% tropical,]$FrT, notch = TRUE, main = 'T', names = c('cold','tropical'))
boxplot(AGG[AGG$Species %in% cold,]$FrG,AGG[AGG$Species %in% tropical,]$FrG, notch = TRUE, main = 'G', names = c('cold','tropical'))
boxplot(AGG[AGG$Species %in% cold,]$FrC,AGG[AGG$Species %in% tropical,]$FrC, notch = TRUE, main = 'C', names = c('cold','tropical'))

dev.off()

### TO DO
# hibernation vs torpor
# effect of generation time (hib are short or long lived...)
# temp of envir of cold-blooded






