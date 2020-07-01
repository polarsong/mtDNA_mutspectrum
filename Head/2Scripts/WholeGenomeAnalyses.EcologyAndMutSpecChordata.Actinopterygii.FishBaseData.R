###################################
rm(list=ls(all=TRUE))
library(ggpubr)

setwd("../../Body/3Results")
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

# SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE)
names(SynNuc)
SynNuc = SynNuc[SynNuc$Gene != 'ND6',]

SynNuc = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species), FUN = sum)
names(SynNuc) = c('Species','NeutralA','NeutralT','NeutralG','NeutralC')
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 

### merge with temperature
TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)
TEMPE = aggregate(Temperature ~ ., median, data = TEMPE); summary(TEMPE$Temperature)
SynNuc = merge(TEMPE,SynNuc, by = 'Species', all = TRUE); summary(SynNuc$Temperature)

###### merge with fishbase maturation

MATUTM = read.table('../../Body/1Raw/FishBaseMaturity_Tm.txt',  header = TRUE)
MATUTM = aggregate(Tm ~ ., median, data = MATUTM) 
SynNuc = merge(MATUTM,SynNuc, by = 'Species', all = TRUE); nrow(SynNuc)

SynNuc$CtoTSkew = (SynNuc$FrC-SynNuc$FrT)/(SynNuc$FrC+SynNuc$FrT); summary(SynNuc$CtoTSkew) 
SynNuc$GtoASkew = (SynNuc$FrG-SynNuc$FrA)/(SynNuc$FrG+SynNuc$FrA); summary(SynNuc$GtoASkew)

### ANALYSES:
summary(SynNuc$Temperature)
summary(SynNuc$Tm)
summary(lm(FrT ~ scale(Temperature)+scale(Tm), data = SynNuc))
summary(lm(FrT ~ log2(Temperature + 2)*log2(Tm), data = SynNuc))  # keep it for presentation!!!
summary(lm(FrG ~ log2(Temperature + 2)+log2(Tm), data = SynNuc)) # strong
summary(lm(FrA ~ log2(Temperature + 2)+log2(Tm), data = SynNuc)) # strong
summary(lm(GtoASkew ~ log2(Temperature + 2)+log2(Tm), data = SynNuc)) # the highest R^2 = 0.17 
summary(lm(CtoTSkew ~ log2(Temperature + 2)+log2(Tm), data = SynNuc))

#### Rank corr
cor.test(log2(SynNuc$Temperature),SynNuc$FrA, method = "spearman")
cor.test(log2(SynNuc$Temperature),SynNuc$FrT, method = "spearman")
cor.test(log2(SynNuc$Temperature),SynNuc$FrG, method = "spearman")
cor.test(log2(SynNuc$Temperature),SynNuc$FrC, method = "spearman")

pdf("../../Body/4Figures/WholeGenomeAnalyses.NucContent.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE2A.pdf")
ggscatter(SynNuc, x = "Temperature", y = "FrT",
          color = "#0918e6", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n"), xlab="Temperature, C", ylab="FrA")
ggscatter(SynNuc, x = "Temperature", y = "FrG",
          color = "#009414", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n"), xlab="Temperature, C", ylab="FrC")
dev.off()
