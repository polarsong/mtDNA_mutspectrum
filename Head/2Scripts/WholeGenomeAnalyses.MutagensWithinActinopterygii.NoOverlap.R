###################################
######
###################################

rm(list=ls(all=TRUE))

setwd("../../Body/3Results/")
getwd()

############

library(gridExtra) # install.packages("gridExtra")
library(grid)

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")){file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")} 

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
VecOfTaxa = unique(SynNuc$TAXON)

########### TEMPERATURE FISHES
TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)
class(TEMPE$Temperature)
class(TEMPE$Species)
TEMPE = aggregate(Temperature ~ ., mean, data = TEMPE)

############ merge (work only with fishes)
SynNucTEMPE = merge(TEMPE,SynNuc)
length(unique(SynNucTEMPE$Species))  # 333 species
table(SynNucTEMPE$TAXON)

########### question 1: which nucleotides better correlate with GT: 
AGG = aggregate(list(SynNucTEMPE$FrA,SynNucTEMPE$FrT,SynNucTEMPE$FrG,SynNucTEMPE$FrC), by = list(SynNucTEMPE$Species,SynNucTEMPE$Temperature), FUN = mean)
names(AGG) = c('Species','Temperature','FrA','FrT','FrG','FrC')



pdf("../../Body/4Figures/WholeGenomeAnalyses.MutagensWithinFishesNoOverlap.R.01.pdf", width = 50, height = 30)


ColG = rgb(0.1,0.1,0.1,0.2)
ColT = rgb(0.1,0.1,1,0.2)
ColC = rgb(0.1,1,0.1,0.2)
ColA = rgb(1,0.1,0.1,0.2)

par(oma = c(2, 2, 0, 0), cex.main = 2, cex.lab = 1.5, cex = 6, pch =16)
plot(log2(AGG$Temperature),AGG$FrA, col = ColA, ylim = c(0,0.6), xlab = '', ylab = '', main = ''); par(new=TRUE);
plot(log2(AGG$Temperature),AGG$FrT, col = ColT, ylim = c(0,0.6), xlab = '', ylab = '', main = ''); par(new=TRUE);
plot(log2(AGG$Temperature),AGG$FrG, col = ColG, ylim = c(0,0.6), xlab = '', ylab = '', main = ''); par(new=TRUE);
plot(log2(AGG$Temperature),AGG$FrC, col = ColC, ylim = c(0,0.6), xlab = 'log2(Temperature)', ylab = 'Nucleotide Content', main = ''); par(new=TRUE);
abline(h =0.25, lt = 1, col = 'red');
legend("topright",legend=c('A','C','T','G'), col = c(ColA,ColC,ColT,ColG), pch = 16, horiz = FALSE)

par(mfrow=c(2,2), oma = c(0, 0, 2, 0),cex.main = 2, cex.lab = 2, cex = 2, pch = 16)
plot(log2(AGG$Temperature),AGG$FrA, col = ColA, ylim = c(0,0.6), xlab = 'log2(T)', main = 'FrA', cex = 2); abline(h =0.25, lt = 2, col = 'red')
a <-lm(AGG$FrA ~ log2(AGG$Temperature)); summary(a); abline(a, lwd = 2, col = 'red')
plot(log2(AGG$Temperature),AGG$FrT, col = ColT, ylim = c(0,0.6), xlab = 'log2(T)', main = 'FrT', cex = 2); abline(h =0.25, lt = 2, col = 'red')
a <-lm(AGG$FrT ~ log2(AGG$Temperature)); summary(a); abline(a, lwd = 5, col = 'dark blue')
plot(log2(AGG$Temperature),AGG$FrG, col = ColG, ylim = c(0,0.6), xlab = 'log2(T)', main = 'FrG', cex = 2); abline(h =0.25, lt = 2, col = 'red')
a <-lm(AGG$FrG ~ log2(AGG$Temperature)); summary(a); abline(a, lwd = 2, col = 'red')
plot(log2(AGG$Temperature),AGG$FrC, col = ColC, ylim = c(0,0.6), xlab = 'log2(T)', main = 'FrC', cex = 2); abline(h =0.25, lt = 2, col = 'red')
a <-lm(AGG$FrC ~ log2(AGG$Temperature)); summary(a); abline(a, lwd = 5, col = 'dark green')

#mtext("log2(GT) = 11.08294 - 0.12 scale(FrT) + 0.45 scale(FrC)", outer = TRUE, cex = 1.5)


dev.off()

########### question 3: pairwise correlations between nucleotide content - which nucleotides anticorrelate with each other better?
### in this way we can try to move from 4 dimensions to 6 dimensions

pdf("../../Body/4Figures/WholeGenomeAnalyses.MutagensWithinFishesNoOverlap.Tables.R.01.pdf")

Res = c()
AT = cor.test(AGG$FrA,AGG$FrT, method = 'spearman'); Res = c('AT', as.numeric(AT[3]), as.numeric(AT[4]));
AG = cor.test(AGG$FrA,AGG$FrG, method = 'spearman'); Res = rbind(Res, c('AG', as.numeric(AG[3]), as.numeric(AG[4])));
AC = cor.test(AGG$FrA,AGG$FrC, method = 'spearman'); Res = rbind(Res, c('AC', as.numeric(AC[3]), as.numeric(AC[4])));
TC = cor.test(AGG$FrT,AGG$FrC, method = 'spearman'); Res = rbind(Res, c('TC', as.numeric(TC[3]), as.numeric(TC[4])));
TG = cor.test(AGG$FrT,AGG$FrG, method = 'spearman'); Res = rbind(Res, c('TG', as.numeric(TG[3]), as.numeric(TG[4])));
CG = cor.test(AGG$FrC,AGG$FrG, method = 'spearman'); Res = rbind(Res, c('CG', as.numeric(CG[3]), as.numeric(CG[4])));
names(Res) = c('Subst', 'Pvalue','SpearmanRho')

grid.newpage()
grid.table(Res)

dev.off()
