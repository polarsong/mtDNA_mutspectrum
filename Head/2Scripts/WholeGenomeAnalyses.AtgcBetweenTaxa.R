################################
################## COMPARE A T G C nucleotides between taxa / genes: Proper taxa; warm/cold-blooded taxa; all genes, 13 genes separately
################################

rm(list=ls(all=TRUE))

############ Syn mut
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsage.txt", header = TRUE)

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

SynNuc$TAXON = SynNuc$Class
VecOfTaxa = unique(SynNuc$TAXON)
SynNucAll = SynNuc

str(SynNuc)
SynNuc$TAXON = ordered(SynNuc$TAXON, levels = c('AncientFish','Actinopterygii','Amphibia','Reptilia','Mammalia','Aves'))
SynNuc$Gene =  ordered(SynNuc$Gene, levels = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB', 'ND1','ND2'))
SynNuc = SynNuc[order(SynNuc$TAXON,SynNuc$Gene),]
str(SynNuc)

#boxplot(FrA ~ TAXON, data = AGG,  notch = TRUE)
pdf("../../Body/4Figures/WholeGenomeAnalyses.AtgcBetweenTaxa.R.01.pdf", height = 20, width = 50)
par(mfrow=c(2,2))
par(cex.main = 4)
par(cex.axis = 2)
boxplot(FrA ~ TAXON, data = SynNuc,  notch = TRUE, outline = FALSE, las = 1, col = c('dark blue','blue','orange','green','red','brown'), main = 'FrA'); abline(h=0.2, col = 'red', lt = 2)
boxplot(FrT ~ TAXON, data = SynNuc,  notch = TRUE, outline = FALSE, las = 1, col = c('dark blue','blue','orange','green','red','brown'), main = 'FrT'); abline(h=0.2, col = 'red', lt = 2)
boxplot(FrG ~ TAXON, data = SynNuc,  notch = TRUE, outline = FALSE, las = 1, col = c('dark blue','blue','orange','green','red','brown'), main = 'FrG'); abline(h=0.2, col = 'red', lt = 2)
boxplot(FrC ~ TAXON, data = SynNuc,  notch = TRUE, outline = FALSE, las = 1, col = c('dark blue','blue','orange','green','red','brown'), main = 'FrC'); abline(h=0.2, col = 'red', lt = 2)

par(mfrow=c(1,1))
par(oma = c(12, 0, 0, 0)) # increae space for bottom lines

boxplot(FrA ~ Gene*TAXON, data = SynNuc,  notch = TRUE, outline = FALSE, las = 2, col = c(rep('dark blue',13),rep('blue',13),rep('orange',13),rep('green',13),rep('red',13),rep('brown',13)), main = 'FrA'); abline(h=0.2, col = 'red', lt = 2)
boxplot(FrT ~ Gene*TAXON, data = SynNuc,  notch = TRUE, outline = FALSE, las = 2, col = c(rep('dark blue',13),rep('blue',13),rep('orange',13),rep('green',13),rep('red',13),rep('brown',13)), main = 'FrT'); abline(h=0.2, col = 'red', lt = 2)
boxplot(FrG ~ Gene*TAXON, data = SynNuc,  notch = TRUE, outline = FALSE, las = 2, col = c(rep('dark blue',13),rep('blue',13),rep('orange',13),rep('green',13),rep('red',13),rep('brown',13)), main = 'FrG'); abline(h=0.2, col = 'red', lt = 2)
boxplot(FrC ~ Gene*TAXON, data = SynNuc,  notch = TRUE, outline = FALSE, las = 2, col = c(rep('dark blue',13),rep('blue',13),rep('orange',13),rep('green',13),rep('red',13),rep('brown',13)), main = 'FrC'); abline(h=0.2, col = 'red', lt = 2)

boxplot(FrA ~ TAXON*Gene, data = SynNuc,  notch = TRUE, outline = FALSE, las = 2, col = c(rep(c('dark blue','blue','orange','green','red','brown'),6)), main = 'FrA'); abline(h=0.2, col = 'red', lt = 2)
boxplot(FrT ~ TAXON*Gene, data = SynNuc,  notch = TRUE, outline = FALSE, las = 2, col = c(rep(c('dark blue','blue','orange','green','red','brown'),6)), main = 'FrT'); abline(h=0.2, col = 'red', lt = 2)
boxplot(FrG ~ TAXON*Gene, data = SynNuc,  notch = TRUE, outline = FALSE, las = 2, col = c(rep(c('dark blue','blue','orange','green','red','brown'),6)), main = 'FrG'); abline(h=0.2, col = 'red', lt = 2)
boxplot(FrC ~ TAXON*Gene, data = SynNuc,  notch = TRUE, outline = FALSE, las = 2, col = c(rep(c('dark blue','blue','orange','green','red','brown'),6)), main = 'FrC'); abline(h=0.2, col = 'red', lt = 2)

dev.off()

## 

