################################
################## COMPARE A T G C nucleotides between taxa / genes: Proper taxa; warm/cold-blooded taxa; all genes, 13 genes separately
################################

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
VecOfTaxa = unique(SynNuc$TAXON)
table(SynNuc$TAXON)/13

### want to use only species with 13 genes (all of them should be already with 13 genes):
SynNuc$Number = 1
AGG = aggregate(SynNuc$Number, by = list(SynNuc$Species), FUN = sum)
VecOfGoodSpecies = AGG[AGG$x == 13,]$Group.1
length(VecOfGoodSpecies)
nrow(SynNuc) # 46670
SynNuc = SynNuc[SynNuc$Species %in% VecOfGoodSpecies,]
nrow(SynNuc) # 46670
table(SynNuc$TAXON)/13
# Actinopterygii       Amphibia    AncientFish           Aves       Mammalia       Reptilia 
# 1770                   205            126              432          788            269

SynNuc$TAXON = ordered(SynNuc$TAXON, levels = c('AncientFish','Actinopterygii','Amphibia','Reptilia','Mammalia','Aves'))
SynNuc$Gene =  ordered(SynNuc$Gene, levels = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB', 'ND1','ND2'))
SynNuc = SynNuc[order(SynNuc$TAXON,SynNuc$Gene),]

#boxplot(FrA ~ TAXON, data = AGG,  notch = TRUE)
pdf("../../Body/4Figures/WholeGenomeAnalyses.AtgcBetweenTaxaNoOverlap.R.01.pdf", height = 20, width = 50)
# G - grey  rgb(0.1,0.1,0.1,0.1)
# A - red   rgb(1,0.1,0.1,0.1)
# T - blue  rgb(0.1,0.1,1,0.1)
# C - green rgb(0.1,1,0.1,0.1)
ColG = rgb(0.1,0.1,0.1,0.2)
ColT = rgb(0.1,0.1,1,0.2)
ColC = rgb(0.1,1,0.1,0.2)
ColA = rgb(1,0.1,0.1,0.2)
TitleActinopterygii = paste('Actinopterygii, N = ',nrow(SynNuc[SynNuc$TAXON == 'Actinopterygii',])/13, sep = '')
TitleAmphibia = paste('Amphibia, N = ',nrow(SynNuc[SynNuc$TAXON == 'Amphibia',])/13, sep = '')
TitleReptilia = paste('Reptilia, N = ',nrow(SynNuc[SynNuc$TAXON == 'Reptilia',])/13, sep = '')
TitleMammalia = paste('Mammalia, N = ',nrow(SynNuc[SynNuc$TAXON == 'Mammalia',])/13, sep = '')
TitleAves = paste('Aves, N = ',nrow(SynNuc[SynNuc$TAXON == 'Aves',])/13, sep = '')

### aggregate nucleotide count per each species and calculate fraction per all 13 protein-coding genes

SYN = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species, SynNuc$Class, SynNuc$TAXON, SynNuc$Taxonomy), FUN = sum)
names(SYN) = c('Species', 'Class', 'TAXON', 'Taxonomy','NeutralA','NeutralT','NeutralG','NeutralC')
SYN$FrA = SYN$NeutralA / (SYN$NeutralA + SYN$NeutralT + SYN$NeutralG + SYN$NeutralC)
SYN$FrT = SYN$NeutralT / (SYN$NeutralA + SYN$NeutralT + SYN$NeutralG + SYN$NeutralC)
SYN$FrG = SYN$NeutralG / (SYN$NeutralA + SYN$NeutralT + SYN$NeutralG + SYN$NeutralC)
SYN$FrC = SYN$NeutralC / (SYN$NeutralA + SYN$NeutralT + SYN$NeutralG + SYN$NeutralC)

par(mfrow=c(1,5))
par(cex.main = 4)
par(cex.axis = 3)
par(oma = c(1, 1, 1, 1)) 
boxplot(SYN[SYN$TAXON == 'Actinopterygii',]$FrG, SYN[SYN$TAXON == 'Actinopterygii',]$FrT, SYN[SYN$TAXON == 'Actinopterygii',]$FrC, SYN[SYN$TAXON == 'Actinopterygii',]$FrA, notch = TRUE, outline = FALSE, las = 1, col = c(ColG, ColT, ColC, ColA), main = TitleActinopterygii, ylim = c(0,0.65), names = c('G','T','C','A')); abline(h=0.25, col = 'red', lt = 2)
boxplot(SYN[SYN$TAXON == 'Amphibia',]$FrG, SYN[SYN$TAXON == 'Amphibia',]$FrT, SYN[SYN$TAXON == 'Amphibia',]$FrC, SYN[SYN$TAXON == 'Amphibia',]$FrA, notch = TRUE, outline = FALSE, las = 1, col = c(ColG, ColT, ColC, ColA), main = TitleAmphibia, ylim = c(0,0.65), names = c('G','T','C','A'));                               abline(h=0.25, col = 'red', lt = 2)
boxplot(SYN[SYN$TAXON == 'Reptilia',]$FrG, SYN[SYN$TAXON == 'Reptilia',]$FrT, SYN[SYN$TAXON == 'Reptilia',]$FrC, SYN[SYN$TAXON == 'Reptilia',]$FrA, notch = TRUE, outline = FALSE, las = 1, col = c(ColG, ColT, ColC, ColA), main = TitleReptilia, ylim = c(0,0.65), names = c('G','T','C','A'));                               abline(h=0.25, col = 'red', lt = 2)
boxplot(SYN[SYN$TAXON == 'Mammalia',]$FrG, SYN[SYN$TAXON == 'Mammalia',]$FrT, SYN[SYN$TAXON == 'Mammalia',]$FrC, SYN[SYN$TAXON == 'Mammalia',]$FrA, notch = TRUE, outline = FALSE, las = 1, col = c(ColG, ColT, ColC, ColA), main = TitleMammalia, ylim = c(0,0.65), names = c('G','T','C','A'));                               abline(h=0.25, col = 'red', lt = 2)
boxplot(SYN[SYN$TAXON == 'Aves',]$FrG, SYN[SYN$TAXON == 'Aves',]$FrT, SYN[SYN$TAXON == 'Aves',]$FrC, SYN[SYN$TAXON == 'Aves',]$FrA, notch = TRUE, outline = FALSE, las = 1, col = c(ColG, ColT, ColC, ColA), main = TitleAves, ylim = c(0,0.65), names = c('G','T','C','A'));                                                   abline(h=0.25, col = 'red', lt = 2)
wilcox.test(SYN[SYN$TAXON == 'Amphibia',]$FrT, SYN[SYN$TAXON == 'Amphibia',]$FrC) # 0.008962 PAPER DEC 2018
wilcox.test(SYN[SYN$TAXON == 'Aves',]$FrC, SYN[SYN$TAXON == 'Aves',]$FrA)         # 5.602e-09 PAPER DEC 2018

par(mfrow=c(1,4))
par(oma = c(14, 3, 0, 0)) 
boxplot(SYN[SYN$TAXON == 'Actinopterygii' | SYN$TAXON == 'Amphibia' | SYN$TAXON == 'Reptilia',]$FrG,SYN[SYN$TAXON == 'Mammalia' | SYN$TAXON == 'Aves',]$FrG, notch = TRUE, outline = FALSE, las = 2, col = c(ColG), main = 'Fraction of G in vertebrate classes', names = c('Cold-blooded','Warm-blooded'), ylim = c(0,0.15)); abline(h=0.25, col = 'red', lt = 2)
boxplot(SYN[SYN$TAXON == 'Actinopterygii',]$FrG,  SYN[SYN$TAXON == 'Amphibia',]$FrG, SYN[SYN$TAXON == 'Reptilia',]$FrG,SYN[SYN$TAXON == 'Mammalia',]$FrG,SYN[SYN$TAXON == 'Aves',]$FrG, notch = TRUE, outline = FALSE, las = 2, col = c(ColG), main = 'Fraction of G in vertebrate classes', names = c('Actinopterygii','Amphibia','Reptilia','Mammalia','Aves'), ylim = c(0,0.15)); abline(h=0.25, col = 'red', lt = 2)
# boxplot(SYN[grepl('Monotremata',SYN$Taxonomy),]$FrG,  SYN[grepl('Eutheria',SYN$Taxonomy),]$FrG, notch = TRUE, outline = FALSE, las = 2, col = c(ColG), main = 'Fraction of G in mammals', ylim = c(0,0.15), names = c('Monotremata','Eutheria')); abline(h=0.25, col = 'red', lt = 2)
# wilcox.test(SYN[grepl('Monotremata',SYN$Taxonomy),]$FrG,SYN[grepl('Eutheria',SYN$Taxonomy),]$FrG, alternative =  'greater')
wilcox.test(SYN[SYN$TAXON == 'Actinopterygii' | SYN$TAXON == 'Amphibia' | SYN$TAXON == 'Reptilia',]$FrG,SYN[SYN$TAXON == 'Mammalia' | SYN$TAXON == 'Aves',]$FrG) # 2.2e-16 PAPER DEC 2018
wilcox.test(SYN[SYN$TAXON == 'Actinopterygii',]$FrG,  SYN[SYN$TAXON == 'Amphibia',]$FrG) # 7.18e-13 PAPER DEC 2018
wilcox.test(SYN[SYN$TAXON == 'Amphibia',]$FrG,  SYN[SYN$TAXON == 'Reptilia',]$FrG)       # 3.845e-05 PAPER DEC 2018
wilcox.test(SYN[SYN$TAXON == 'Reptilia',]$FrG,  SYN[SYN$TAXON == 'Mammalia',]$FrG)       # 6.9e-05 PAPER DEC 2018
wilcox.test(SYN[SYN$TAXON == 'Mammalia',]$FrG,  SYN[SYN$TAXON == 'Aves',]$FrG)           # 0.3084 PAPER DEC 2018

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