################################
################## COMPARE A T G C nucleotides between taxa / genes: Proper taxa; warm/cold-blooded taxa; all genes, 13 genes separately
################################

rm(list=ls(all=TRUE))

SynNuc = read.table("../../Body/3Results/AllGenesCodonUsage.txt", header = TRUE)############ Syn mut

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

VecOfTaxa = unique(SynNuc$TAXON)

Gene = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB', 'ND1','ND2') # ATP6 and ND4 
Timing = seq(1:13)
NewData = data.frame(Gene,Timing)
SynNuc = merge(SynNuc,NewData)

# sd(SynNuc$FrA)/mean(SynNuc$FrA)


CoeffVar <- function(x) {sd(x)/mean(x)} 
AGG = aggregate(list(SynNuc$FrA,SynNuc$FrT,SynNuc$FrG,SynNuc$FrC), by = list(SynNuc$TAXON,SynNuc$Gene, SynNuc$Timing), FUN = CoeffVar)
names(AGG) = c('TAXON','Gene','Timing','CoeffVarA','CoeffVarT','CoeffVarG','CoeffVarC')

pdf("../../Body/4Figures/WholeGenomeAnalyses.AtgcCoeffVarTaxaGenes.R.01.pdf", height = 10, width = 15)
for (i in 1:length(VecOfTaxa))
{
TAX = VecOfTaxa[i]

plot(NA, xlim=c(1,13), ylim=c(0,1), xlab='', ylab="Coefficient of Variation", main = TAX, xaxt="n")
axis(side = 1, at=c(1:13), labels=c(Gene), las = 2); par(new = TRUE) 
lines(AGG[AGG$TAXON == TAX,]$Timing,AGG[AGG$TAXON == TAX,]$CoeffVarA, xaxt="n", xlim=c(1,13), ylim=c(0,1), xlab = '', ylab = '', pch = 16, col = rgb(1,0.1,0.1,1), lwd = 2); par(new = TRUE) 
lines(AGG[AGG$TAXON == TAX,]$Timing,AGG[AGG$TAXON == TAX,]$CoeffVarT, xaxt="n", xlim=c(1,13), ylim=c(0,1), xlab = '', ylab = '', pch = 16, col = rgb(0.1,1,0.1,1), lwd = 2); par(new = TRUE) 
lines(AGG[AGG$TAXON == TAX,]$Timing,AGG[AGG$TAXON == TAX,]$CoeffVarG, xaxt="n", xlim=c(1,13), ylim=c(0,1), xlab = '', ylab = '', pch = 16, col = rgb(0.1,1,1,1), lwd = 2); par(new = TRUE) 
lines(AGG[AGG$TAXON == TAX,]$Timing,AGG[AGG$TAXON == TAX,]$CoeffVarC, xaxt="n", xlim=c(1,13), ylim=c(0,1), xlab = '', ylab = '', pch = 16, col = rgb(0.1,0.1,0.1,1), lwd = 2)
legend(13,1,legend=c('A','T','G','C'), col = c(rgb(1,0.1,0.1,1),rgb(0.1,1,0.1,1),rgb(0.1,1,1,1),rgb(0.1,0.1,0.1,1)), pch = 16, horiz = FALSE)
}
dev.off()



