################################
################## COMPARE A T G C nucleotides between taxa / genes: Proper taxa; warm/cold-blooded taxa; all genes, 13 genes separately
################################

rm(list=ls(all=TRUE))

### set working directory as a directory where R script was open
wd <- getwd()
wd = gsub('Head/2Scripts','',wd)
setwd(wd)

############ Syn mut
SynNuc = read.table("./Body/2Derived/GcAtSkewNucContBig.csv", header = TRUE)

NotND6 = SynNuc[SynNuc$Gene != 'ND6',]
ND6 = SynNuc[SynNuc$Gene == 'ND6',]
NotND6$FrA = NotND6$A / NotND6$SitesNumber
NotND6$FrT = NotND6$T / NotND6$SitesNumber
NotND6$FrG = NotND6$G / NotND6$SitesNumber
NotND6$FrC = NotND6$C / NotND6$SitesNumber
ND6$FrA = ND6$T / ND6$SitesNumber
ND6$FrT = ND6$A / ND6$SitesNumber
ND6$FrG = ND6$C / ND6$SitesNumber
ND6$FrC = ND6$G / ND6$SitesNumber
SynNuc = rbind(NotND6,ND6)

VecOfTaxa = unique(SynNuc$TAXON)

Gene = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB', 'ND1','ND2') # ATP6 and ND4 
Timing = seq(1:13)
NewData = data.frame(Gene,Timing)
SynNuc = merge(SynNuc,NewData)

# sd(SynNuc$FrA)/mean(SynNuc$FrA)


CoeffVar <- function(x) {sd(x)/mean(x)} 
AGG = aggregate(list(SynNuc$FrA,SynNuc$FrT,SynNuc$FrG,SynNuc$FrC), by = list(SynNuc$TAXON,SynNuc$Gene, SynNuc$Timing), FUN = CoeffVar)
names(AGG) = c('TAXON','Gene','Timing','CoeffVarA','CoeffVarT','CoeffVarG','CoeffVarC')

pdf("./Body/4Figures/APMutSpectrum.SynSites.AtgcCoeffVarTaxaGenes.R.01.pdf", height = 10, width = 15)
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



