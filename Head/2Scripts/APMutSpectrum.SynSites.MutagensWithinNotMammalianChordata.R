###################################
###### 
###################################

rm(list=ls(all=TRUE))
wd <- getwd()
wd = gsub('Head/2_Scripts','',wd)
setwd(wd)

############# Longevity of birds
AA = read.table("./Body/1Raw/anage_data.txt", header = TRUE, sep = '\t')
AA$Species = paste(AA$Genus,AA$Species,sep = '_')

############ SYN MUT
SynNuc = read.table("./Body/2Derived/GcAtSkewNucContBig.csv", header = TRUE)
SynNuc$FrA = SynNuc$A / SynNuc$SitesNumber
SynNuc$FrT = SynNuc$T / SynNuc$SitesNumber
SynNuc$FrG = SynNuc$G / SynNuc$SitesNumber
SynNuc$FrC = SynNuc$C / SynNuc$SitesNumber

SynNuc = SynNuc[SynNuc$TAXON == 'Aves',]

############ merge with ecology
SynNucAA = merge(AA,SynNuc)
length(unique(SynNucAA$Species))  # 192 species

########### question 1: which nucleotides better correlate with GT: log2(GT) = 11 - 0.29*scale(FrT) + 0.33*scale(FrC) (in line with our mutational spectrum result that T->C correlates with generation time)
AGG = aggregate(list(SynNucAA$FrA,SynNucAA$FrT,SynNucAA$FrG,SynNucAA$FrC), by = list(SynNucAA$Species,SynNucAA$Female.maturity..days.), FUN = mean)
names(AGG) = c('Species','FemaleMaturityDays','FrA','FrT','FrG','FrC')

###### start from pairwise correlations and go to multiple linear model:
## it is opposite s compared to mammals

cor.test(log2(AGG$FemaleMaturityDays),AGG$FrA) # nothing
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrT) # weak positive
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrG) # nothing
cor.test(log2(AGG$FemaleMaturityDays),AGG$FrC) # weak negative

A <- lm(log2(AGG$FemaleMaturityDays) ~ AGG$FrT + AGG$FrC)
summary(A)  # 
A <- lm(log2(AGG$FemaleMaturityDays) ~ scale(AGG$FrT) + scale(AGG$FrC))
summary(A)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     8.97856    0.07987 112.410   <2e-16 ***
#  scale(AGG$FrT)  0.04965    0.10183   0.488   0.6267    
#scale(AGG$FrC) -0.21923    0.10183  -2.153   0.0332 *  
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.9211 on 130 degrees of freedom
#Multiple R-squared:  0.07109,	Adjusted R-squared:  0.05679 
#F-statistic: 4.974 on 2 and 130 DF,  p-value: 0.008287

# plot it:

pdf("./Body/4_Figures/APMutSpectrum.SynSites.MutagensWithinMammals.R.01.pdf", width = 50, height = 30)
par(mfrow=c(2,2),oma = c(0, 0, 2, 0),cex.main = 2, cex.lab = 2)
plot(log2(AGG$FemaleMaturityDays),AGG$FrA, col = 'gray', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrA'); abline(h =0.4, lt = 2, col = 'red')
plot(log2(AGG$FemaleMaturityDays),AGG$FrT, col = 'blue', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrT'); abline(h =0.4, lt = 2, col = 'red')
plot(log2(AGG$FemaleMaturityDays),AGG$FrG, col = 'green', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrG'); abline(h =0.4, lt = 2, col = 'red')
plot(log2(AGG$FemaleMaturityDays),AGG$FrC, col = 'cyan', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrC'); abline(h =0.4, lt = 2, col = 'red')
mtext("log2(GT) = 11 - 0.29*scale(FrT) + 0.33*scale(FrC)", outer = TRUE, cex = 1.5)


########### question 2: which genes better correlate with GT (why T in ATP6,COX3 and ND4 do not correlate with GT and high absolute value - fast replication, no tRNA before them?)
## T is negatively and C is positively (T->C)
VecOfGenes = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB','ND1','ND2') # ND1 ND2 ND6
length(VecOfGenes)
par(mfcol=c(4,13))
for (i in 1:length(VecOfGenes))
{ # i = 1
  OneGene = SynNucAA[SynNucAA$Gene == VecOfGenes[i],]
  main = VecOfGenes[i]
  plot(log2(OneGene$Female.maturity..days.),OneGene$FrA, col = 'gray', main = main, ylim = c(0,1), xlab = 'log2(FM)') # a bit negative
  plot(log2(OneGene$Female.maturity..days.),OneGene$FrT, col = 'blue', main = main, ylim = c(0,1), xlab = 'log2(FM)') # a bit negative
  plot(log2(OneGene$Female.maturity..days.),OneGene$FrG, col = 'green', main = main, ylim = c(0,0.6), xlab = 'log2(FM)') # a bit negative
  plot(log2(OneGene$Female.maturity..days.),OneGene$FrC, col = 'cyan', main = main, ylim = c(0,0.6), xlab = 'log2(FM)') # a bit negative
}
dev.off()

