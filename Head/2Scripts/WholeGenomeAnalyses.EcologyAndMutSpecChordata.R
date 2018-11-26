###################################

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
SynNucAll = SynNuc
VecOfTaxa = unique(SynNuc$TAXON)
VecOfTaxaShort = c('Actinopterygii','Reptilia','Aves','Mammalia','Amphibia')

############# Longevity 
AA = read.table("../../Body/1Raw/anage_data.txt", header = TRUE, sep = '\t')
AA$Species = paste(AA$Genus,AA$Species,sep = '_')

pdf("../../Body/4Figures/WholeGenomeAnalyses.EcologyAndMutSpecChordata.R.01.pdf", width = 25, height = 25)
###########
for (i in 1:length(VecOfTaxaShort))
{ # i = 1
  VecOfTaxaShort[i]
  SynNuc = SynNucAll
  SynNuc = SynNuc[SynNuc$TAXON == VecOfTaxaShort[i],]

  ############ merge with ecology
  SynNucAA = merge(AA,SynNuc)
  length(unique(SynNucAA$Species))  # 192 species

  ########### which nucleotides better correlate with GT
  AGG = aggregate(list(SynNucAA$FrA,SynNucAA$FrT,SynNucAA$FrG,SynNucAA$FrC), by = list(SynNucAA$Species,SynNucAA$Female.maturity..days.), FUN = mean)
  names(AGG) = c('Species','FemaleMaturityDays','FrA','FrT','FrG','FrC')

  par(mfrow=c(2,2),oma = c(0, 0, 2, 0),cex.main = 2, cex.lab = 2)
  plot(log2(AGG$FemaleMaturityDays),AGG$FrA, col = 'gray', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrA', cex = 2); abline(h =0.4, lt = 2, col = 'red')
  plot(log2(AGG$FemaleMaturityDays),AGG$FrT, col = 'blue', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrT', cex = 2); abline(h =0.4, lt = 2, col = 'red')
  plot(log2(AGG$FemaleMaturityDays),AGG$FrG, col = 'green', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrG', cex = 2); abline(h =0.4, lt = 2, col = 'red')
  plot(log2(AGG$FemaleMaturityDays),AGG$FrC, col = 'cyan', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrC', cex = 2); abline(h =0.4, lt = 2, col = 'red')
  mtext(VecOfTaxaShort[i], outer = TRUE, cex = 1.5)

  ########### which nucleotides better correlate with MLS
  AGG = aggregate(list(SynNucAA$FrA,SynNucAA$FrT,SynNucAA$FrG,SynNucAA$FrC), by = list(SynNucAA$Species,SynNucAA$Maximum.longevity..yrs.), FUN = mean)
  names(AGG) = c('Species','MaxLifeSpan','FrA','FrT','FrG','FrC')
  
  par(mfrow=c(2,2),oma = c(0, 0, 2, 0),cex.main = 2, cex.lab = 2)
  plot(log2(AGG$MaxLifeSpan),AGG$FrA, col = 'gray', ylim = c(0,1), xlab = 'log2(MLS)', main = 'FrA', cex = 2); abline(h =0.4, lt = 2, col = 'red')
  plot(log2(AGG$MaxLifeSpan),AGG$FrT, col = 'blue', ylim = c(0,1), xlab = 'log2(MLS)', main = 'FrT', cex = 2); abline(h =0.4, lt = 2, col = 'red')
  plot(log2(AGG$MaxLifeSpan),AGG$FrG, col = 'green', ylim = c(0,1), xlab = 'log2(MLS)', main = 'FrG', cex = 2); abline(h =0.4, lt = 2, col = 'red')
  plot(log2(AGG$MaxLifeSpan),AGG$FrC, col = 'cyan', ylim = c(0,1), xlab = 'log2(MLS)', main = 'FrC', cex = 2); abline(h =0.4, lt = 2, col = 'red')
  mtext(VecOfTaxaShort[i], outer = TRUE, cex = 1.5)
  
  ########### question 2: which genes better correlate with GT
  ## T is negatively and C is positively (T->C)
  VecOfGenes = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB','ND1','ND2')
  length(VecOfGenes)
  par(mfcol=c(4,13),oma = c(0, 0, 3, 0),cex.main = 2, cex.lab = 2)
  for (j in 1:length(VecOfGenes))
    { # j = 1
    OneGene = SynNucAA[SynNucAA$Gene == VecOfGenes[j],]
    main = VecOfGenes[j]
    plot(log2(OneGene$Female.maturity..days.),OneGene$FrA, col = 'gray', main = main, ylim = c(0,1), xlab = 'log2(MLS)') # a bit negative
    plot(log2(OneGene$Female.maturity..days.),OneGene$FrT, col = 'blue', main = main, ylim = c(0,1), xlab = 'log2(MLS)') # a bit negative
    plot(log2(OneGene$Female.maturity..days.),OneGene$FrG, col = 'green', main = main, ylim = c(0,0.6), xlab = 'log2(MLS)') # a bit negative
    plot(log2(OneGene$Female.maturity..days.),OneGene$FrC, col = 'cyan', main = main, ylim = c(0,0.6), xlab = 'log2(MLS)') # a bit negative
    }
    mtext(VecOfTaxaShort[i], outer = TRUE, cex = 1.5)
}  
dev.off()

