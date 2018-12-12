########################
################ ATGC FRACTIONS ALONG GENOMES
########################
## plot segments of ATGC along genomes
## for ND6 we count complementary nucleotides (T instead of A and so on) so that we analyze only light chain content.
## find species with unsusual nucleotide content:

rm(list=ls(all=TRUE))

############ Syn mut

unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip", exdir = "../../Body/3Results/")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")

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
VecOfTaxa = c('Mammalia', 'Actinopterygii', 'Aves', 'Reptilia', 'Amphibia')
SynNucAll = SynNuc

pdf("../../Body/4Figures/WholeGenomeAnalyses.AtgcAlongGenomesNoOverlap.R.01.pdf", height = 20, width = 40)
#par(mfrow=c(2,1))
par(cex = 3) # par(mfrow=c(2,2),oma = c(0, 0, 2, 0),cex.main = 2, cex.lab = 2)

for (taxa in 1:length(VecOfTaxa))
{ # taxa = 1
  TAX = as.character(VecOfTaxa[taxa])
  SynNuc = SynNucAll
  #ColG = rainbow(4, alpha = 0.1)[4]
  #ColT = rainbow(4, alpha = 0.1)[2]
  #ColC = rainbow(4, alpha = 0.1)[1]
  #ColA = rainbow(4, alpha = 0.1)[3]
  
  ColG = rgb(0.1,0.1,0.1,0.1)
  ColT = rgb(0.1,0.1,1,0.1)
  ColC = rgb(0.1,1,0.1,0.1)
  ColA = rgb(1,0.1,0.1,0.1)

  Gene = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB', 'ND1','ND2') # ATP6 and ND4 
  Timing = seq(1:13)
  NewData = data.frame(Gene,Timing)
  SynNuc = merge(SynNuc,NewData)
  # SynNuc = SynNuc[SynNuc$Gene != 'ND6' & SynNuc$Gene != 'ATP8' ,] # !!!!!!!!!!! - in this case it is similar
  
  AGG = aggregate(list(SynNuc$FrA,SynNuc$FrT,SynNuc$FrG,SynNuc$FrC), by = list(SynNuc$Class,SynNuc$Timing), FUN = mean)
  names(AGG)=c('Class','Timing','FrA','FrT','FrG','FrC')
  
  cor.test(AGG[AGG$Class == 'Actinopterygii',]$FrA,AGG[AGG$Class == 'Actinopterygii',]$Timing, method = 'spearman') # -0.5769231, p = 0.04
  cor.test(AGG[AGG$Class == 'Actinopterygii',]$FrT,AGG[AGG$Class == 'Actinopterygii',]$Timing, method = 'spearman') # -0.2142857, p = 0.4819
  cor.test(AGG[AGG$Class == 'Actinopterygii',]$FrG,AGG[AGG$Class == 'Actinopterygii',]$Timing, method = 'spearman') #  0.2802198, p = 0.3534
  cor.test(AGG[AGG$Class == 'Actinopterygii',]$FrC,AGG[AGG$Class == 'Actinopterygii',]$Timing, method = 'spearman') #  0.7802198, p = 0.002621
  
  cor.test(AGG[AGG$Class == 'Mammalia',]$FrA,AGG[AGG$Class == 'Mammalia',]$Timing, method = 'spearman') # -0.1483516, p = 0.6298
  cor.test(AGG[AGG$Class == 'Mammalia',]$FrT,AGG[AGG$Class == 'Mammalia',]$Timing, method = 'spearman') # -0.3736264, p = 0.2094
  cor.test(AGG[AGG$Class == 'Mammalia',]$FrG,AGG[AGG$Class == 'Mammalia',]$Timing, method = 'spearman') # -0.5714286, p = 0.04489
  cor.test(AGG[AGG$Class == 'Mammalia',]$FrC,AGG[AGG$Class == 'Mammalia',]$Timing, method = 'spearman') #  0.6153846, p = 0.02854
  
  cor.test(AGG[AGG$Class == 'Aves',]$FrA,AGG[AGG$Class == 'Aves',]$Timing, method = 'spearman') # -0.6208791, p = 0.02687
  cor.test(AGG[AGG$Class == 'Aves',]$FrT,AGG[AGG$Class == 'Aves',]$Timing, method = 'spearman') # 0.01648352, p = 0.9639
  cor.test(AGG[AGG$Class == 'Aves',]$FrG,AGG[AGG$Class == 'Aves',]$Timing, method = 'spearman') # 0.1538462, p = 0.6168
  cor.test(AGG[AGG$Class == 'Aves',]$FrC,AGG[AGG$Class == 'Aves',]$Timing, method = 'spearman') #  0.6318681, p = 0.02374
  
  cor.test(AGG[AGG$Class == 'Amphibia',]$FrA,AGG[AGG$Class == 'Amphibia',]$Timing, method = 'spearman') # -0.2087912, p = 0.4935
  cor.test(AGG[AGG$Class == 'Amphibia',]$FrT,AGG[AGG$Class == 'Amphibia',]$Timing, method = 'spearman') # -0.1703297, p = 0.5785
  cor.test(AGG[AGG$Class == 'Amphibia',]$FrG,AGG[AGG$Class == 'Amphibia',]$Timing, method = 'spearman') # 0.08241758, p = 0.7925
  cor.test(AGG[AGG$Class == 'Amphibia',]$FrC,AGG[AGG$Class == 'Amphibia',]$Timing, method = 'spearman') #  0.7307692, p = 0.006323
  
  cor.test(AGG[AGG$Class == 'Reptilia',]$FrA,AGG[AGG$Class == 'Reptilia',]$Timing, method = 'spearman') # -0.3241758 , p = 0.2799
  cor.test(AGG[AGG$Class == 'Reptilia',]$FrT,AGG[AGG$Class == 'Reptilia',]$Timing, method = 'spearman') # -0.08791209, p = 0.7785
  cor.test(AGG[AGG$Class == 'Reptilia',]$FrG,AGG[AGG$Class == 'Reptilia',]$Timing, method = 'spearman') # -0.3571429, p = 0.2315
  cor.test(AGG[AGG$Class == 'Reptilia',]$FrC,AGG[AGG$Class == 'Reptilia',]$Timing, method = 'spearman') #  0.6593407, p = 0.01713
  
  VecOfSpecies  = as.character(unique(SynNuc[SynNuc$TAXON == TAX,]$Species))
  
  plot(NA, xlim=c(1,13), ylim=c(0,0.7), xlab='', ylab="Nucleotide Fractions", main = TAX, xaxt="n")
  axis(side = 1, at=c(1:13), labels=c(Gene), las = 2) 
    for (i in 1:length(VecOfSpecies))
    { # i = 18
      Temp = SynNuc[SynNuc$Species == VecOfSpecies[i],]
      Temp = Temp[order(Temp$Timing),]
      if (nrow(Temp) == 13) # 10
      {
        for (count in 1:(nrow(Temp)-1))
        {
          segments(count, Temp$FrA[count], count+1, Temp$FrA[count+1], col = ColA, lwd = 1) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
          segments(count, Temp$FrT[count], count+1, Temp$FrT[count+1], col = ColT, lwd = 1) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
          segments(count, Temp$FrG[count], count+1, Temp$FrG[count+1], col = ColG, lwd = 1) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
          segments(count, Temp$FrC[count], count+1, Temp$FrC[count+1], col = ColC, lwd = 1) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
        }
      }
    }
  
  ColG = rgb(0.1,0.1,0.1,0.5)
  ColT = rgb(0.1,0.1,1,0.5)
  ColC = rgb(0.1,1,0.1,0.5)
  ColA = rgb(1,0.1,0.1,0.5)
  
  legend("topright",legend=c('A','T','G','C'), col = c(ColA,ColT,ColG,ColC), pch = 16, horiz = FALSE)
}

dev.off()

