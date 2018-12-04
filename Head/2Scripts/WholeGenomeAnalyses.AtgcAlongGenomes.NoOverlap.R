########################
################ ATGC FRACTIONS ALONG GENOMES
########################
## plot segments of ATGC along genomes
## for ND6 we count complementary nucleotides (T instead of A and so on) so that we analyze only light chain content.
## find species with unsusual nucleotide content:

rm(list=ls(all=TRUE))

############ Syn mut
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE)

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

pdf("../../Body/4Figures/WholeGenomeAnalyses.AtgcAlongGenomesNoOverlap.R.01.pdf", height = 10, width = 15)
par(mfrow=c(2,1))


for (taxa in 1:length(VecOfTaxa))
{ # taxa = 1
  ColG = rgb(0.1,0.1,0.1,0.1)
  ColT = rgb(0.1,0.1,1,0.1)
  ColC = rgb(0.1,1,0.1,0.1)
  ColA = rgb(1,0.1,0.1,0.1)
  
  TAX = as.character(VecOfTaxa[taxa])
  SynNuc = SynNucAll
    
  Gene = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB', 'ND1','ND2') # ATP6 and ND4 
  Timing = seq(1:13)
  NewData = data.frame(Gene,Timing)
  SynNuc = merge(SynNuc,NewData)
  
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
  legend(13,0.6,legend=c('A','T','G','C'), col = c(ColA,ColT,ColG,ColC), pch = 16, horiz = FALSE)
}

dev.off()

