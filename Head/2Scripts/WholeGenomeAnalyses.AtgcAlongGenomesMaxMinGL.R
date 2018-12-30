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

############ Generation length
GenLength = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', sep='\t', header=TRUE)
GenLength$Species = gsub(' ','_',GenLength$Scientific_name)
GenLength = GenLength[,c(11,13)]

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

SynNucAll = SynNuc

####### GenLength
SynNucAll = merge(SynNucAll, GenLength, by='Species')

min_gl = SynNucAll[SynNucAll$GenerationLength_d < quantile(SynNucAll$GenerationLength_d, 0.25),]
max_gl = SynNucAll[SynNucAll$GenerationLength_d > quantile(SynNucAll$GenerationLength_d, 0.75),]


Gene = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB', 'ND1','ND2') # ATP6 and ND4 
Timing = seq(1:13)
NewData = data.frame(Gene,Timing)
SynNuc_max = merge(max_gl, NewData)
SynNuc_min = merge(min_gl, NewData)
  # SynNuc = SynNuc[SynNuc$Gene != 'ND6' & SynNuc$Gene != 'ATP8' ,] # !!!!!!!!!!! - in this case it is similar

VecOfSpecies_max  = as.character(unique(SynNuc_max$Species))
VecOfSpecies_min  = as.character(unique(SynNuc_min$Species))

pdf("../../Body/4Figures/WholeGenomeAnalyses.AtgcAlongGenomesMaxMinGL.R.01.pdf", height = 20, width = 40)
#par(mfrow=c(2,1))
# par(cex = 3) # par(mfrow=c(2,2),oma = c(0, 0, 2, 0),cex.main = 2, cex.lab = 2)
par(cex = 2)

###### q1 GL
ColG = rgb(0.1,0.1,0.1,0.1)
ColT = rgb(0.1,0.1,1,0.1)
ColC = rgb(0.1,1,0.1,0.1)
ColA = rgb(1,0.1,0.1,0.1)

plot(NA, xlim=c(1,13), ylim=c(0,0.8), xlab='', ylab="Nucleotide Fractions", main = 'Q1 generation length', xaxt="n")
axis(side = 1, at=c(1:13), labels=c(Gene), las = 2) 

for (i in 1:length(VecOfSpecies_min))
{ # i = 18
  Temp = SynNuc_min[SynNuc_min$Species == VecOfSpecies_min[i],]
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

legend("topright",legend=c('A','T','G','C'), col = c(ColA,ColT,ColG,ColC), pch = 16, horiz = FALSE)

####### q4 GL

plot(NA, xlim=c(1,13), ylim=c(0,0.8), xlab='', ylab="Nucleotide Fractions", main = 'Q4 generation length', xaxt="n")
axis(side = 1, at=c(1:13), labels=c(Gene), las = 2) 

for (i in 1:length(VecOfSpecies_max))
  { # i = 18
    Temp = SynNuc_max[SynNuc_max$Species == VecOfSpecies_max[i],]
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

legend("topright",legend=c('A','T','G','C'), col = c(ColA,ColT,ColG,ColC), pch = 16, horiz = FALSE)


dev.off()

NeutralMinGL = SynNuc_min[, c("Species", "Gene", "NeutralA", "NeutralT", "NeutralG", "NeutralC",
                              "GenerationLength_d")]
NeutralMaxGL = SynNuc_max[, c("Species", "Gene", "NeutralA", "NeutralT", "NeutralG", "NeutralC",
                              "GenerationLength_d")]

summary(NeutralMinGL)
summary(NeutralMaxGL)

write.table(NeutralMaxGL, '../../Body/3Results/ATGCinHighGLspecies.txt', sep='\t',
            row.names = FALSE, quote = FALSE)
write.table(NeutralMinGL, '../../Body/3Results/ATGCinLowGLspecies.txt', sep='\t',
            row.names = FALSE, quote = FALSE)
