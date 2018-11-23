###################################
###### 
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
VecOfTaxa = unique(SynNuc$TAXON)

############# GENERATION LENGTH FOR ALL MAMMALS

GT = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

############ merge (work only with mammals since GT is only for mammals)
SynNucGT = merge(GT,SynNuc)
length(unique(SynNucGT$Species))  # 705 species
table(SynNucGT$TAXON)

########### question 1: which nucleotides better correlate with GT: log2(GT) = 11 - 0.29*scale(FrT) + 0.33*scale(FrC) (in line with our mutational spectrum result that T->C correlates with generation time)
AGG = aggregate(list(SynNucGT$FrA,SynNucGT$FrT,SynNucGT$FrG,SynNucGT$FrC), by = list(SynNucGT$Species,SynNucGT$GenerationLength_d), FUN = mean)
names(AGG) = c('Species','GenerationLength_d','FrA','FrT','FrG','FrC')

###### start from pairwise correlations and go to multiple linear model:

cor.test(log2(AGG$GenerationLength_d),AGG$FrA) # negative
cor.test(log2(AGG$GenerationLength_d),AGG$FrT) # negative
cor.test(log2(AGG$GenerationLength_d),AGG$FrG) # positive (the less strong)
cor.test(log2(AGG$GenerationLength_d),AGG$FrC) # positive

A <- lm(log2(AGG$GenerationLength_d) ~ AGG$FrA + AGG$FrT + AGG$FrC); summary(A)
A <- lm(log2(AGG$GenerationLength_d) ~ AGG$FrT + AGG$FrC); summary(A)
A <- lm(log2(AGG$GenerationLength_d) ~ scale(AGG$FrT) + scale(AGG$FrC)); summary(A) # log2(GT) = 11.06737 -0.11 scale(FrT)+0.46(FrC)


# plot it:

pdf("../../Body/4Figures/WholeGenomeAnalyses.MutagensWithinMammals.R.01.pdf", width = 50, height = 30)
par(mfrow=c(2,2),oma = c(0, 0, 2, 0),cex.main = 2, cex.lab = 2)
plot(log2(AGG$GenerationLength_d),AGG$FrA, col = 'gray', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrA', cex = 2); abline(h =0.25, lt = 2, col = 'red')
plot(log2(AGG$GenerationLength_d),AGG$FrT, col = 'blue', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrT', cex = 2); abline(h =0.25, lt = 2, col = 'red')
plot(log2(AGG$GenerationLength_d),AGG$FrG, col = 'green', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrG', cex = 2); abline(h =0.25, lt = 2, col = 'red')
plot(log2(AGG$GenerationLength_d),AGG$FrC, col = 'cyan', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrC', cex = 2); abline(h =0.25, lt = 2, col = 'red')
mtext("log2(GT) = 11.06737 -0.11*scale(FrT)+0.46*(FrC)", outer = TRUE, cex = 1.5)


########### question 2: which genes better correlate with GT (why T in ATP6,COX3 and ND4 do not correlate with GT and high absolute value - fast replication, no tRNA before them?)
## T is negatively and C is positively (T->C)
VecOfGenes = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB','ND1','ND2') # ND1 ND2 ND6
length(VecOfGenes)
par(mfcol=c(4,13))
for (i in 1:length(VecOfGenes))
{ # i = 1
  OneGene = SynNucGT[SynNucGT$Gene == VecOfGenes[i],]
  main = VecOfGenes[i]
  plot(log2(OneGene$GenerationLength_d),OneGene$FrA, col = 'gray', main = main, ylim = c(0,1), xlab = 'log2(GT)') # a bit negative
  plot(log2(OneGene$GenerationLength_d),OneGene$FrT, col = 'blue', main = main, ylim = c(0,1), xlab = 'log2(GT)') # a bit negative
  plot(log2(OneGene$GenerationLength_d),OneGene$FrG, col = 'green', main = main, ylim = c(0,0.6), xlab = 'log2(GT)') # a bit negative
  plot(log2(OneGene$GenerationLength_d),OneGene$FrC, col = 'cyan', main = main, ylim = c(0,0.6), xlab = 'log2(GT)') # a bit negative
}
dev.off()

