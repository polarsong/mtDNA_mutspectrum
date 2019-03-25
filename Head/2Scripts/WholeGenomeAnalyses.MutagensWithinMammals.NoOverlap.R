###################################
###### 
###################################

rm(list=ls(all=TRUE))

############

library(gridExtra) # install.packages("gridExtra")
library(grid)

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip", exdir = "../../Body/3Results/")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")

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
length(unique(SynNucGT$Species))  # 649 species
table(SynNucGT$TAXON)

############ data for PICs
library(ape)

tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")
data = SynNucGT[which(as.character(SynNucGT$Species) %in% tree$tip.label),]

df_vec <- as.character(SynNucGT$Species)
tree_vec <- tree$tip.label

a <- setdiff(df_vec, tree_vec)
b <- setdiff(tree_vec, df_vec)

# some errors below:
# row.names(data) = data$Species # it gives an error!!
# tree2 <- drop.tip(tree, b)
# TempData = data[, -c(1, 5)]
# contrasts <- as.data.frame(apply(TempData, 2, pic, tree2))
# names(contrasts) = names(TempData)

########### question 1: which nucleotides better correlate with GT: log2(GT) = 11 - 0.29*scale(FrT) + 0.33*scale(FrC) (in line with our mutational spectrum result that T->C correlates with generation time)
AGG = aggregate(list(SynNucGT$FrA,SynNucGT$FrT,SynNucGT$FrG,SynNucGT$FrC), by = list(SynNucGT$Species,SynNucGT$GenerationLength_d), FUN = mean)
names(AGG) = c('Species','GenerationLength_d','FrA','FrT','FrG','FrC')

############ data for PICs
library(ape)
tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")
data = AGG[which(as.character(AGG$Species) %in% tree$tip.label),]
data$GenerationLength_d = log2(data$GenerationLength_d)  # I added log of generation length
df_vec <- as.character(AGG$Species)
tree_vec <- tree$tip.label
a <- setdiff(df_vec, tree_vec)
b <- setdiff(tree_vec, df_vec)
data = data[-597,] # is it duplicate!!!???
row.names(data) = data$Species
tree2 <- drop.tip(tree, b)
TempData = data[, -1]
contrasts <- as.data.frame(apply(TempData, 2, pic, tree2))
names(contrasts) = names(TempData)

cor.test(contrasts$GenerationLength_d, contrasts$FrT, method= 'spearman')  # p = 0.02491, rho = -0.08810244  # PAPER
cor.test(contrasts$GenerationLength_d, contrasts$FrC, method= 'spearman')  # p = 0.02121, rho =  0.09050626  # PAPER
cor.test(contrasts$GenerationLength_d, contrasts$FrA, method= 'spearman')  # p = 0.09415
cor.test(contrasts$GenerationLength_d, contrasts$FrG, method= 'spearman')  # p = 0.06404

### start from pairwise correlations and go to multiple linear model:
cor.test(log2(AGG$GenerationLength_d),AGG$FrA, method = 'spearman') # rho -0.2681362; p = 3.635e-12
cor.test(log2(AGG$GenerationLength_d),AGG$FrT, method = 'spearman') # rho -0.3066279; p = 1.287e-15
cor.test(log2(AGG$GenerationLength_d),AGG$FrG, method = 'spearman') # rho  0.1804395; p = 3.665e-06
cor.test(log2(AGG$GenerationLength_d),AGG$FrC, method = 'spearman') # rho  0.4717114;  p < 2.2e-16

cor.test(log2(contrasts$GenerationLength_d), contrasts$FrA, method = 'spearman') # rho -0.07057586; p = 0.1659
cor.test(log2(contrasts$GenerationLength_d), contrasts$FrT, method = 'spearman') # rho -0.1002072; p = 0.04885
cor.test(log2(contrasts$GenerationLength_d), contrasts$FrG, method = 'spearman') # rho  0.03871029; p = 0.4476
cor.test(log2(contrasts$GenerationLength_d), contrasts$FrC, method = 'spearman') # rho  0.1187421;  p = 0.01946

### for multiple linear backward model we choosed A,T & C, on the second step we delete A from the model and the final model is just with T and C

A <- lm(log2(AGG$GenerationLength_d) ~ AGG$FrA + AGG$FrT + AGG$FrC); summary(A) # A is not significant - delete it
A <- lm(log2(AGG$GenerationLength_d) ~ AGG$FrT + AGG$FrC); summary(A)
A <- lm(log2(AGG$GenerationLength_d) ~ scale(AGG$FrT) + scale(AGG$FrC)); summary(A) # log2(GT) = 11.08294 - 0.12 scale(FrT) + 0.45 scale(FrC)
# Adjusted R-squared:  0.2321 - not bad!!!

A <- lm(scale(contrasts$GenerationLength_d) ~ contrasts$FrA + contrasts$FrT + contrasts$FrC); summary(A) # A is not significant - delete it
A <- lm(scale(contrasts$GenerationLength_d) ~ contrasts$FrT + contrasts$FrC); summary(A) # only C is significant
A <- lm(scale(contrasts$GenerationLength_d) ~ contrasts$FrC); summary(A) # only C is significant: 0.00596 **
#(Intercept)   0.001477   0.039088   0.038  0.96987   
#contrasts$FrC 1.091700   0.395708   2.759  0.00596 **
A <- lm(scale(contrasts$GenerationLength_d) ~ 0 + contrasts$FrC); summary(A) # only C is significant: 0.00596 **

# plot it:

pdf("../../Body/4Figures/WholeGenomeAnalyses.MutagensWithinMammalsNoOverlap.R.01.pdf", width = 50, height = 30)

#ColG = rainbow(4)[1]
#ColT = rainbow(4)[2]
#ColC = rainbow(4)[3]
#ColA = rainbow(4)[4]

ColG = rgb(0.1,0.1,0.1,0.2)
ColT = rgb(0.1,0.1,1,0.2)
ColC = rgb(0.1,1,0.1,0.2)
ColA = rgb(1,0.1,0.1,0.2)

par(oma = c(2, 2, 0, 0), cex.main = 2, cex.lab = 1.5, cex = 6, pch =16)
plot(log2(AGG$GenerationLength_d),AGG$FrA, col = ColA, ylim = c(0,0.6), xlab = '', ylab = '', main = ''); par(new=TRUE);
plot(log2(AGG$GenerationLength_d),AGG$FrT, col = ColT, ylim = c(0,0.6), xlab = '', ylab = '', main = ''); par(new=TRUE);
plot(log2(AGG$GenerationLength_d),AGG$FrG, col = ColG, ylim = c(0,0.6), xlab = '', ylab = '', main = ''); par(new=TRUE);
plot(log2(AGG$GenerationLength_d),AGG$FrC, col = ColC, ylim = c(0,0.6), xlab = 'log2(Generation Time in days)', ylab = 'Nucleotide Content', main = ''); par(new=TRUE);
abline(h =0.25, lt = 1, col = 'red');
legend("topright",legend=c('A','C','T','G'), col = c(ColA,ColC,ColT,ColG), pch = 16, horiz = FALSE)

grid.newpage()
par(mfrow=c(2,2))
grid.table(res)

par(mfrow=c(2,2), oma = c(0, 0, 2, 0),cex.main = 2, cex.lab = 2, cex = 2, pch = 16)
plot(log2(AGG$GenerationLength_d),AGG$FrA, col = ColA, ylim = c(0,0.6), xlab = 'log2(GT)', main = 'FrA', cex = 2); abline(h =0.25, lt = 2, col = 'red')
a <-lm(AGG$FrA ~ log2(AGG$GenerationLength_d)); summary(a); abline(a, lwd = 2, col = 'red')
plot(log2(AGG$GenerationLength_d),AGG$FrT, col = ColT, ylim = c(0,0.6), xlab = 'log2(GT)', main = 'FrT', cex = 2); abline(h =0.25, lt = 2, col = 'red')
a <-lm(AGG$FrT ~ log2(AGG$GenerationLength_d)); summary(a); abline(a, lwd = 5, col = 'dark blue')
plot(log2(AGG$GenerationLength_d),AGG$FrG, col = ColG, ylim = c(0,0.6), xlab = 'log2(GT)', main = 'FrG', cex = 2); abline(h =0.25, lt = 2, col = 'red')
a <-lm(AGG$FrG ~ log2(AGG$GenerationLength_d)); summary(a); abline(a, lwd = 2, col = 'red')
plot(log2(AGG$GenerationLength_d),AGG$FrC, col = ColC, ylim = c(0,0.6), xlab = 'log2(GT)', main = 'FrC', cex = 2); abline(h =0.25, lt = 2, col = 'red')
a <-lm(AGG$FrC ~ log2(AGG$GenerationLength_d)); summary(a); abline(a, lwd = 5, col = 'dark green')
mtext("log2(GT) = 11.08294 - 0.12 scale(FrT) + 0.45 scale(FrC)", outer = TRUE, cex = 1.5)

########### question 2: which genes better correlate with GT (why T in ATP6,COX3 and ND4 do not correlate with GT and high absolute value - fast replication, no tRNA before them?)
## T is negatively and C is positively (T->C)
VecOfGenes = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB','ND1','ND2') # ND1 ND2 ND6
length(VecOfGenes)
par(mfcol=c(4,13), cex.main = 2, cex.lab = 2, cex = 2, pch = 16)
for (i in 1:length(VecOfGenes))
{ # i = 1
  OneGene = SynNucGT[SynNucGT$Gene == VecOfGenes[i],]
  main = VecOfGenes[i]
  plot(log2(OneGene$GenerationLength_d),OneGene$FrA, col = ColA, main = main, ylim = c(0,0.7), xlab = '', ylab = '')
  plot(log2(OneGene$GenerationLength_d),OneGene$FrT, col = ColT, main = main, ylim = c(0,0.7), xlab = '', ylab = '')
  plot(log2(OneGene$GenerationLength_d),OneGene$FrG, col = ColG, main = main, ylim = c(0,0.7), xlab = '', ylab = '')
  plot(log2(OneGene$GenerationLength_d),OneGene$FrC, col = ColC, main = main, ylim = c(0,0.7), xlab = '', ylab = '')
}
dev.off()

