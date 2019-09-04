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
#library(ape)

#tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")
#data = SynNucGT[which(as.character(SynNucGT$Species) %in% tree$tip.label),]

#df_vec <- as.character(SynNucGT$Species)
#tree_vec <- tree$tip.label

#a <- setdiff(df_vec, tree_vec)
#b <- setdiff(tree_vec, df_vec)

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
library(pacman) # install.packages("pacman")

tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")
data = AGG[which(as.character(AGG$Species) %in% tree$tip.label),]
# data$GenerationLength_d = log2(data$GenerationLength_d)  # I added log of generation length
df_vec <- as.character(AGG$Species)
tree_vec <- tree$tip.label
a <- setdiff(df_vec, tree_vec)
b <- setdiff(tree_vec, df_vec)
data = data[-597,] # It is duplicate
row.names(data) = data$Species
tree2 <- drop.tip(tree, b)

matchedData = data[match(tree2$tip.label, data$Species),]

# Old
# Phylogenetic tree: tree2
#
#   Number of tips: 705
#   Number of nodes: 704
#   Branch lengths:
#     mean: 0.02279781
#     variance: 0.0007592831
#     distribution summary:
#        Min.     1st Qu.      Median     3rd Qu.        Max.
# 0.000002000 0.005374533 0.013691633 0.029818258 0.207282537
#   No root edge.
#   First ten tip labels: Acinonyx_jubatus
#                         Caracal_caracal
#                         Leptailurus_serval
#                         Felis_silvestris
#                         Felis_margarita
#                         Felis_chaus
#                         Felis_nigripes
#                         Prionailurus_bengalensis
#                         Prionailurus_viverrinus
#                         Prionailurus_rubiginosus
#   No node labels.
# New
# Phylogenetic tree: .
#
#   Number of tips: 649
#   Number of nodes: 648
#   Branch lengths:
#     mean: 0.02377763
#     variance: 0.000807141
#     distribution summary:
#        Min.     1st Qu.      Median     3rd Qu.        Max.
# 0.000002000 0.005784191 0.014175031 0.031247906 0.207282537
#   No root edge.
#   First ten tip labels: Acinonyx_jubatus
#                         Caracal_caracal
#                         Leptailurus_serval
#                         Felis_silvestris
#                         Felis_margarita
#                         Felis_chaus
#                         Felis_nigripes
#                         Prionailurus_bengalensis
#                         Prionailurus_viverrinus
#                         Prionailurus_rubiginosus
#   No node labels.
# TempData = data[, -1]
 # TempData %>% skim()
# n obs: 649
# n variables: 5
#
# -- Variable type:numeric -------------------------------------------------------
#           variable missing complete   n     mean       sd       p0      p25      p50      p75     p100     hist
#                FrA       0      649 649    0.49     0.044   0.34      0.46     0.49     0.52      0.64 ▁▁▃▇▇▃▁▁
#                FrC       0      649 649    0.27     0.049   0.12      0.24     0.27     0.29      0.44 ▁▂▃▇▆▂▁▁
#                FrG       0      649 649    0.046    0.018   0.0097    0.033    0.042    0.056     0.1  ▁▅▇▅▃▂▁▁
#                FrT       0      649 649    0.19     0.044   0.096     0.16     0.19     0.22      0.33 ▁▃▇▆▅▂▁▁
# GenerationLength_d       0      649 649 2799.22  2116.86  341.27   1458.08  2190     3650     18980    ▇▃▁▁▁▁▁▁


# generlen       0      705 705 2790.31 2074.87 341.27 1460 2227.96 3650 18980 ▇▃▁▁▁▁▁▁




# p_load(tibble, dplyr, magrittr, purrr, skimr)
# contrasts <- TempData %>%
#  select(GenerationLength_d, FrA, FrT, FrG, FrC) %>%
#   mutate_if(is.numeric, log2) %>%
#   map(pic, tree2)

matchedData$GenerationLength_d = log2(matchedData$GenerationLength_d)

p_load(tibble, dplyr, magrittr, purrr, skimr)
contrasts <- matchedData[, -1] %>%
#  select(GenerationLength_d, FrA, FrT, FrG, FrC) %>%
#  mutate_if(is.numeric, log2) %>%
  map(pic, tree2)

  # n obs: 648
  # n variables: 5
  #
  # -- Variable type:numeric -------------------------------------------------------
  #           variable missing complete   n    mean   sd     p0   p25    p50     p75     p100     hist
  #                FrA       0      648 648  0.0086 0.82  -4.22 -0.32  0.028  0.36       3.15 ▁▁▁▂▇▂▁▁
  #                FrC       0      648 648 -0.028  1.7  -10.13 -0.65 -0.035  0.59       9.15 ▁▁▁▃▇▁▁▁
  #                FrG       0      648 648  0.24   3.98 -17.33 -1.62  0.11   1.95      23.51 ▁▁▂▇▂▁▁▁
  #                FrT       0      648 648 -0.042  2.26 -13.84 -1    -0.069  0.72      14.23 ▁▁▁▇▆▁▁▁
  # GenerationLength_d       0      648 648 -0.2    0.46  -4.74 -0.19 -0.056 -0.0013 2e-14    ▁▁▁▁▁▁▁▇

# contrasts2 <- as.data.frame(apply(TempData, 2, pic, tree2))
# names(contrasts2) = names(TempData)

# summary(pic(log2(matchedData$GenerationLength_d), tree2)) == summary(contrasts$GenerationLength_d)
summary(pic(matchedData$GenerationLength_d, tree2)) == summary(contrasts$GenerationLength_d)

############################################################################

library(caper)


comparTable = comparative.data(tree2, matchedData, Species)

crunch(scale(GenerationLength_d) ~ scale(FrT), data=comparTable) # est -0.10025, p 0.0172
crunch(scale(GenerationLength_d) ~ scale(FrC), data=comparTable) # est 0.12110 p 0.00591
crunch(GenerationLength_d ~ FrA, data=comparTable)
crunch(GenerationLength_d ~ FrG, data=comparTable)

# crunch(scale(GenerationLength_d) ~ scale(FrT), data=comparTable, node.depth = 2)
# crunch(scale(GenerationLength_d) ~ scale(FrC), data=comparTable, node.depth = 2)

crunch(GenerationLength_d ~ FrT + FrC + FrA, data=comparTable) # the same for contrasts
crunch(scale(GenerationLength_d) ~ FrT + FrC, data=comparTable)


cor.test(contrasts$GenerationLength_d, contrasts$FrT, method= 'spearman')  # p = 0.02491, rho = -0.08810244  # PAPER
cor.test(contrasts$GenerationLength_d, contrasts$FrC, method= 'spearman')  # p = 0.02121, rho =  0.09050626  # PAPER
cor.test(contrasts$GenerationLength_d, contrasts$FrA, method= 'spearman')  # p = 0.09415
cor.test(contrasts$GenerationLength_d, contrasts$FrG, method= 'spearman')  # p = 0.06404

### start from pairwise correlations and go to multiple linear model:
cor.test(log2(AGG$GenerationLength_d),AGG$FrA, method = 'spearman') # rho -0.2681362; p = 3.635e-12
cor.test(log2(AGG$GenerationLength_d),AGG$FrT, method = 'spearman') # rho -0.3066279; p = 1.287e-15
cor.test(log2(AGG$GenerationLength_d),AGG$FrG, method = 'spearman') # rho  0.1804395; p = 3.665e-06
cor.test(log2(AGG$GenerationLength_d),AGG$FrC, method = 'spearman') # rho  0.4717114;  p < 2.2e-16

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

########### question 3: pairwise correlations between nucleotide content - which nucleotides anticorrelate with each other better?
### in this way we can try to move from 4 dimensions to 6 dimensions

pdf("../../Body/4Figures/WholeGenomeAnalyses.MutagensWithinMammalsNoOverlap.Tables.R.01.pdf")

Res = c()
AT = cor.test(AGG$FrA,AGG$FrT, method = 'spearman'); Res = c('AT', as.numeric(AT[3]), as.numeric(AT[4]));
AG = cor.test(AGG$FrA,AGG$FrG, method = 'spearman'); Res = rbind(Res, c('AG', as.numeric(AG[3]), as.numeric(AG[4])));
AC = cor.test(AGG$FrA,AGG$FrC, method = 'spearman'); Res = rbind(Res, c('AC', as.numeric(AC[3]), as.numeric(AC[4])));
TC = cor.test(AGG$FrT,AGG$FrC, method = 'spearman'); Res = rbind(Res, c('TC', as.numeric(TC[3]), as.numeric(TC[4])));
TG = cor.test(AGG$FrT,AGG$FrG, method = 'spearman'); Res = rbind(Res, c('TG', as.numeric(TG[3]), as.numeric(TG[4])));
CG = cor.test(AGG$FrC,AGG$FrG, method = 'spearman'); Res = rbind(Res, c('CG', as.numeric(CG[3]), as.numeric(CG[4])));
names(Res) = c('Subst', 'Pvalue','SpearmanRho')

grid.newpage()
grid.table(Res)

# the same with contrasts:
ResC = c()
AT = cor.test(contrasts$FrA,contrasts$FrT, method = 'spearman'); ResC = c('AT', as.numeric(AT[3]), as.numeric(AT[4]));
AG = cor.test(contrasts$FrA,contrasts$FrG, method = 'spearman'); ResC = rbind(ResC, c('AG', as.numeric(AG[3]), as.numeric(AG[4])));
AC = cor.test(contrasts$FrA,contrasts$FrC, method = 'spearman'); ResC = rbind(ResC, c('AC', as.numeric(AC[3]), as.numeric(AC[4])));
TC = cor.test(contrasts$FrT,contrasts$FrC, method = 'spearman'); ResC = rbind(ResC, c('TC', as.numeric(TC[3]), as.numeric(TC[4])));
TG = cor.test(contrasts$FrT,contrasts$FrG, method = 'spearman'); ResC = rbind(ResC, c('TG', as.numeric(TG[3]), as.numeric(TG[4])));
CG = cor.test(contrasts$FrC,contrasts$FrG, method = 'spearman'); ResC = rbind(ResC, c('CG', as.numeric(CG[3]), as.numeric(CG[4])));
names(ResC) = c('Subst', 'Pvalue','SpearmanRho')

grid.newpage()
grid.table(ResC)
dev.off()
