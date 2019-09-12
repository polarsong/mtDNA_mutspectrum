################################
################################

rm(list=ls(all=TRUE))

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")){file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

names(SynNuc)

### make ND6 complementary:
#NotND6 = SynNuc[SynNuc$Gene != 'ND6',]
#ND6 = SynNuc[SynNuc$Gene == 'ND6',]
#A = ND6$NeutralT
#T = ND6$NeutralA
#G = ND6$NeutralC
#C = ND6$NeutralG
#ND6$NeutralA = A
#ND6$NeutralT = T
#ND6$NeutralG = G
#ND6$NeutralC = C
#SynNuc = rbind(NotND6,ND6)

VecOfTaxa = unique(SynNuc$Class)
table(SynNuc$Class)/13

AGG = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species,SynNuc$Class,SynNuc$Gene), FUN = sum)
names(AGG) = c('Species','Class','Gene','NeutralA','NeutralT','NeutralG','NeutralC')
AGG$AGSkew = (AGG$NeutralA - AGG$NeutralG)/(AGG$NeutralA + AGG$NeutralG)



pdf("../../Body/4Figures/WholeGenomeAnalyses.NoOverlap.AGSkew.Gradient.R.pdf", height = 20, width = 40)

############# GENERATION LENGTH FOR MAMMALS
GT = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

M = merge(AGG,GT, by ='Species')
summary(M$GenerationLength_d) # median = 2190.0
ShortLived = unique(M[M$GenerationLength_d <= median(M$GenerationLength_d),]$Species); length(ShortLived)
LongLived = unique(M[M$GenerationLength_d  > median(M$GenerationLength_d),]$Species);  length(LongLived)
MShort = M[M$Species %in% ShortLived,]; MShort$GT = 'short'; MLong = M[M$Species %in% LongLived,]; MLong$GT = 'long';
M = rbind(MShort,MLong)

M$Gene =  ordered(M$Gene, levels = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB', 'ND1','ND2'))
# M$Gene =  ordered(M$Gene, levels = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB'))
M = M[order(M$Gene),]

par(mfrow=c(2,1), oma = c(3, 1, 1, 1), cex = 2)
boxplot(AGSkew ~ GT*Gene, data = M,  notch = TRUE, outline = FALSE, las = 2, col = c('red','green'), ylim = c(-1,0), main = 'Mammalia, AG skew')


dev.off()

