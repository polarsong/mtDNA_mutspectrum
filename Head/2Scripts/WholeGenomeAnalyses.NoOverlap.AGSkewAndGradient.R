################################
################################

rm(list=ls(all=TRUE))

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")){file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

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

VecOfTaxa = unique(SynNuc$Class)
table(SynNuc$Class)/13

AGG = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species,SynNuc$Class,SynNuc$Gene), FUN = sum)
names(AGG) = c('Species','Class','Gene','NeutralA','NeutralT','NeutralG','NeutralC')
AGG$TCSkew = (AGG$NeutralT - AGG$NeutralC)/(AGG$NeutralT + AGG$NeutralC)



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
boxplot(TCSkew ~ GT*Gene, data = M,  notch = TRUE, outline = FALSE, las = 2, col = c('red','green'), ylim = c(-1,0), main = 'Mammalia, AG skew')


dev.off()
pdf("../../Body/4Figures/WholeGenomeAnalyses.NoOverlap.AGSkew.R.pdf", height = 20, width = 40)
par(mfcol=c(2,3), cex = 2)
boxplot(AGG[AGG$Class == 'Actinopterygii',]$TCSkew,AGG[AGG$Class == 'Amphibia',]$TCSkew,AGG[AGG$Class == 'Reptilia',]$TCSkew,AGG[AGG$Class == 'Mammalia',]$TCSkew,AGG[AGG$Class == 'Aves',]$TCSkew, notch = TRUE, outline = FALSE, names = c('Actinopterygii','Amphibia','Reptilia','Mammalia','Aves'), main = 'AG skew', las = 1)


############# GENERATION LENGTH FOR MAMMALS
GT = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

Mam = merge(AGG,GT, by ='Species')
plot(log2(Mam$GenerationLength_d),Mam$TCSkew, main = 'Mammalia', xlab = 'log2(Generation Length)', ylab = 'AG skew', ylim=c(-0.12,0.7))
a<-lm(Mam$TCSkew ~ log2(Mam$GenerationLength_d)); summary(a)
abline(a, col = 'red', lwd = 4)


cor.test(log2(Mam$GenerationLength_d),Mam$TCSkew, method = 'spearman') 
median(log2(Mam$GenerationLength_d)) # 11.09672
cor.test(log2(Mam[log2(Mam$GenerationLength_d) > 11.09672,]$GenerationLength_d),Mam[log2(Mam$GenerationLength_d) > 11.09672,]$TCSkew, method = 'spearman') 
cor.test(log2(Mam[log2(Mam$GenerationLength_d) < 11.09672,]$GenerationLength_d),Mam[log2(Mam$GenerationLength_d) < 11.09672,]$TCSkew, method = 'spearman') 


dev.off()



