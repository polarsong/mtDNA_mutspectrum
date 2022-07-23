rm(list=ls(all=TRUE))

##### 1: read table with fractions of all mutations for all vertebrates and filter out species with TsTv == Inf (Tv == 0)

MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegCytb.txt', header = TRUE) # OnlyFourFoldDegCytb
MUT$TsTv = (MUT$T_C + MUT$C_T + MUT$G_A + MUT$A_G) / (MUT$T_A + MUT$A_T + MUT$G_C + MUT$C_G + MUT$G_T + MUT$T_G + MUT$C_A + MUT$A_C)
MUT$TC_TATGTC = (MUT$T_C) / (MUT$T_A + MUT$T_G + MUT$T_C)
summary(MUT$TC_TATGTC)
names(MUT)
summary(MUT$TsTv)
dim(MUT)
MUT = MUT[MUT$TsTv < Inf,]
MUT = MUT[!is.na(MUT$TC_TATGTC),]
dim(MUT)

##### 2: read GenLength table (only mammals) and merge with Mut to get MamGt

GT = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', sep = '\t', header = TRUE)
GT$Species = gsub(' ','_',GT$Scientific_name)

dim(MUT)
MamGt = merge(GT,MUT, by = 'Species')
dim(MamGt)
Qual = data.frame(table(MamGt$Species))
UniqueSpecies = Qual[Qual$Freq == 1,]$Var1; length(UniqueSpecies);
MamGt = MamGt[MamGt$Species %in% UniqueSpecies,]
dim(MamGt)
MamGt$TC_TCGA = MamGt$T_C / (MamGt$T_C + MamGt$G_A)

##### 3: read the number of mutations per species per gene and keep in MamGt only species with at least 15 fourfold mutations in Cytb
Number = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.txt')
Number = data.frame(Number$Species,Number$NumOfFourFoldMutInCytB); names(Number)=c('Species','NumOfFourFoldMutInCytB')
Number = unique(Number)
nrow(MamGt)
MamGt = merge(MamGt,Number, by ='Species')
nrow(MamGt)

summary(MamGt$NumOfFourFoldMutInCytB); # dev.off()
MamGt = MamGt[MamGt$NumOfFourFoldMutInCytB >= 15,]
dim(MamGt) 
hist(MamGt$NumOfFourFoldMutInCytB,breaks=100)

##### 4: key cor.test analyses (T_C ~ GT) for the whole dataset and 3 quartiles: all are significant and positive
names(MamGt)
summary(MamGt$T_C)
summary(MamGt$GenerationLength_d)

# MoreOrEqual36
nrow(MamGt[MamGt$NumOfFourFoldMutInCytB >= 36,])
cor.test(MamGt[MamGt$NumOfFourFoldMutInCytB >= 36,]$T_C,MamGt[MamGt$NumOfFourFoldMutInCytB >= 36,]$GenerationLength_d, method = 'spearman')

# MoreOrEqual60
nrow(MamGt[MamGt$NumOfFourFoldMutInCytB >= 60,])
cor.test(MamGt[MamGt$NumOfFourFoldMutInCytB >= 60,]$T_C,MamGt[MamGt$NumOfFourFoldMutInCytB >= 60,]$GenerationLength_d, method = 'spearman')

# MoreOrEqual84
nrow(MamGt[MamGt$NumOfFourFoldMutInCytB >= 84,])
cor.test(MamGt[MamGt$NumOfFourFoldMutInCytB >= 84,]$T_C,MamGt[MamGt$NumOfFourFoldMutInCytB >= 84,]$GenerationLength_d, method = 'spearman')

# MoreOrEqual108
nrow(MamGt[MamGt$NumOfFourFoldMutInCytB >= 108,])
cor.test(MamGt[MamGt$NumOfFourFoldMutInCytB >= 108,]$T_C,MamGt[MamGt$NumOfFourFoldMutInCytB >= 108,]$GenerationLength_d, method = 'spearman')

# MoreOrEqual132
nrow(MamGt[MamGt$NumOfFourFoldMutInCytB >= 132,])
cor.test(MamGt[MamGt$NumOfFourFoldMutInCytB >= 132,]$T_C,MamGt[MamGt$NumOfFourFoldMutInCytB >= 132,]$GenerationLength_d, method = 'spearman')

# MoreOrEqual156
nrow(MamGt[MamGt$NumOfFourFoldMutInCytB >= 156,])
cor.test(MamGt[MamGt$NumOfFourFoldMutInCytB >= 156,]$T_C,MamGt[MamGt$NumOfFourFoldMutInCytB >= 156,]$GenerationLength_d, method = 'spearman')

##### 5: key PGLS analyses (T_C ~ GT) for the dataset with different thresholds
##### PICs, Alina+Kostya

library(ape)
library(geiger)
library(caper)

names(MamGt)
summary(MamGt$GenerationLength_d)
tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")
row.names(MamGt) = MamGt$Species
tree_w = treedata(tree, MamGt[, c('Species', 'TsTv', 'T_C', 'TC_TCGA', 'G_A','TC_TATGTC','GenerationLength_d','NumOfFourFoldMutInCytB')], 
                sort=T, warnings=T)$phy
data<-as.data.frame(treedata(tree_w, MamGt[, c('Species', 'TsTv', 'T_C', 'TC_TCGA','G_A','TC_TATGTC','GenerationLength_d','NumOfFourFoldMutInCytB')], 
                 sort=T, warnings=T)$data)
nrow(data)
data$Species = as.character(data$Species)
data$TsTv = as.numeric(as.character(data$TsTv))
data$T_C = as.numeric(as.character(data$T_C))
data$G_A = as.numeric(as.character(data$G_A))
data$TC_TCGA = as.numeric(as.character(data$TC_TCGA))
data$TC_TATGTC = as.numeric(as.character(data$TC_TATGTC))
data$GenerationLength_d = as.numeric(as.character(data$GenerationLength_d))
data$NumOfFourFoldMutInCytB = as.numeric(as.character(data$NumOfFourFoldMutInCytB))
nrow(data) # 211

MutCompMoreOrEqual36 = comparative.data(tree_w, data[data$NumOfFourFoldMutInCytB >= 36,], Species, vcv=TRUE)
MutCompMoreOrEqual60 = comparative.data(tree_w, data[data$NumOfFourFoldMutInCytB >= 60,], Species, vcv=TRUE)
MutCompMoreOrEqual84 = comparative.data(tree_w, data[data$NumOfFourFoldMutInCytB >= 84,], Species, vcv=TRUE)
MutCompMoreOrEqual108 = comparative.data(tree_w, data[data$NumOfFourFoldMutInCytB >= 108,], Species, vcv=TRUE)
MutCompMoreOrEqual132 = comparative.data(tree_w, data[data$NumOfFourFoldMutInCytB >= 132,], Species, vcv=TRUE)
MutCompMoreOrEqual156 = comparative.data(tree_w, data[data$NumOfFourFoldMutInCytB >= 156,], Species, vcv=TRUE)

## analyses

summary(pgls(T_C ~ log2(GenerationLength_d), MutCompMoreOrEqual36, lambda="ML")) # intercepts everythere are non-significant, so we run regression through the zero
summary(pgls(T_C ~ 0+log2(GenerationLength_d), MutCompMoreOrEqual36, lambda="ML"))
summary(pgls(T_C ~ 0+log2(GenerationLength_d), MutCompMoreOrEqual60, lambda="ML"))
summary(pgls(T_C ~ 0+log2(GenerationLength_d), MutCompMoreOrEqual84, lambda="ML"))
summary(pgls(T_C ~ 0+log2(GenerationLength_d), MutCompMoreOrEqual108, lambda="ML"))
summary(pgls(T_C ~ 0+log2(GenerationLength_d), MutCompMoreOrEqual132, lambda="ML"))
summary(pgls(T_C ~ 0+log2(GenerationLength_d), MutCompMoreOrEqual156, lambda="ML"))

