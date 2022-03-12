rm(list=ls(all=TRUE))
library(dplyr)

SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
SynNuc = SynNuc[SynNuc$Gene != 'ND6',]

names(SynNuc)

GT = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GT$Species = gsub(' ','_',GT$Scientific_name)

GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

Mam = merge(SynNuc,GT, by ='Species')

names(Mam)

# calculate freq of each nucleotide in each position
vec_all = c('TTC','TTT','TCC','TCT','TAC','TAT','TGC','TGT',
            'TTA','TTG','TCA','TCG','TAA','TAG','TGA','TGG',
            'CTC','CTT','CCC','CCT','CAC','CAT','CGC','CGT',
            'CTA','CTG','CCA','CCG','CAA','CAG','CGA','CGG',
            'ATC','ATT','ACC','ACT','AAC','AAT','AGC','AGT',
            'ATA','ATG','ACA','ACG','AAA','AAG','AGA','AGG',
            'GTC','GTT','GCC','GCT','GAC','GAT','GGC','GGT',
            'GTA','GTG','GCA','GCG','GAA','GAG','GGA','GGG')


sp_sum_gen = data.frame(unique(Mam$Species))

for (codon in vec_all)
{
  
  sum_of_codon = aggregate(Mam[ ,codon], by = list(Mam$Species), FUN = 'sum')[2]
  sp_sum_gen = cbind(sp_sum_gen, sum_of_codon)
  
}

names(sp_sum_gen) = c('Species', vec_all)

test = sp_sum_gen %>% 
  mutate(AXX = select(., c('ATC','ATT','ACC','ACT','AAC','AAT','AGC','AGT',
                         'ATA','ATG','ACA','ACG','AAA','AAG','AGA','AGG')) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(TXX = select(., c('TTC','TTT','TCC','TCT','TAC','TAT','TGC','TGT',
                           'TTA','TTG','TCA','TCG','TAA','TAG','TGA','TGG')) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(CXX = select(., c('CTC','CTT','CCC','CCT','CAC','CAT','CGC','CGT',
                           'CTA','CTG','CCA','CCG','CAA','CAG','CGA','CGG')) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(GXX = select(., c('GTC','GTT','GCC','GCT','GAC','GAT','GGC','GGT',
                           'GTA','GTG','GCA','GCG','GAA','GAG','GGA','GGG')) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(XAX = select(., c('TAC','TAT','TAA','TAG','CAC','CAT','CAA','CAG',
                           'AAC','AAT','AAA','AAG','GAC','GAT','GAA','GAG')) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(XTX = select(., c('TTC','TTT','TTA','TTG','CTC','CTT','CTA','CTG',
                           'ATC','ATT','ATA','ATG','GTC','GTT','GTA','GTG')) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(XCX = select(., c('TCC','TCT','TCA','TCG','CCC','CCT','CCA','CCG',
                           'ACC','ACT','ACA','ACG','GCC','GCT','GCA','GCG')) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(XGX = select(., c('TGC','TGT','TGA','TGG','CGC','CGT','CGA','CGG',
                           'AGC','AGT','AGA','AGG','GGC','GGT','GGA','GGG')) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(XXA = select(., c('TTA','TCA','TAA','TGA','CTA','CCA','CAA','CGA',
                           'ATA','ACA','AAA','AGA','GTA','GCA','GAA','GGA')) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(XXT = select(., c('TTT','TCT','TAT','TGT','CTT','CCT','CAT','CGT',
                           'ATT','ACT','AAT','AGT','GTT','GCT','GAT','GGT')) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(XXC = select(., c('TTC','TCC','TAC','TGC','CTC','CCC','CAC','CGC',
                           'ATC','ACC','AAC','AGC','GTC','GCC','GAC','GGC')) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(XXG = select(., c('TTG','TCG','TAG','TGG','CTG','CCG','CAG','CGG',
                           'ATG','ACG','AAG','AGG','GTG','GCG','GAG','GGG')) %>% rowSums(na.rm = TRUE))

test = test[,-c(2:65)]

test$triplet_sum_1 = test$AXX + test$TXX + test$CXX + test$GXX
test$triplet_sum_2 = test$XAX + test$XTX + test$XCX + test$XGX  
test$triplet_sum_3 = test$XXA + test$XXT + test$XXC + test$XXG    

test$norm_AXX = test$AXX / test$triplet_sum_1
test$norm_TXX = test$TXX / test$triplet_sum_1
test$norm_CXX = test$CXX / test$triplet_sum_1
test$norm_GXX = test$GXX / test$triplet_sum_1

test$norm_XAX = test$XAX / test$triplet_sum_2
test$norm_XTX = test$XTX / test$triplet_sum_2
test$norm_XCX = test$XCX / test$triplet_sum_2
test$norm_XGX = test$XGX / test$triplet_sum_2

test$norm_XXA = test$XXA / test$triplet_sum_3
test$norm_XXT = test$XXT / test$triplet_sum_3
test$norm_XXC = test$XXC / test$triplet_sum_3
test$norm_XXG = test$XXG / test$triplet_sum_3


final = merge(test, GT, by ='Species')


library(ape)
library(geiger)
library(caper)

tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")

Mam = final
Mam = Mam[Mam$Species != 'Neophocaena_phocaenoides',]
row.names(Mam) = Mam$Species


tree_w = treedata(tree, Mam, sort=T, warnings=T)$phy

data<-as.data.frame(treedata(tree_w, Mam, sort=T, warnings=T)$data)

data$Species = as.character(data$Species)

data$TXX = as.numeric(as.character(data$TXX))
data$XTX = as.numeric(as.character(data$XTX))
data$XXT = as.numeric(as.character(data$XXT))

data$CXX = as.numeric(as.character(data$CXX))
data$XCX = as.numeric(as.character(data$XCX))
data$XXC = as.numeric(as.character(data$XXC))

data$GXX = as.numeric(as.character(data$GXX))
data$XGX = as.numeric(as.character(data$XGX))
data$XXG = as.numeric(as.character(data$XXG))

data$AXX = as.numeric(as.character(data$AXX))
data$XAX = as.numeric(as.character(data$XAX))
data$XXA = as.numeric(as.character(data$XXA))


data$GenerationLength_d = as.numeric(as.character(data$GenerationLength_d))

### PGLS WITH posotions

MutComp = comparative.data(tree_w, data, Species, vcv=TRUE)

summary(pgls(TXX ~ log2(GenerationLength_d), MutComp, lambda="ML"))
summary(pgls(XTX ~ log2(GenerationLength_d), MutComp, lambda="ML"))
summary(pgls(XXT ~ log2(GenerationLength_d), MutComp, lambda="ML"))


summary(pgls(CXX ~ log2(GenerationLength_d), MutComp, lambda="ML"))
summary(pgls(XCX ~ log2(GenerationLength_d), MutComp, lambda="ML"))
summary(pgls(XXC ~ log2(GenerationLength_d), MutComp, lambda="ML"))


summary(pgls(GXX ~ log2(GenerationLength_d), MutComp, lambda="ML"))
summary(pgls(XGX ~ log2(GenerationLength_d), MutComp, lambda="ML"))
summary(pgls(XXG ~ log2(GenerationLength_d), MutComp, lambda="ML"))


summary(pgls(AXX ~ log2(GenerationLength_d), MutComp, lambda="ML"))
summary(pgls(XAX ~ log2(GenerationLength_d), MutComp, lambda="ML"))
summary(pgls(XXA ~ log2(GenerationLength_d), MutComp, lambda="ML"))


summary(pgls(log2(GenerationLength_d) ~ TXX + XXT, MutComp, lambda="ML"))
summary(pgls(log2(GenerationLength_d) ~ TXX * XXT, MutComp, lambda="ML"))

summary(pgls(log2(GenerationLength_d) ~ CXX + XXC, MutComp, lambda="ML"))
summary(pgls(log2(GenerationLength_d) ~ CXX * XXC, MutComp, lambda="ML"))
