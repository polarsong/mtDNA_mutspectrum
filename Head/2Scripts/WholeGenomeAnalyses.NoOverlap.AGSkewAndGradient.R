################################
################################

rm(list=ls(all=TRUE))

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
#if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")){file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

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

########################################## GENOME WIDE SKEW FOR EACH SPECIES

AGG = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species,SynNuc$Class), FUN = sum)
names(AGG) = c('Species','Class','NeutralA','NeutralT','NeutralG','NeutralC')

### count fraction of nucleotides
AGG$FrA = AGG$NeutralA / (AGG$NeutralA + AGG$NeutralT + AGG$NeutralG + AGG$NeutralC)
AGG$FrT = AGG$NeutralT / (AGG$NeutralA + AGG$NeutralT + AGG$NeutralG + AGG$NeutralC) 
AGG$FrG = AGG$NeutralG / (AGG$NeutralA + AGG$NeutralT + AGG$NeutralG + AGG$NeutralC) 
AGG$FrC = AGG$NeutralC / (AGG$NeutralA + AGG$NeutralT + AGG$NeutralG + AGG$NeutralC) 

## all six different skews
AGG$CTSkew = (AGG$NeutralC - AGG$NeutralT)/(AGG$NeutralC + AGG$NeutralT); summary(AGG$CTSkew) # GA on heavy
AGG$CGSkew = (AGG$NeutralC - AGG$NeutralG)/(AGG$NeutralC + AGG$NeutralG); summary(AGG$CGSkew) # 
AGG$CASkew = (AGG$NeutralC - AGG$NeutralA)/(AGG$NeutralC + AGG$NeutralA); summary(AGG$CASkew) # 
AGG$TGSkew = (AGG$NeutralT - AGG$NeutralG)/(AGG$NeutralT + AGG$NeutralG); summary(AGG$TGSkew) # 
AGG$TASkew = (AGG$NeutralT - AGG$NeutralA)/(AGG$NeutralT + AGG$NeutralA); summary(AGG$TASkew) # 
AGG$GASkew = (AGG$NeutralG - AGG$NeutralA)/(AGG$NeutralG + AGG$NeutralA); summary(AGG$GASkew) # 

AGG$TCSkew = (AGG$NeutralT - AGG$NeutralC)/(AGG$NeutralT + AGG$NeutralC); summary(AGG$CTSkew) # AG on heavy. Added it for fun, just to be sure, that it is opposite to CT (GA on heavy) 

#######opening MutSpec data
MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)
MutCTskew = merge(MUT, AGG)
######mut spec with nuc content analyses
cor.test(MutCTskew$T_C, MutCTskew$CTSkew, method = "spearman")


GT = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

Mam = merge(AGG,GT, by ='Species')
MutCTskew = merge(MutCTskew, GT)
summary(lm(CTSkew ~ T_C, data=MutCTskew))
summary(lm(T_C ~ FrT + FrC, data=MutCTskew))
cor.test(MutCTskew$T_C, MutCTskew$FrT, method = "spearman")
cor.test(MutCTskew$T_C, MutCTskew$FrC, method = "spearman")
MutCTskew$FrCFrT = MutCTskew$FrC/MutCTskew$FrT
cor.test(MutCTskew$T_C, MutCTskew$FrCFrT, method = "spearman")

############### generation length versus all six skews: 6 rank correlations and one multiple model:
###### pairwise rank corr:
cor.test(log2(Mam$GenerationLength_d),Mam$CTSkew, method = 'spearman')  # rho = 0.4117201; p-value < 2.2e-16. I would add PIC
cor.test(log2(Mam$GenerationLength_d),Mam$CGSkew, method = 'spearman')  # rho = 0.05774663, p = 0.1414
cor.test(log2(Mam$GenerationLength_d),Mam$CASkew, method = 'spearman')  # rho = 0.4691996,  p-value < 2.2e-16
cor.test(log2(Mam$GenerationLength_d),Mam$TGSkew, method = 'spearman')  # rho = -0.24,  p-value = 5.183e-10
cor.test(log2(Mam$GenerationLength_d),Mam$TASkew, method = 'spearman')  # rho = -0.1183936,  p-value = 0.0025
cor.test(log2(Mam$GenerationLength_d),Mam$GASkew, method = 'spearman')  # rho = 0.2176591,  p-value = 2.072e-08

######### multiple Linear Model (it is reasonable to add only 4 significant parameters - without CG and TA. But If I add everything, results are unusual and close to opposite...). Interactions?
A<-lm(log2(Mam$GenerationLength_d) ~ Mam$CTSkew + Mam$CGSkew + Mam$CASkew + Mam$TGSkew + Mam$TASkew + Mam$GASkew); summary(A)
# (Intercept)   10.394      1.192   8.722  < 2e-16 ***
# Mam$CTSkew     2.076      3.054   0.680   0.4969    
# Mam$CGSkew   -10.172      2.391  -4.254 2.42e-05 ***
# Mam$CASkew     7.169      3.450   2.078   0.0381 *  
# Mam$TGSkew     3.531      1.567   2.253   0.0246 *  
# Mam$TASkew    -1.858      3.502  -0.531   0.5958    
# Mam$GASkew    -7.859      4.083  -1.925   0.0547 .  

A<-lm(log2(Mam$GenerationLength_d) ~ Mam$CTSkew +  Mam$CASkew + Mam$TGSkew + Mam$TASkew + Mam$GASkew); summary(A)
A<-lm(log2(Mam$GenerationLength_d) ~ Mam$CTSkew +  Mam$CASkew + Mam$TGSkew + Mam$GASkew); summary(A)
#  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   12.917      1.047  12.342  < 2e-16 ***
#  Mam$CTSkew     3.187      1.045   3.049  0.00239 ** 
#  Mam$CASkew     1.268      1.116   1.136  0.25619    
# Mam$TGSkew     2.700      1.558   1.733  0.08361 .  
# Mam$GASkew     4.658      2.874   1.621  0.10559    
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.9818 on 645 degrees of freedom
# Multiple R-squared:  0.2435,	Adjusted R-squared:  0.2388 
# F-statistic: 51.91 on 4 and 645 DF,  p-value: < 2.2e-16


##################### PICs

library(ape)
library(geiger)
library(caper)

tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")

Mam = Mam[-404,]
row.names(Mam) = Mam$Species

tree_w = treedata(tree, Mam, sort=T, warnings=T)$phy

data<-as.data.frame(treedata(tree_w, Mam, sort=T, warnings=T)$data)

data$Species = as.character(data$Species)

data$CTSkew = as.numeric(as.character(data$CTSkew))
data$GenerationLength_d = as.numeric(as.character(data$GenerationLength_d))

cor.test(pic(data$CTSkew, tree_w), pic(data$GenerationLength_d, tree_w), method = 'spearman')

# rho 
# 0.0942479 
# p-value = 0.0164

MutComp = comparative.data(tree_w, data, Species, vcv=TRUE)

model = pgls(scale(GenerationLength_d) ~ scale(CTSkew), MutComp, lambda="ML")
summary(model)

# lambda [ ML]  : 0.920
# (Intercept)   0.037040   0.387705  0.0955 0.9239175    
# scale(CTSkew) 0.134009   0.036792  3.6423 0.0002918 ***

crunch(scale(GenerationLength_d) ~ scale(CTSkew), MutComp)

# S = 41138230, p-value = 0.01806
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.09286304 

###############################3### more PICs
library(ape)
library(geiger)
library(caper)

MutCTskew = MutCTskew[-156,]
MutCTskew$FrCFrT = MutCTskew$FrC/MutCTskew$FrT
tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")
row.names(MutCTskew) = MutCTskew$Species
tree_w = treedata(tree, MutCTskew, sort=T, warnings=T)$phy


data<-as.data.frame(treedata(tree_w, MutCTskew, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)
data$T_C = as.numeric(as.character(data$T_C))
data$FrT = as.numeric(as.character(data$FrT))
data$FrC = as.numeric(as.character(data$FrC))
data$FrCFrT = as.numeric(as.character(data$FrCFrT))
data$CTSkew = as.numeric(as.character(data$CTSkew))

MutComp = comparative.data(tree_w, data, Species, vcv=TRUE)
summary(pgls(T_C ~ FrT + FrC, MutComp, lambda="ML"))
summary(pgls(T_C ~ FrCFrT, MutComp, lambda="ML"))
summary(pgls(T_C ~ FrC, MutComp, lambda="ML"))
summary(pgls(T_C ~ FrT, MutComp, lambda="ML"))
summary(pgls(T_C ~ FrT * FrC, MutComp, lambda="ML"))
summary(pgls(T_C ~ CTSkew, MutComp, lambda="ML"))

str(MutComp)

###### plot 
pdf("../../Body/4Figures/WholeGenomeAnalyses.NoOverlap.AGSkew.R.01.pdf", height = 10, width = 20)

plot(log2(Mam$GenerationLength_d),Mam$CTSkew, main = 'Mammalia', xlab = 'log2(Generation Length)', ylab = 'GA skew')
a<-lm(Mam$CTSkew ~ log2(Mam$GenerationLength_d)); summary(a)
abline(a, col = 'red', lwd = 4)
abline(h=0, col = 'grey', lwd = 2, lt = 2)

#################################### GENE_SPECIFIC SKEW FOR EACH SPECIES

AGG = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species,SynNuc$Class,SynNuc$Gene), FUN = sum)
names(AGG) = c('Species','Class','Gene','NeutralA','NeutralT','NeutralG','NeutralC')
AGG$CTSkew = (AGG$NeutralC - AGG$NeutralT)/(AGG$NeutralC + AGG$NeutralT); summary(AGG$TCSkew) # GA on heavy

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
MShort = M[M$Species %in% ShortLived,]; MShort$GT = 'short lived'; MLong = M[M$Species %in% LongLived,]; MLong$GT = 'long lived';
M = rbind(MShort,MLong)

M$Gene =  ordered(M$Gene, levels = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB', 'ND1','ND2'))
# M$Gene =  ordered(M$Gene, levels = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB'))
M = M[order(M$Gene),]

# par(mfrow=c(2,1), oma = c(3, 1, 1, 1), cex = 2)
boxplot(CTSkew ~ GT*Gene, data = M,  notch = TRUE, outline = FALSE, las = 2, col = c('red','green'), main = 'Mammalia, GA skew')

M = M[!M$Gene %in% c('ND6','ND1','ND2'),]
M$Gene =  ordered(M$Gene, levels = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB'))
boxplot(CTSkew ~ GT*Gene, data = M,  notch = TRUE, outline = FALSE, las = 2, col = c('red','green'), main = 'Mammalia, GA skew')


########By Order analyses######################
###############################################
library(dplyr)
Taxa = read.table("../../Body/2Derived/MammalianTaxonomy.txt", header = TRUE)
TaxaSpOrd=data.frame(Taxa$Species, Taxa$Order)
names(TaxaSpOrd) = c("Species", "Order")
M = merge(M, TaxaSpOrd)
Orders = unique(M$Order)
table(M$Order)

#Afrotheria       Caniformia          Cetacea   Dasyuromorphia  Didelphimorphia    Diprotodontia 
#110              240              430               50               50              120 
#Euarchontoglires       Feliformia      Haplorrhini   Laurasiatheria   Megachiroptera   Microbiotheria 
#150              370              840             1260               70               10 
#Microchiroptera      Monotremata Notoryctemorphia Paucituberculata  Peramelemorphia         Rodentia 
#160               20               10               20               30              910 
#Ruminantia    Strepsirrhini        Xenarthra 
#1120              250              230 

for ( i in unique(M$Order)){
  out = ''
  out = paste(i, IQR(M[M$Order == i,]$GenerationLength_d), sep = ' ', nrow(M[M$Order == i,]))
  print(out)
}

#[1] "Feliformia 773.05315 370"
#[1] "Rodentia 478.73066123767 910"
#[1] "Ruminantia 799.3417708285 1120"
#[1] "Laurasiatheria 2447.0747625 1260"
#[1] "Haplorrhini 1150.105875 840"
#[1] "Microchiroptera 662.2593 160"
#[1] "Caniformia 832.217431938 240"
#[1] "Strepsirrhini 1371.624375 250"
#[1] "Cetacea 2664.5 430"
#[1] "Xenarthra 336.6637 230"
#[1] "Paucituberculata 0 20"
#[1] "Afrotheria 6993.4424166715 110"
#[1] "Megachiroptera 1120.91840714286 70"
#[1] "Diprotodontia 896.200516485402 120"
#[1] "Didelphimorphia 103.149975 50"
#[1] "Microbiotheria 0 10"
#[1] "Euarchontoglires 60.1993 150"
#[1] "Peramelemorphia 560.1198 30"
#[1] "Dasyuromorphia 291.8336 50"
#[1] "Notoryctemorphia 0 10"
#[1] "Monotremata 2808.2654 20"

Laurasia = M[M$Order == "Laurasiatheria",]
Laurasia$GL = "Middle"
Laurasia[Laurasia$GenerationLength_d <= quantile(Laurasia$GenerationLength_d)["25%"],]$GL = "ShortLived"
Laurasia[Laurasia$GenerationLength_d >= quantile(Laurasia$GenerationLength_d)["75%"],]$GL = "LongLived"
Rumi = M[M$Order == "Ruminantia",]
Rumi$GL = "Middle"
Rumi[Rumi$GenerationLength_d <= quantile(Rumi$GenerationLength_d)["25%"],]$GL = "ShortLived"
Rumi[Rumi$GenerationLength_d >= quantile(Rumi$GenerationLength_d)["75%"],]$GL = "LongLived"

### Made MU test inside families
mu_df = rbind(Laurasia, Rumi)
mu_df = mu_df[mu_df$GL != 'Middle',]

out_m = data.frame()

for (ord in unique(mu_df$Order))
{
  for (gen in unique(mu_df$Gene))
  {
    subs = mu_df %>% filter(Order == ord, Gene == gen)
    test = wilcox.test(subs[subs$GL == 'ShortLived',]$CTSkew, subs[subs$GL == 'LongLived',]$CTSkew, paired = FALSE)
    med_sh = median(subs[subs$GL == 'ShortLived',]$CTSkew)
    med_ln = median(subs[subs$GL == 'LongLived',]$CTSkew)
    out_m = rbind(out_m, data.frame(ord, gen, test$p.value, med_sh, med_ln))
    
  }
}

names(out_m) = c('Order', 'Gene', 'p-value','median_short','median_long')

### Made MU test between families
out_m = data.frame()
for (gen in unique(mu_df$Gene)){
  subs = mu_df %>% filter(Gene == gen)
  test = wilcox.test(subs[subs$GL == 'ShortLived',]$CTSkew, subs[subs$GL == 'LongLived',]$CTSkew, paired = FALSE)
  med_sh = median(subs[subs$GL == 'ShortLived',]$CTSkew)
  med_ln = median(subs[subs$GL == 'LongLived',]$CTSkew)
  out_m = rbind(out_m, data.frame(gen, test$p.value, med_sh, med_ln)) 
}

names(out_m) = c('Gene', 'p-value','median_short','median_long')


####### naive multiple linear model, assigning numbers to order of genes. We can improve it substituting integers by real time or using more correct stat
FromGenesToNumbers = data.frame(c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB'),seq(1:10)); names(FromGenesToNumbers)=c('Gene','TSSS') # time spend single stranded
FromGenesToNumbers = cbind(FromGenesToNumbers, rep(c(0,1), each = 5))
names(FromGenesToNumbers)=c('Gene','TBSS','DummyHighTbss') # time spend single stranded
gl_and_tbss = merge(M,FromGenesToNumbers) 
gl_and_tbss$GT = factor(gl_and_tbss$GT, levels = c('short lived','long lived'))

gl_and_tbss = gl_and_tbss %>% mutate(DummyLonglived = as.numeric(GT)) %>% mutate(DummyLonglived = DummyLonglived - 1)

summary(lm(gl_and_tbss$CTSkew ~ log2(gl_and_tbss$GenerationLength_d) + gl_and_tbss$TBSS))
summary(lm(gl_and_tbss$CTSkew ~ log2(gl_and_tbss$GenerationLength_d) * gl_and_tbss$TBSS))

summary(lm(gl_and_tbss$CTSkew ~ gl_and_tbss$DummyLonglived + gl_and_tbss$TBSS))
summary(lm(gl_and_tbss$CTSkew ~ gl_and_tbss$DummyLonglived * gl_and_tbss$TBSS))

summary(lm(gl_and_tbss$CTSkew ~ gl_and_tbss$DummyLonglived + gl_and_tbss$DummyHighTbss))
summary(lm(gl_and_tbss$CTSkew ~ gl_and_tbss$DummyLonglived * gl_and_tbss$DummyHighTbss))

## complete lm for big orders

Laura = gl_and_tbss[gl_and_tbss$Order == 'Laurasiatheria',]

summary(lm(Laura$CTSkew ~ log2(Laura$GenerationLength_d) + Laura$TBSS))
summary(lm(Laura$CTSkew ~ Laura$DummyLonglived + Laura$TBSS))
summary(lm(Laura$CTSkew ~ Laura$DummyLonglived + Laura$DummyHighTbss))


Rumi = gl_and_tbss[gl_and_tbss$Order == 'Ruminantia',]


summary(lm(Rumi$CTSkew ~ log2(Rumi$GenerationLength_d) + Rumi$TBSS))
summary(lm(Rumi$CTSkew ~ Rumi$DummyLonglived + Rumi$TBSS))
summary(lm(Laura$CTSkew ~ Laura$DummyLonglived + Laura$DummyHighTbss))


Rodent = gl_and_tbss[gl_and_tbss$Order == 'Rodentia',]

summary(lm(Rodent$CTSkew ~ Rodent$DummyLonglived + Rodent$TBSS))
summary(lm(Rodent$CTSkew ~ Rodent$DummyLonglived * Rodent$TBSS))

Haplo = gl_and_tbss[gl_and_tbss$Order == 'Haplorrhini',] 

summary(lm(Haplo$CTSkew ~ Haplo$DummyLonglived + Haplo$TBSS))
summary(lm(Haplo$CTSkew ~ Haplo$DummyLonglived * Haplo$TBSS))


dev.off()

####beautiful boxplots###################

library("ggpubr")
pdf("../../Body/4Figures/WholeGenomeAnalyses.NoOverlap.AGSkew.R.02.pdf", height = 10, width = 20)
MM = M
MM$GT<- factor(MM$GT, levels = c("short lived", "long lived"))
ggboxplot(MM, "Gene", "CTSkew",
          fill = "GT", palette = c("#22AF1D", "#0B350B"), xlab="Genes", ylab="GA skew", title = "GA skew in long- vs shortlived mammals", legend.title = "Mammals", width = 0.7, notch = TRUE)
dev.off()
