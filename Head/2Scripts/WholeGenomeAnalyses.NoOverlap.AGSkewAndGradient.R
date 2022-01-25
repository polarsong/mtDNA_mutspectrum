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

########################################## GENOME WIDE SKEW FOR EACH SPECIES

AGG = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species,SynNuc$Class), FUN = sum)
names(AGG) = c('Species','Class','NeutralA','NeutralT','NeutralG','NeutralC')

## all six different skews
AGG$CTSkew = (AGG$NeutralC - AGG$NeutralT)/(AGG$NeutralC + AGG$NeutralT); summary(AGG$CTSkew) # GA on heavy
AGG$CGSkew = (AGG$NeutralC - AGG$NeutralG)/(AGG$NeutralC + AGG$NeutralG); summary(AGG$CGSkew) # 
AGG$CASkew = (AGG$NeutralC - AGG$NeutralA)/(AGG$NeutralC + AGG$NeutralA); summary(AGG$CASkew) # 
AGG$TGSkew = (AGG$NeutralT - AGG$NeutralG)/(AGG$NeutralT + AGG$NeutralG); summary(AGG$TGSkew) # 
AGG$TASkew = (AGG$NeutralT - AGG$NeutralA)/(AGG$NeutralT + AGG$NeutralA); summary(AGG$TASkew) # 
AGG$GASkew = (AGG$NeutralG - AGG$NeutralA)/(AGG$NeutralG + AGG$NeutralA); summary(AGG$GASkew) # 

AGG$TCSkew = (AGG$NeutralT - AGG$NeutralC)/(AGG$NeutralT + AGG$NeutralC); summary(AGG$CTSkew) # AG on heavy. Added it for fun, just to be sure, that it is opposite to CT (GA on heavy) 

GT = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

Mam = merge(AGG,GT, by ='Species')

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
Taxa = read.table("../../Body/2Derived/MammalianTaxonomy.txt", header = TRUE)
TaxaSpOrd=data.frame(Taxa$Species, Taxa$Order)
names(TaxaSpOrd) = c("Species", "Order")
M = merge(M, TaxaSpOrd)
Orders = unique(M$Order)
M$GT = NULL
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
  out = paste(i, IQR(M[M$Order == i,]$GenerationLength_d), sep = ' ')
  print(out)
}
#[1] "Feliformia 773.05315"
#[1] "Rodentia 478.73066123767"
#[1] "Ruminantia 799.3417708285"
#[1] "Laurasiatheria 2447.0747625"
#[1] "Haplorrhini 1150.105875"
#[1] "Microchiroptera 662.2593"
#[1] "Caniformia 832.217431938"
#[1] "Strepsirrhini 1371.624375"
#[1] "Cetacea 2664.5"
#[1] "Xenarthra 336.6637"
#[1] "Paucituberculata 0"
#[1] "Afrotheria 6993.4424166715"
#[1] "Megachiroptera 1120.91840714286"
#[1] "Diprotodontia 896.200516485402"
#[1] "Didelphimorphia 103.149975"
#[1] "Microbiotheria 0"
#[1] "Euarchontoglires 60.1993"
#[1] "Peramelemorphia 560.1198"
#[1] "Dasyuromorphia 291.8336"
#[1] "Notoryctemorphia 0"
#[1] "Monotremata 2808.2654"

Laurasia = M[M$Order == "Laurasiatheria",]
Laurasia$GL = "Middle"
Laurasia[Laurasia$GenerationLength_d <= quantile(Laurasia$GenerationLength_d)["25%"],]$GL = "ShortLived"
Laurasia[Laurasia$GenerationLength_d >= quantile(Laurasia$GenerationLength_d)["75%"],]$GL = "LongLived"
Afro = M[M$Order == "Afrotheria",]
Afro$GL = "Middle"
Afro[Afro$GenerationLength_d <= quantile(Afro$GenerationLength_d)["25%"],]$GL = "ShortLived"
Afro[Afro$GenerationLength_d >= quantile(Afro$GenerationLength_d)["75%"],]$GL = "LongLived"
Ceta =  M[M$Order == "Cetacea",]
Ceta$GL = "Middle"
Ceta[Ceta$GenerationLength_d <= quantile(Ceta$GenerationLength_d)["25%"],]$GL = "ShortLived"
Ceta[Ceta$GenerationLength_d >= quantile(Ceta$GenerationLength_d)["75%"],]$GL = "LongLived"
Hapl =  M[M$Order == "Haplorrhini",]
Hapl$GL = "Middle"
Hapl[Hapl$GenerationLength_d <= quantile(Hapl$GenerationLength_d)["25%"],]$GL = "ShortLived"
Hapl[Hapl$GenerationLength_d >= quantile(Hapl$GenerationLength_d)["75%"],]$GL = "LongLived"
Strep = M[M$Order == "Strepsirrhini",]
Strep$GL = "Middle"
Strep[Strep$GenerationLength_d <= quantile(Strep$GenerationLength_d)["25%"],]$GL = "ShortLived"
Strep[Strep$GenerationLength_d >= quantile(Strep$GenerationLength_d)["75%"],]$GL = "LongLived"





####### naive multiple linear model, assigning numbers to order of genes. We can improve it substituting integers by real time or using more correct stat
FromGenesToNumbers = data.frame(c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB'),seq(1:10)); names(FromGenesToNumbers)=c('Gene','TSSS') # time spend single stranded
M=merge(M,FromGenesToNumbers) 

A<-lm(M$CTSkew ~ log2(M$GenerationLength_d) +  M$TSSS); summary(A)
#  (Intercept)                -0.8477990  0.0283319  -29.92   <2e-16 ***
#  log2(M$GenerationLength_d)  0.0722798  0.0025027   28.88   <2e-16 ***
#  M$TSSS                      0.0387429  0.0009798   39.54   <2e-16 ***

A<-lm(M$CTSkew ~ log2(M$GenerationLength_d) +  scale(M$TSSS)); summary(A) # unique(scale(M$TSSS))
# (Intercept)                -0.634713   0.027815  -22.82   <2e-16 ***
#  log2(M$GenerationLength_d)  0.072280   0.002503   28.88   <2e-16 ***
#  scale(M$TSSS)               0.111289   0.002815   39.54   <2e-16 ***

A<-lm(M$CTSkew ~ scale(M$GenerationLength_d) +  scale(M$TSSS)); summary(A)
# (Intercept)                 0.164476   0.002865   57.41   <2e-16 ***
#  scale(M$GenerationLength_d) 0.068882   0.002865   24.04   <2e-16 ***
#  scale(M$TSSS)               0.111289   0.002865   38.84   <2e-16 ***

A<-lm(M$CTSkew ~ scale(M$GenerationLength_d)*scale(M$TSSS)); summary(A)  # interaction is not significant for CT.
# (Intercept)                                0.164476   0.002865  57.412   <2e-16 ***
#  scale(M$GenerationLength_d)                0.068882   0.002865  24.042   <2e-16 ***
#  scale(M$TSSS)                              0.111289   0.002865  38.844   <2e-16 ***
#  scale(M$GenerationLength_d):scale(M$TSSS) -0.002585   0.002865  -0.902    0.367  

dev.off()

####beautiful boxplots###################

library("ggpubr")
pdf("../../Body/4Figures/WholeGenomeAnalyses.NoOverlap.AGSkew.R.02.pdf", height = 10, width = 20)
MM = M
MM$GT<- factor(MM$GT, levels = c("short lived", "long lived"))
ggboxplot(MM, "Gene", "CTSkew",
          fill = "GT", palette = c("#22AF1D", "#0B350B"), xlab="Genes", ylab="GA skew", title = "GA skew in long- vs shortlived mammals", legend.title = "Mammals", width = 0.7, notch = TRUE)
dev.off()
