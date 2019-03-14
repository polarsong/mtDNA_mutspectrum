# In GTEx we expect more strong selection (as compared to cancer data) - (i) play with annotation (ii) play with VAF (the rarier - the more neutral)
# We can check Samuels tissue-specific paper

rm(list=ls(all=TRUE))

pdf("../../Body/4Figures/HealthySomaticHumanMutations.GTEx.R.01.pdf")
par(mfrow=c(2,2))

############ read human ref seq (note that it is hg19 - not cambridge reference!):
RefSeq = read.table("../../Body/1Raw/chrM_refAllele.txt", header = FALSE, sep = '')
names(RefSeq)=c('Position','AncestralAllele')

############ read GTEx mutations and merge them with Ref Seq:
SomMut = read.table("../../Body/1Raw/TissueSpecificMtDNAMutationsGTEx.txt", header = TRUE, sep = '\t')
SomMut$Position = gsub('_(.*)','',SomMut$Mutation)  # gsub("_(.*)",'',READS$patient)
SomMut$DerivedAllele = gsub('(.*)_','',SomMut$Mutation)  # gsub("_(.*)",'',READS$patient)

Som = merge(SomMut,RefSeq, by = 'Position')

Som = Som[Som$DerivedAllele != Som$AncestralAllele,] # correct - in all lines Derived != Ancestral
Som$Substitution = paste(Som$AncestralAllele,Som$DerivedAllele, sep = '_')
table(Som$Substitution)
# A_C A_G A_T C_A C_G C_T G_A G_C G_T T_A T_C T_G 
# 43 473  49  49   6 336 728  63  57  42 693  26

## I don't need common varaints because of potential selection => it is better to go to more rare
hist(Som$AF, breaks = 100) 
summary(Som$AF) # median = 0.0437
# Som = Som[Som$AF <= median(Som$AF),]
# Som = Som[Som$AF <= 0.0687,]
Som = Som[Som$AF <= 0.0437,]

############ read GTEx expression data:
Expr = read.table("../../Body/1Raw/TissueSpecificMtDNACoverageGTEx.txt", header = TRUE, sep = '\t')

# merge by SRRID (lost 2565 - 2523 observations, but got expression level (percentMito) for each sample and alternative name of the Tissue (more simple))
nrow(Som) # 2565
Som = merge(Som,Expr, by.x = 'SRR', by.y = 'SRRID')
nrow(Som) # 2523
VecOfGTExTissues = unique(Som$Tissue)
VecOfGTExTissues
# Skin            Esophagus       Colon           Blood Vessel    Breast          Blood           Liver           Muscle          Adrenal Gland   Vagina         
# Stomach         Pancreas        Nerve           Uterus          Adipose Tissue  Prostate        Spleen          Testis          Kidney          Lung           
# Ovary           Small Intestine Salivary Gland  Brain           Thyroid         Pituitary       Heart 

########### analyses: derive for each Tissue mean percentMito (approximation of the level of metabolism) and MutSpec 

### derive percentMito:
Agg1 = aggregate(Som$percentMito, by = list(Som$Tissue), FUN = mean);
names(Agg1)=c('Tissue','percentMito')
Agg1 = Agg1[order(Agg1$percentMito),]

### derive MutSpec:
Som$T_C = 0; Som$G_A = 0;
Som$C_T = 0; Som$A_G = 0;
Som$Tv = 0; Som$Ts = 0;
for (i in 1:nrow(Som))
{
  if (Som$Substitution[i] == 'T_C') {Som$T_C[i] = 1;}
  if (Som$Substitution[i] == 'C_T') {Som$C_T[i] = 1;}
  if (Som$Substitution[i] == 'G_A') {Som$G_A[i] = 1;}
  if (Som$Substitution[i] == 'A_G') {Som$A_G[i] = 1;}
  if (Som$Substitution[i] == 'T_C' | Som$Substitution[i] == 'C_T' | Som$Substitution[i] == 'G_A' | Som$Substitution[i] == 'A_G') {Som$Ts[i] = 1;}
  if (Som$Substitution[i] == 'A_T' | Som$Substitution[i] == 'T_A' | Som$Substitution[i] == 'G_C' | Som$Substitution[i] == 'C_G' | Som$Substitution[i] == 'A_C' | Som$Substitution[i] == 'C_A'  | Som$Substitution[i] == 'G_T' | Som$Substitution[i] == 'T_G') {Som$Tv[i] = 1;}
}
Agg2 = aggregate(list(Som$T_C,Som$C_T,Som$G_A,Som$A_G,Som$Ts, Som$Tv), by = list(Som$Tissue), FUN = mean);
names(Agg2) = c('Tissue','T_C','C_T','G_A','A_G','Ts','Tv')

### derive total number of all mutations:
Som$Mut = 1
Agg3 = aggregate(Som$Mut, by = list(Som$Tissue), FUN = sum);
names(Agg3) = c('Tissue','TotalMut')

### merge by tissue three Agg dataframes 

Agg = merge(Agg1,Agg2, by = 'Tissue'); Agg = merge(Agg,Agg3, by = 'Tissue')

### compare percentMito with MutSpec (a trend that T>C is higher in tissues with high percentMito)

Agg = Agg[order(Agg$percentMito),]
Agg$TsTv = Agg$Ts/Agg$Tv
Agg$TC_GA = Agg$T_C/Agg$G_A
cor.test(Agg$percentMito,Agg$T_C, method = 'spearman');
cor.test(Agg$percentMito,Agg$T_C, method = 'spearman', alternative = 'greater'); # p = 0.08 
cor.test(Agg[Agg$TotalMut >= 10,]$percentMito,Agg[Agg$TotalMut >= 10,]$T_C, method = 'spearman', alternative = 'greater');  # 0.048
cor.test(Agg[Agg$TotalMut >= 10,]$percentMito,Agg[Agg$TotalMut >= 10,]$TC_GA, method = 'spearman', alternative = 'greater'); # 0.078
boxplot(Agg[Agg$percentMito < median(Agg$percentMito),]$T_C,Agg[Agg$percentMito >= median(Agg$percentMito),]$T_C)
wilcox.test(Agg[Agg$percentMito < median(Agg$percentMito),]$T_C,Agg[Agg$percentMito >= median(Agg$percentMito),]$T_C, alternative = 'less')
boxplot(Agg[Agg$TotalMut >= 10 & Agg$percentMito < median(Agg$percentMito),]$T_C,Agg[Agg$TotalMut >= 10 & Agg$percentMito >= median(Agg$percentMito),]$T_C)
wilcox.test(Agg[Agg$percentMito < median(Agg$percentMito),]$T_C,Agg[Agg$percentMito >= median(Agg$percentMito),]$T_C, alternative = 'less')

cor.test(Agg$percentMito,Agg$G_A, method = 'spearman')
cor.test(Agg$percentMito,Agg$Tv, method = 'spearman')
cor.test(Agg$percentMito,Agg$TsTv, method = 'spearman')
cor.test(Agg$percentMito,Agg$Ts, method = 'spearman')
cor.test(Agg$percentMito,Agg$TotalMut, method = 'spearman')

plot(Agg$percentMito,Agg$T_C, pch = ''); text(Agg$percentMito,Agg$T_C,Agg$Tissue)
plot(Agg[Agg$TotalMut >=10,]$percentMito,Agg[Agg$TotalMut >=10,]$T_C, pch = ''); text(Agg[Agg$TotalMut >=10,]$percentMito,Agg[Agg$TotalMut >=10,]$T_C,Agg[Agg$TotalMut >=10,]$Tissue)
plot(Agg[Agg$TotalMut > 30,]$percentMito,Agg[Agg$TotalMut > 30,]$T_C);

### lm

Som$Position = as.numeric(as.character(Som$Position))
a<-lm(Som$T_C ~ Som$percentMito + Som$AF); summary(a)
a<-lm(Som$G_A ~ Som$percentMito + Som$AF); summary(a)
a<-lm(Som$C_T ~ Som$percentMito + Som$AF); summary(a)
a<-lm(Som$A_G ~ Som$percentMito + Som$AF); summary(a) # negative with percentMito 
a<-lm(Som$Ts ~ Som$percentMito + Som$AF + Som$Position); summary(a)
a<-lm(Som$Ts ~ Som$percentMito + Som$AF); summary(a) # Ts are more rare in high percentMito (due to A>G)
a<-lm(Som$Tv ~ Som$percentMito + Som$AF); summary(a) # Tv are more often in high percentMito

###### derive cell lifespan: https://docs.google.com/document/d/1UECub1DIdmuXwPqDLK8WRcZ6uIjEQiVmHXN1_XNVvUU/edit?usp=sharing 
## Gladyshev: https://www.nature.com/articles/npjamd201614
# some changes in order to fit Tissue names from Som:
# Muscle <= sceletal muscle
# Heart <= heart Muscle
# Brain <= Neuron (neocortex)
# Thyroid <= Thyroid gland
# Skin <= Keratinocytes (skin epidermis)
VecOfTissuesGladyshev = c('Adipose Tissue', 'Adrenal Gland', 'Bone Marrow', 'Monocytes', 'Colon',  'Endometrium', 'Esophagus', 'Heart',  'Skin', 'Kidney', 'Liver', 'Lung', 'Brain', 'Osteoblasts (bone)', 'Rectum', 'Salivary Gland', 'Muscle', 'Smooth muscle', 'Spleen',  'Thyroid',  'Bladder')  
VecOfLifespansGladyshev = c(2448,455,3.2,2,3.5,13,10,25300,64,270,327, 200, 32850, 8.3,  3.5, 60, 5510, 67.5, 7.8, 3180, 49) 
First = data.frame(VecOfTissuesGladyshev,VecOfLifespansGladyshev); names(First) = c('Tissue','TurnOverRate')

intersect(VecOfGTExTissues,VecOfTissuesGladyshev)
# Skin Esophagus"      "Colon"          "Liver"          "Muscle"         "Adrenal Gland"  "Adipose Tissue" "Spleen"         "Kidney"         "Lung"          
# "Salivary Gland" "Brain"          "Thyroid"        "Heart"   
setdiff(VecOfGTExTissues,VecOfTissuesGladyshev)
# "Blood Vessel"    "Breast"          "Blood"           "Vagina"          "Stomach"         "Pancreas"        "Nerve"           "Uterus"         
# "Prostate"        "Testis"          "Ovary"           "Small Intestine" "Pituitary"      

### get more from Tomacetti and others: 
VecOfTissuesTomacettiAndOthers = c("Breast", "Stomach", "Pancreas", "Uterus", 'Prostate','Testis', 'Ovary', 'Small Intestine')
VecOfLifespansTomacettiAndOthers = c(84.5, 5.5, 360, 4, 120.0, 60, 11000, 3)
Second = data.frame(VecOfTissuesTomacettiAndOthers,VecOfLifespansTomacettiAndOthers); names(Second) = c('Tissue','TurnOverRate')
All = rbind(First,Second)

setdiff(VecOfGTExTissues,All$Tissue)
# "Blood Vessel" "Blood"        "Vagina"       "Nerve"        "Pituitary"    - to find it

####### merge Som with All
nrow(Som) # 2523
Som = merge(Som,All, by = 'Tissue')
nrow(Som) # 1608 = lost a lot due to these several tissues

Agg4 = aggregate(list(Som$T_C,Som$G_A,Som$C_T,Som$A_G,Som$Ts,Som$Tv,Som$percentMito), by = list(Som$Tissue,Som$TurnOverRate), FUN = mean); 
names(Agg4)=c('Tissue','TurnOverRate','T_C','G_A','C_T','A_G','Ts','Tv','percentMito')
Agg4$TC_GA = Agg4$T_C/Agg4$G_A
Agg4$TC_Tv = Agg4$T_C/Agg4$Tv

cor.test(Agg4$TurnOverRate,Agg4$percentMito, method = 'spearman') # nothing
cor.test(Agg4$TurnOverRate,Agg4$T_C, method = 'spearman')   # nothing
cor.test(Agg4$TurnOverRate,Agg4$G_A, method = 'spearman')   # negative trend
cor.test(Agg4$TurnOverRate,Agg4$TC_GA, method = 'spearman') # nothing
cor.test(Agg4$TurnOverRate,Agg4$TC_Tv, method = 'spearman')
plot(Agg4$TurnOverRate,Agg4$G_A,pch = ''); text(Agg4$TurnOverRate,Agg4$G_A,Agg4$Tissue)

Agg = merge(Agg[,c(1,9)],Agg4, by = 'Tissue')
cor.test(Agg$TurnOverRate,Agg$TC_GA, method = 'spearman') # positive trend
cor.test(Agg[Agg$TotalMut >= 10,]$TurnOverRate,Agg[Agg$TotalMut >= 10,]$TC_GA, method = 'spearman') # positive - and even better!
# if I take rare and get T>C/G>A I have good positive correlation
plot(Agg[Agg$TotalMut >= 10,]$TurnOverRate,Agg[Agg$TotalMut >= 10,]$TC_GA,pch = ''); text(Agg[Agg$TotalMut >= 10,]$TurnOverRate,Agg[Agg$TotalMut >= 10,]$TC_GA,Agg[Agg$TotalMut >= 10,]$Tissue)
par(mfrow=c(1,1))
plot(log2(Agg[Agg$TotalMut >= 10,]$TurnOverRate),Agg[Agg$TotalMut >= 10,]$TC_GA,pch = ''); text(log2(Agg[Agg$TotalMut >= 10,]$TurnOverRate),Agg[Agg$TotalMut >= 10,]$TC_GA,Agg[Agg$TotalMut >= 10,]$Tissue)

dev.off()
