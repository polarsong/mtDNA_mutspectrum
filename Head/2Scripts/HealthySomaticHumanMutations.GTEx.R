# In GTEx we expect more strong selection (as compared to cancer data) - (i) play with annotation (ii) play with VAF (the rarier - the more neutral) (iii) take care of tissue-specificity in each patient
# We can check Samuels tissue-specific paper

rm(list=ls(all=TRUE))

pdf("../../Body/4Figures/HealthySomaticHumanMutations.GTEx.R.01.pdf")

############ A: read human ref seq (note that it is hg19 - not cambridge reference!):
RefSeq = read.table("../../Body/1Raw/chrM_refAllele.txt", header = FALSE, sep = '')
names(RefSeq)=c('Position','AncestralAllele')

############ B: read GTEx mutations and merge them with Ref Seq:
SomMut = read.table("../../Body/1Raw/TissueSpecificMtDNAMutationsGTEx.txt", header = TRUE, sep = '\t')
SomMut$Position = gsub('_(.*)','',SomMut$Mutation)  # gsub("_(.*)",'',READS$patient)
SomMut$DerivedAllele = gsub('(.*)_','',SomMut$Mutation)  # gsub("_(.*)",'',READS$patient)

Som = merge(SomMut,RefSeq, by = 'Position')

Som = Som[Som$DerivedAllele != Som$AncestralAllele,] # correct - in all lines Derived != Ancestral
Som$Substitution = paste(Som$AncestralAllele,Som$DerivedAllele, sep = '_')
table(Som$Substitution)
# A_C A_G A_T C_A C_G C_T G_A G_C G_T T_A T_C T_G 
# 43 473  49  49   6 336 728  63  57  42 693  26

write.table(Som, file = "../../Body/2Derived/GTExRef.txt", quote = FALSE, row.names = FALSE)

### It is reasonable to remove variants, shared between at least two tissues in the same individuum (because in this case we don't know the tisse of origin). 
Som = Som[order(Som$subject, Som$Mutation),]


## I don't need common varaints because of potential selection => it is better to go to more rare
hist(Som$AF, breaks = 100) 
summary(Som$AF) # median = 0.0437
# Som = Som[Som$AF <= median(Som$AF),]
# Som = Som[Som$AF <= 0.0687,]
# Som = Som[Som$AF <= 0.0437,]
# Som = Som[Som$AF <= 0.1,]

Som$Position = as.numeric(as.character(Som$Position))
# Som = Som[Som$Position > 5000,]

############ C: read GTEx expression data:
Expr = read.table("../../Body/1Raw/TissueSpecificMtDNACoverageGTEx.txt", header = TRUE, sep = '\t')
ExprTissues = data.frame(table(Expr$Tissue))
ExprTissues = ExprTissues[order(ExprTissues$Var1),]

# merge by SRRID (lost 2565 - 2523 observations, but got expression level (percentMito) for each sample and alternative name of the Tissue (more simple))
nrow(Som) # 2565
Som = merge(Som,Expr, by.x = 'SRR', by.y = 'SRRID')
nrow(Som) # 2523
VecOfGTExTissues = unique(Som$Tissue)
VecOfGTExTissues
# Skin            Esophagus       Colon           Blood Vessel    Breast          Blood           Liver           Muscle          Adrenal Gland   Vagina         
# Stomach         Pancreas        Nerve           Uterus          Adipose Tissue  Prostate        Spleen          Testis          Kidney          Lung           
# Ovary           Small Intestine Salivary Gland  Brain           Thyroid         Pituitary       Heart 

############ D: derive MutSpec
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

########### E: Cell lifespans from the table here: https://docs.google.com/document/d/1UECub1DIdmuXwPqDLK8WRcZ6uIjEQiVmHXN1_XNVvUU/edit?usp=sharing 
VecOfTissues = c('Adipose Tissue','Adrenal Gland','Bladder','Blood','Blood Vessel','Bone Marrow','Brain','Breast','Cervix Uteri','Colon','Esophagus','Heart','Kidney','Liver','Lung','Muscle','Nerve','Ovary','Pancreas','Prostate','Salivary Gland','Skin','Small Intestine','Spleen','Stomach','Testis','Thyroid','Uterus','Vagina')
VecOfTurnOvers = c(  2448,  455,  49,  30,  67.5,  3.2,  24637.5,  65.75,  5.7,  4.25,  10.5,  25300,  270,  363.5,  126,  5510,  24637.5,  30000,  315,  120,  60,  64,  7.05,  7.8,  1.4,  63.5,  3679.5,  13,  3.9)
All = data.frame(VecOfTissues,VecOfTurnOvers); names(All) = c('Tissue','TurnOverRate')
Som = merge(Som,All, by = 'Tissue') 

setdiff(VecOfGTExTissues,All$Tissue) # Pituitary


# derive number of total somatic mutations per each individuum
Som$Number = 1
agg = aggregate(Som$Number, by = list(Som$subject), FUN = sum); names(agg) = c('subject','NumOfMutPerIndiv')
Som = merge(Som,agg, by = 'subject')

### COverage is very important; The higher the coverage the higher probability to see G>A. Why?
## logistic regression with a set of confounders:
# nDonorTissue = number of tissues from a given donor
# COV = coverage if a given region in a donor
# AF = heteroplasmy level in a given sample (minimum is 3%)
# NumOfMutPerIndiv = number of mutations per a given individuum
# Position = position => sibstitute it by the time spent single stranded (from 6000 (COX1) to 16000(CYTB) time being single stranded is increasing, the rest is unknown (for me) - may be we can approximate it somehow)
# position is singificant, but there are two explanations: time being single-stranded and selection (there are many mutation in D loop???)
# repeat it with synonymous only!!!!!!!!!!! KG help????

# backward multiple logistic regression for G>A
a<-glm(Som$G_A ~ Som$TurnOverRate + Som$COV + Som$NumOfMutPerIndiv + Som$AF + Som$Position + Som$nDonorTissue + Som$percentMito + Som$uniqueMappedReads + Som$uniqueMappedReadsPercent, family = binomial()); summary(a)
a<-glm(Som$G_A ~ Som$TurnOverRate + Som$COV + Som$NumOfMutPerIndiv + Som$Position + Som$nDonorTissue + Som$percentMito + Som$uniqueMappedReads + Som$uniqueMappedReadsPercent, family = binomial()); summary(a)
a<-glm(Som$G_A ~ Som$TurnOverRate + Som$COV + Som$NumOfMutPerIndiv + Som$Position + Som$percentMito + Som$uniqueMappedReads + Som$uniqueMappedReadsPercent, family = binomial()); summary(a)
a<-glm(Som$G_A ~ Som$TurnOverRate + Som$COV + Som$NumOfMutPerIndiv + Som$Position + Som$uniqueMappedReads + Som$uniqueMappedReadsPercent, family = binomial()); summary(a)
a<-glm(Som$G_A ~ Som$TurnOverRate + Som$COV + Som$NumOfMutPerIndiv + Som$Position + Som$uniqueMappedReads, family = binomial()); summary(a)
a<-glm(Som$G_A ~ Som$TurnOverRate + Som$COV + Som$NumOfMutPerIndiv + Som$Position, family = binomial()); summary(a)
a<-glm(Som$G_A ~ Som$TurnOverRate + Som$COV + Som$Position, family = binomial()); summary(a) # final step
# I don't understand by heart effect of position and coverage

MajorArc = Som[Som$Position > 6000 & Som$Position < 16000,]
a<-glm(MajorArc$G_A ~ MajorArc$TurnOverRate + MajorArc$COV + MajorArc$NumOfMutPerIndiv + MajorArc$AF + MajorArc$Position, family = binomial())
summary(a)

# backward multiple logistic regression for T>C
a<-glm(Som$T_C ~ Som$TurnOverRate + Som$COV + Som$NumOfMutPerIndiv + Som$AF + Som$Position + Som$nDonorTissue + Som$percentMito + Som$uniqueMappedReads + Som$uniqueMappedReadsPercent, family = binomial()); summary(a)
a<-glm(Som$T_C ~ Som$TurnOverRate + Som$NumOfMutPerIndiv + Som$AF + Som$Position + Som$nDonorTissue + Som$percentMito + Som$uniqueMappedReads + Som$uniqueMappedReadsPercent, family = binomial()); summary(a)
a<-glm(Som$T_C ~ Som$TurnOverRate + Som$NumOfMutPerIndiv + Som$AF + Som$nDonorTissue + Som$percentMito + Som$uniqueMappedReads + Som$uniqueMappedReadsPercent, family = binomial()); summary(a)
a<-glm(Som$T_C ~ Som$TurnOverRate + Som$NumOfMutPerIndiv + Som$AF + Som$nDonorTissue + Som$uniqueMappedReads + Som$uniqueMappedReadsPercent, family = binomial()); summary(a)
a<-glm(Som$T_C ~ Som$TurnOverRate + Som$NumOfMutPerIndiv + Som$AF + Som$uniqueMappedReads + Som$uniqueMappedReadsPercent, family = binomial()); summary(a)
a<-glm(Som$T_C ~ Som$TurnOverRate + Som$NumOfMutPerIndiv + Som$AF + Som$uniqueMappedReads, family = binomial()); summary(a)
a<-glm(Som$T_C ~ Som$TurnOverRate + Som$NumOfMutPerIndiv + Som$AF, family = binomial()); summary(a)
a<-glm(Som$T_C ~ Som$TurnOverRate + Som$AF, family = binomial()); summary(a) # final step
# weak and NEGATIVE (EXPECT POSITIVE) correlation with TurnOver and correlation with AF ~ heteroplasmy level: the higher the heteroplasmy the more chances that mutation is T>C (T>C are "old" alleles)

## multiple model: TurnOver as a function of MutSpec:
a<-lm(Som$TurnOverRate ~  Som$G_A + Som$A_G + Som$T_C + Som$C_T + Som$AF + Som$COV + Som$NumOfMutPerIndiv + Som$Position + Som$nDonorTissue + Som$percentMito + Som$uniqueMappedReads + Som$uniqueMappedReadsPercent); summary(a)
a<-lm(Som$TurnOverRate ~  Som$G_A + Som$A_G + Som$T_C + Som$C_T + Som$AF + Som$COV + Som$NumOfMutPerIndiv + Som$Position + Som$nDonorTissue + Som$percentMito + Som$uniqueMappedReadsPercent); summary(a)
a<-lm(Som$TurnOverRate ~  Som$G_A + Som$A_G + Som$T_C + Som$C_T + Som$AF + Som$COV + Som$NumOfMutPerIndiv + Som$nDonorTissue + Som$percentMito + Som$uniqueMappedReadsPercent); summary(a)
a<-lm(Som$TurnOverRate ~  Som$G_A + Som$A_G + Som$T_C + Som$AF + Som$COV + Som$NumOfMutPerIndiv + Som$nDonorTissue + Som$percentMito + Som$uniqueMappedReadsPercent); summary(a)
a<-lm(Som$TurnOverRate ~  Som$G_A + Som$A_G + Som$T_C + Som$COV + Som$NumOfMutPerIndiv + Som$nDonorTissue + Som$percentMito + Som$uniqueMappedReadsPercent); summary(a) # final
a<-lm(scale(Som$TurnOverRate) ~ 0 +  scale(Som$G_A) + scale(Som$A_G) + scale(Som$T_C) + scale(Som$COV) + scale(Som$NumOfMutPerIndiv) + scale(Som$nDonorTissue) + scale(Som$percentMito) + scale(Som$uniqueMappedReadsPercent)); summary(a) # final
# all transitions are decreasing with cell longevity!???
# PercentMito (~expression level of mtDNA) positively correlate with cell longevity => neurons and muscles are active! 
#  scale(Som$G_A)                      -0.14796    0.02432  -6.084 1.36e-09 ***
#  scale(Som$A_G)                      -0.06438    0.02316  -2.779 0.005486 ** 
#  scale(Som$T_C)                      -0.11350    0.02397  -4.736 2.31e-06 ***
#  scale(Som$COV)                      -0.05521    0.02027  -2.723 0.006509 ** 
#  scale(Som$NumOfMutPerIndiv)          0.10628    0.02026   5.246 1.69e-07 ***
#  scale(Som$nDonorTissue)             -0.07082    0.01964  -3.605 0.000318 ***
#  scale(Som$percentMito)               0.16913    0.02102   8.047 1.30e-15 ***
#  scale(Som$uniqueMappedReadsPercent) -0.07080    0.02123  -3.336 0.000864 ***

# continue with more stringent p value threshold:
a<-lm(scale(Som$TurnOverRate) ~ 0 +  scale(Som$G_A) + scale(Som$A_G) + scale(Som$T_C) + scale(Som$NumOfMutPerIndiv) + scale(Som$nDonorTissue) + scale(Som$percentMito) + scale(Som$uniqueMappedReadsPercent)); summary(a) # final
a<-lm(scale(Som$TurnOverRate) ~ 0 +  scale(Som$G_A) + scale(Som$T_C) + scale(Som$NumOfMutPerIndiv) + scale(Som$nDonorTissue) + scale(Som$percentMito) + scale(Som$uniqueMappedReadsPercent)); summary(a) # final
a<-lm(scale(Som$TurnOverRate) ~ 0 +  scale(Som$G_A) + scale(Som$T_C) + scale(Som$NumOfMutPerIndiv) + scale(Som$nDonorTissue) + scale(Som$percentMito)); summary(a) # final
a<-lm(scale(Som$TurnOverRate) ~ 0 +  scale(Som$G_A) + scale(Som$T_C) + scale(Som$NumOfMutPerIndiv) + scale(Som$percentMito)); summary(a) # final
# scale(Som$G_A)              -0.13273    0.02132  -6.224 5.67e-10 ***
# scale(Som$T_C)              -0.09113    0.02130  -4.278 1.96e-05 ***
# scale(Som$NumOfMutPerIndiv)  0.12505    0.02005   6.238 5.19e-10 ***
# scale(Som$percentMito)       0.13891    0.02000   6.947 4.77e-12 ***

########### G: analyses: derive for each Tissue mean percentMito (approximation of the level of metabolism), MutSpec and merge with 

### G1: percentMito versus TurnOver
Agg1 = aggregate(Som$percentMito, by = list(Som$Tissue), FUN = mean);
names(Agg1)=c('Tissue','percentMito')
Agg1 = Agg1[order(Agg1$percentMito),]

### G2: derive total number of all mutations:
Som$Mut = 1
Agg2 = aggregate(Som$Mut, by = list(Som$Tissue), FUN = sum);
names(Agg2) = c('Tissue','TotalMut')

### percentMito versus MutSpec: nothing
summary(Som$AF) # 0.03000 0.03520 0.04370 0.07916 0.06893 0.99120 
# FOR ALL:
Agg3 = aggregate(list(Som$T_C,Som$C_T,Som$G_A,Som$A_G,Som$Ts, Som$Tv), by = list(Som$Tissue,Som$TurnOverRate), FUN = mean);
names(Agg3) = c('Tissue','TurnOverRate','T_C','C_T','G_A','A_G','Ts','Tv')

AGG = merge(Agg1,Agg2, by = 'Tissue'); AGG = merge(AGG,Agg3, by = 'Tissue');  
AGG$TC_GA =AGG$T_C / AGG$G_A
AGG$TsTv =AGG$Ts/ AGG$Tv

######### analyses Cor-tests:

## three main rate-related properties: nothing
cor.test(AGG$percentMito,AGG$TotalMut, method = 'spearman')
cor.test(AGG$percentMito,AGG$TurnOverRate, method = 'spearman')
cor.test(AGG$TotalMut,AGG$TurnOverRate, method = 'spearman')

## mut spec and TurnOver (when I decrease VAF - effect is getting weaker)

cor.test(AGG$T_C,AGG$TurnOverRate, method = 'spearman')   # 0.8826
cor.test(AGG$C_T,AGG$TurnOverRate, method = 'spearman')   # a bit positive???
cor.test(AGG$G_A,AGG$TurnOverRate, method = 'spearman')   # -0.5130049, 0.007362
cor.test(AGG$A_G,AGG$TurnOverRate, method = 'spearman')   # a bit positive???
cor.test(AGG$TsTv,AGG$TurnOverRate, method = 'spearman')  # nothing.. TsTv doesn't work, but TC/GA works well
cor.test(AGG$TC_GA,AGG$TurnOverRate, method = 'spearman') # 0.5165042, 0.006904 - if we decrease FAV - correlation is robust more or less.

nrow(AGG[AGG$TotalMut >= 10,]) # 25
cor.test(AGG[AGG$TotalMut >= 10,]$T_C,AGG[AGG$TotalMut >= 10,]$TurnOverRate, method = 'spearman') # 00.9723 PAPER
cor.test(AGG[AGG$TotalMut >= 10,]$TC_GA,AGG[AGG$TotalMut >= 10,]$TurnOverRate, method = 'spearman') # 0.5791803, 0.002415;  PAPER
cor.test(AGG[AGG$TotalMut >= 10,]$G_A,AGG[AGG$TotalMut >= 10,]$TurnOverRate, method = 'spearman') # -0.532615, 0.006125 PAPER

plot(AGG[AGG$TotalMut >= 10,]$TC_GA,log2(AGG[AGG$TotalMut >= 10,]$TurnOverRate), pch = 16, col ='grey'); text(AGG[AGG$TotalMut >= 10,]$TC_GA,log2(AGG[AGG$TotalMut >= 10,]$TurnOverRate),AGG[AGG$TotalMut >= 10,]$Tissue);   # dev.off()
plot(AGG[AGG$TotalMut >= 10,]$G_A,log2(AGG[AGG$TotalMut >= 10,]$TurnOverRate), pch = 16, col ='grey'); text(AGG[AGG$TotalMut >= 10,]$G_A,log2(AGG[AGG$TotalMut >= 10,]$TurnOverRate),AGG[AGG$TotalMut >= 10,]$Tissue);   # dev.off()

######## analyses lm:

a<-lm(Som$G_A ~ Som$TurnOverRate + Som$AF); summary(a)
a<-lm(Som$T_C ~ Som$TurnOverRate + Som$AF); summary(a)
a<-lm(Som$TurnOverRate ~ Som$AF); summary(a) # negative a bit -> in slowly dividing tissues (with high number of days) AF are low.

######## distribution along the genome T>C and G>A

par(mfrow=c(2,1))
hist(Som[Som$G_A == 1,]$Position, breaks = 100) # decreasing a bit
hist(Som[Som$T_C == 1,]$Position, breaks = 100) # goes 
hist(Som$Position, breaks = 100) # goes 

dev.off()
