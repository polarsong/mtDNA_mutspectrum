# our table: one line - one mutation
# (мутации мтДНК) ~ гипоксия (Buffa score AND/OR скорость деления клеток) + VAF (время появления) + ткань (тканеспецифичное что-то)

#################################################
##### 1: DERIVE TABLE
#################################################

rm(list=ls(all=TRUE))
library(dplyr)
library(tidyr)
library(ggplot2)

Mut = read.table("../../Body/1Raw/mtDNA_snv_Oct2016.txt", head = TRUE, sep = '\t')  
Decode = read.table("../../Body/1Raw/PancancerInfoMetadata.txt", head = TRUE, sep = '\t')  
Decode = Decode[,c(2,4)]
Mut = merge(Mut,Decode, by.x = 'sample', by.y = 'submitter_donor_id', all.x = TRUE)
Hyp = read.table("../../Body/1Raw/hypoxicCancers.txt", head = TRUE, sep = '\t')  
HypMut = merge(Mut,Hyp, by = 'aliquot_id')

HypMut$Subst = paste(HypMut$ref,HypMut$var,sep = '>')
names(HypMut)
HypMut = select(HypMut,aliquot_id,sample,position,Subst,tissue,Tier2,tumor_var_freq,hypoxia_score_winter,hypoxia_score_ragnum,hypoxia_score_buffa)
HypMut$tumor_var_freq = as.numeric(gsub('%','',HypMut$tumor_var_freq))
nrow(HypMut) # 3110 mutations with known hypoxia
VecOfSamples = unique(HypMut$sample); length(VecOfSamples) # 828 samples

CancerTissue = c('Bladder','Bone/SoftTissue','Breast','Biliary','Cervix','Lymphoid','Myeloid','Colon/Rectum','Prostate','Esophagus','Stomach','CNS','Head/Neck','Kidney','Liver','Lung','Ovary','Pancreas','Skin','Thyroid','Uterus')  
TurnOverDays = c(200,5373,84.5,200,6,30,30,5,120,11,5.5,10000,16,1000,400,5143,11000,360,147,4138,4); length(TurnOverDays)
Turn = data.frame(CancerTissue,TurnOverDays)
Turn = Turn[order(Turn$TurnOverDays),]
HypMut = merge(HypMut,Turn,by.x = 'Tier2', by.y = 'CancerTissue')

HypMut$AhGhDummy = 0
for (i in 1:nrow(HypMut)) {if(HypMut$Subst[i] == 'T>C') {HypMut$AhGhDummy[i] = 1}}

#################################################
####### ANALYSES:
#################################################

##### 1:  if tissue-specific hypoxia is associated with turnover rate? YES, there is a trend

Agg = aggregate(HypMut$hypoxia_score_buffa, by = list(HypMut$Tier2,HypMut$TurnOverDays), FUN = median)
names(Agg) = c('CancerTissue','TurnOverDays','MedianHypoxiaScoreBuffa')
plot(Agg$TurnOverDays,Agg$MedianHypoxiaScoreBuffa)
cor.test(Agg$TurnOverDays,Agg$MedianHypoxiaScoreBuffa, method = 'spearman', alternative = 'less') # rho = -0.44, p = 0.02704

##### 2:  if slow- middle- and fast- dividing tissues (as in BioRxiv) have different hypoxia (N = 19, N = 828)?


##### 3A:  if there is correlation between VAF and hypoxia?  # YES!!! advanced cancers (with high VAF) are more hypoxic 
cor.test(HypMut$tumor_var_freq,HypMut$hypoxia_score_buffa, method = 'spearman') # very positive
Agg = aggregate(list(HypMut$hypoxia_score_buffa,HypMut$tumor_var_freq), by = list(HypMut$sample,HypMut$Tier2), FUN = median)
names(Agg) = c('sample','Tier2','MedianHypoxiaScoreBuffa','MedianVaf')
cor.test(Agg$MedianHypoxiaScoreBuffa,Agg$MedianVaf, method = 'spearman') # still positive

##### 3B:  if there is correlation between VAF and hypoxia within numerous cancer types? (the strongest correlation (minimal p) is in the most numerous kidney)
names(HypMut)
Patients = HypMut[,c(1,3,10,11)]; Patients = unique(Patients); nrow(Patients) # 828
Tissues = data.frame(table(Patients$Tier2)) # 19 vs 21?
Tissues = Tissues[order(-Tissues$Freq),]
NumerousTissues = Tissues[Tissues$Freq >= 50,]$Var1; length(NumerousTissues) 
NumerousTissues # Breast   Kidney   Liver    Lung     Ovary    Pancreas

cor.test(Agg[Agg$Tier2 == NumerousTissues[1],]$MedianHypoxiaScoreBuffa,Agg[Agg$Tier2 == NumerousTissues[1],]$MedianVaf, method = 'spearman') # 
cor.test(Agg[Agg$Tier2 == NumerousTissues[2],]$MedianHypoxiaScoreBuffa,Agg[Agg$Tier2 == NumerousTissues[2],]$MedianVaf, method = 'spearman') #  !!!
cor.test(Agg[Agg$Tier2 == NumerousTissues[3],]$MedianHypoxiaScoreBuffa,Agg[Agg$Tier2 == NumerousTissues[3],]$MedianVaf, method = 'spearman') #  !!!
cor.test(Agg[Agg$Tier2 == NumerousTissues[4],]$MedianHypoxiaScoreBuffa,Agg[Agg$Tier2 == NumerousTissues[4],]$MedianVaf, method = 'spearman') #
cor.test(Agg[Agg$Tier2 == NumerousTissues[5],]$MedianHypoxiaScoreBuffa,Agg[Agg$Tier2 == NumerousTissues[5],]$MedianVaf, method = 'spearman') #
cor.test(Agg[Agg$Tier2 == NumerousTissues[6],]$MedianHypoxiaScoreBuffa,Agg[Agg$Tier2 == NumerousTissues[6],]$MedianVaf, method = 'spearman') # 

##### 4:  if there is correlation between A>G and hypoxia (Whole dataset and within numerous cancer types)?
##### 5:  if there is correlation between Absolute number of mtDNA mutationsand hypoxia (Whole dataset and within numerous cancer types)?

names(HypMut)
Patients = HypMut[,c(1,3,10,11)]; Patients = unique(Patients); nrow(Patients) # 828
Tissues = data.frame(table(Patients$Tier2)) # 19 vs 21?
NumerousTissues = Tissues[Tissues$Freq >= 50,]$Var1; length(NumerousTissues) 
NumerousTissues # Breast   Kidney   Liver    Lung     Ovary    Pancreas

i = 1
temp = HypMut[HypMut$Tier2 == NumerousTissues[i],]
summary(glm(temp$AhGhDummy ~ temp$hypoxia_score_buffa, family = binomial()))










#################################################
####### OLD CODE:
#################################################


## AhGhfr is expected to be higher among high VAF and lower among hypoxic

TvVec = c('A>T','A>C','C>A','C>G','T>A','T>G','G>C','G>T')

table(HypMut$Subst) # light chain
str(HypMut$tumor_var_freq)
HypMut$AhGhDummy = 0
for (i in 1:nrow(HypMut)) {if(HypMut$Subst[i] == 'T>C') {HypMut$AhGhDummy[i] = 1}}
table(HypMut$AhGhDummy)
summary(glm(HypMut$AhGhDummy ~ HypMut$hypoxia_score_ragnum + HypMut$tumor_var_freq, family = binomial()))
summary(glm(HypMut$AhGhDummy ~ HypMut$tumor_var_freq, family = binomial())) # a bit

### frequencies of all four transitions positively correlate with hypoxic score (the higher the score => the higher VAF => early origin and/or more relaxed mtDNA selection in hypoxic cancers)
summary(lm(HypMut$hypoxia_score_ragnum ~ HypMut$tumor_var_freq)) # very positive => the higher the hypoxia the higher VAF (the older all mutations => originated at healthy tissues)
summary(lm(HypMut[HypMut$Subst == 'T>C',]$hypoxia_score_ragnum ~ HypMut[HypMut$Subst == 'T>C',]$tumor_var_freq)) # positive
summary(lm(HypMut[HypMut$Subst == 'C>T',]$hypoxia_score_ragnum ~ HypMut[HypMut$Subst == 'C>T',]$tumor_var_freq)) # positive
summary(lm(HypMut[HypMut$Subst == 'G>A',]$hypoxia_score_ragnum ~ HypMut[HypMut$Subst == 'G>A',]$tumor_var_freq)) # positive
summary(lm(HypMut[HypMut$Subst == 'A>G',]$hypoxia_score_ragnum ~ HypMut[HypMut$Subst == 'A>G',]$tumor_var_freq)) # positive
summary(lm(HypMut[HypMut$Subst %in% TvVec,]$hypoxia_score_ragnum ~ HypMut[HypMut$Subst  %in% TvVec,]$tumor_var_freq)) # positive

# emerging PolG signature? check it!

## can we associate hypoxia with cell division rate of each tissue - YES (correlation is weak, but expected direction)
# from Cancer.DifferencesBetweenCancerTypes.R


# analyses


# if we repeat the same with rare substitutions only - effect is even better - why?
summary(HypMut$tumor_var_freq) # 5.3700
Agg = aggregate(HypMut[HypMut$tumor_var_freq <1.79,]$hypoxia_score_buffa, by = list(HypMut[HypMut$tumor_var_freq <1.79,]$Tier2,HypMut[HypMut$tumor_var_freq <1.79,]$TurnOverDays), FUN = median)
names(Agg) = c('CancerTissue','TurnOverDays','MedianHypoxiaScoreBuffa')
plot(Agg$TurnOverDays,Agg$MedianHypoxiaScoreBuffa)
cor.test(Agg$TurnOverDays,Agg$MedianHypoxiaScoreBuffa, method = 'spearman', alternative = 'less') # 0.02704


## can we rerun the same lm where instead of T>C there is a hypoxia
summary(lm(HypMut$hypoxia_score_buffa ~ HypMut$tumor_var_freq+HypMut$TurnOverDays)) # positive both coefficients!!!
summary(lm(HypMut$hypoxia_score_buffa ~ 0 +  HypMut$tumor_var_freq+HypMut$TurnOverDays)) # positive both coefficients!!!

Final$OtherMut = 1-Final$AhGhfr
colors = c("red","black")

FinalNew = data.frame(Final$AhGhfr, Final$OtherMut, rownames(FinalNew), Final$hypoxia_score_buffa)    
FinalHypo = FinalNew[FinalNew$Final.hypoxia_score_buffa>=32,]
FinalHypo$Final.hypoxia_score_buffa = NULL
FinalNorm = FinalNew[FinalNew$Final.hypoxia_score_buffa<32,]
FinalNorm$Final.hypoxia_score_buffa = NULL

a = FinalHypo %>% gather(subtype, freq,  Final.AhGhfr:Final.OtherMut)

b = FinalNorm %>% gather(subtype, freq,  Final.AhGhfr:Final.OtherMut)

pdf("../../Body/4Figures/Cancer.Hypoxia.pdf", width=7, height=3)
ggplot(a, aes(fill=subtype, y=freq, x=rownames.FinalNew.)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_manual(values=c("red","black"))

ggplot(b, aes(fill=subtype, y=freq, x=rownames.FinalNew.)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_manual(values=c("red","black"))
dev.off()

