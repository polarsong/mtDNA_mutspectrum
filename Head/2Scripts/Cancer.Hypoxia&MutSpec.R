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
VecOfSamples = unique(HypMut$sample); length(VecOfSamples) # 828 samples

Final = data.frame()
for (i in 1:length(VecOfSamples))
{ # i = 1
  temp = HypMut[HypMut$sample == VecOfSamples[i],]
  AhGhfr = nrow(temp[temp$Subst == 'T>C',])/nrow(temp)
  TotalMut = nrow(temp)
  TotalMutOld = nrow(temp[temp$tumor_var_freq > 5.7,])
  TotalMutYoung = nrow(temp[temp$tumor_var_freq <= 5.7,])
  Final = rbind(Final,c(VecOfSamples[i],AhGhfr,TotalMut,TotalMutOld,TotalMutYoung,temp$hypoxia_score_winter[1],temp$hypoxia_score_buffa[1],temp$hypoxia_score_ragnum[1]))
}  
names(Final)=c('SampleId','AhGhfr','TotalMut','TotalMutOld','TotalMutYoung','hypoxia_score_winter','hypoxia_score_buffa','hypoxia_score_ragnum')  

##### all hypoxia scores well correlate with each other
cor.test(Final$hypoxia_score_buffa,Final$hypoxia_score_winter,method = 'spearman')
cor.test(Final$hypoxia_score_buffa,Final$hypoxia_score_ragnum,method = 'spearman')
cor.test(Final$hypoxia_score_winter,Final$hypoxia_score_ragnum,method = 'spearman')

##### hypoxia scores negatively correlate with TotalMutYoung (new mutations do not come)) and positively with TotalMutOld (old mutations have time to go to expansion) !!!!!?????
# or, probably, because fastly-dividing healthy tissues (those, that become hypoxic with a time) divide fast and get more old mutations. 
cor.test(Final$TotalMut,Final$hypoxia_score_ragnum,method = 'spearman')
cor.test(Final$TotalMut,Final$hypoxia_score_winter,method = 'spearman')
cor.test(Final$TotalMut,Final$hypoxia_score_buffa,method = 'spearman')  # PAPER
cor.test(Final$TotalMutOld,Final$hypoxia_score_buffa,method = 'spearman') # positive !!! PAPER
cor.test(Final$TotalMutOld,Final$hypoxia_score_winter,method = 'spearman') # positive !!!
cor.test(Final$TotalMutOld,Final$hypoxia_score_ragnum,method = 'spearman') # positive !!!
nrow(Final)  # PAPER
cor.test(Final$TotalMutYoung,Final$hypoxia_score_buffa,method = 'spearman') # negative !!! PAPER
cor.test(Final$TotalMutYoung,Final$hypoxia_score_winter,method = 'spearman') # negative !!!
cor.test(Final$TotalMutYoung,Final$hypoxia_score_ragnum,method = 'spearman') # negative !!!

cor.test(Final$TotalMutYoung,Final$TotalMutOld,method = 'spearman') # nothing


## hypoxia shows wweak trend to decreased A>G
cor.test(Final$AhGhfr,Final$hypoxia_score_ragnum,method = 'spearman')
cor.test(Final$AhGhfr,Final$hypoxia_score_winter,method = 'spearman')
cor.test(Final$AhGhfr,Final$hypoxia_score_buffa,method = 'spearman')

## compare AhGhfr among normoxic (hypoxia scores < 75%) and hypoxic (hypoxic scores >= 75%) cancers:  !!!!!
# we follow the buffa as the main score (as authors did in NatCom)  !!!!!! THIS IS GOOD !!!!!

wilcox.test(Final[Final$hypoxia_score_buffa < quantile(Final$hypoxia_score_buffa,0.9),]$AhGhfr,Final[Final$hypoxia_score_buffa >= quantile(Final$hypoxia_score_buffa,0.9),]$AhGhfr)
boxplot(Final[Final$hypoxia_score_buffa < quantile(Final$hypoxia_score_buffa,0.9),]$AhGhfr,Final[Final$hypoxia_score_buffa >= quantile(Final$hypoxia_score_buffa,0.9),]$AhGhfr, notch = TRUE, col =c('blue','light blue'))   
# ADD VIOLIN HERE (same code as one line above)! https://github.com/polarsong/mtDNA_mutspectrum/blob/Cancer/Head/2Scripts/Cancer.Hypoxia%26MutSpec.R

summary(Final[Final$hypoxia_score_buffa < quantile(Final$hypoxia_score_buffa,0.9),]$AhGhfr) # 0.28
summary(Final[Final$hypoxia_score_buffa >= quantile(Final$hypoxia_score_buffa,0.9),]$AhGhfr) # 0.20
nrow(Final[Final$hypoxia_score_buffa < quantile(Final$hypoxia_score_buffa,0.9),])
nrow(Final[Final$hypoxia_score_buffa >= quantile(Final$hypoxia_score_buffa,0.9),])

boxplot(Final[Final$hypoxia_score_buffa < quantile(Final$hypoxia_score_buffa,0.75),]$AhGhfr,Final[Final$hypoxia_score_buffa >= quantile(Final$hypoxia_score_buffa,0.75),]$AhGhfr, notch = TRUE)
wilcox.test(Final[Final$hypoxia_score_buffa < quantile(Final$hypoxia_score_buffa,0.75),]$AhGhfr,Final[Final$hypoxia_score_buffa >= quantile(Final$hypoxia_score_buffa,0.75),]$AhGhfr, alternative = 'greater')

## what about all mutations in general? weak decrease
boxplot(Final[Final$hypoxia_score_buffa < quantile(Final$hypoxia_score_buffa,0.9),]$TotalMut,Final[Final$hypoxia_score_buffa >= quantile(Final$hypoxia_score_buffa,0.9),]$TotalMut, notch = TRUE)
wilcox.test(Final[Final$hypoxia_score_buffa < quantile(Final$hypoxia_score_buffa,0.9),]$TotalMut,Final[Final$hypoxia_score_buffa >= quantile(Final$hypoxia_score_buffa,0.9),]$TotalMut)
boxplot(Final[Final$hypoxia_score_buffa < quantile(Final$hypoxia_score_buffa,0.75),]$TotalMut,Final[Final$hypoxia_score_buffa >= quantile(Final$hypoxia_score_buffa,0.75),]$TotalMut, notch = TRUE)
wilcox.test(Final[Final$hypoxia_score_buffa < quantile(Final$hypoxia_score_buffa,0.75),]$TotalMut,Final[Final$hypoxia_score_buffa >= quantile(Final$hypoxia_score_buffa,0.75),]$TotalMut, alternative = 'greater')
cor.test(Final$hypoxia_score_buffa,Final$TotalMut, method = 'spearman')

## try to add effect of heteroplasmy:  AhGhfr is expected to be higher among high VAF and lower among hypoxic

TvVec = c('A>T','A>C','C>A','C>G','T>A','T>G','G>C','G>T')

str(HypMut$tumor_var_freq)
HypMut$AhGhDummy = 0
for (i in 1:nrow(HypMut)) {if(HypMut$Subst[i] == 'A>G') {HypMut$AhGhDummy[i] = 1}}
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

CancerTissue = c('Bladder','Bone/SoftTissue','Breast','Biliary','Cervix','Lymphoid','Myeloid','Colon/Rectum','Prostate','Esophagus','Stomach','CNS','Head/Neck','Kidney','Liver','Lung','Ovary','Pancreas','Skin','Thyroid','Uterus')  
TurnOverDays = c(200,5373,84.5,200,6,30,30,5,120,11,5.5,10000,16,1000,400,5143,11000,360,147,4138,4); length(TurnOverDays)
Turn = data.frame(CancerTissue,TurnOverDays)
Turn = Turn[order(Turn$TurnOverDays),]

HypMut = merge(HypMut,Turn,by.x = 'Tier2', by.y = 'CancerTissue')
Agg = aggregate(HypMut$hypoxia_score_buffa, by = list(HypMut$Tier2,HypMut$TurnOverDays), FUN = median)
names(Agg) = c('CancerTissue','TurnOverDays','MedianHypoxiaScoreBuffa')
plot(Agg$TurnOverDays,Agg$MedianHypoxiaScoreBuffa)
cor.test(Agg$TurnOverDays,Agg$MedianHypoxiaScoreBuffa, method = 'spearman', alternative = 'less') # rho = -0.44, p = 0.02704

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

