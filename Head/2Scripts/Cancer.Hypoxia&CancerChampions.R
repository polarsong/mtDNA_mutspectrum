rm(list=ls(all=TRUE))
library(dplyr)

Mut = read.table("../../Body/1Raw/mtDNA_snv_Oct2016.txt", head = TRUE, sep = '\t')  
length(unique(Mut$sample)) # 2177
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
  Final = rbind(Final,c(VecOfSamples[i],temp$tissue[1],temp$Tier2[1],AhGhfr,TotalMut,TotalMutOld,TotalMutYoung,temp$hypoxia_score_winter[1],temp$hypoxia_score_buffa[1],temp$hypoxia_score_ragnum[1]))
}  
names(Final)=c('sample','tissue','Tier2','AhGhfr','TotalMut','TotalMutOld','TotalMutYoung','hypoxia_score_winter','hypoxia_score_buffa','hypoxia_score_ragnum')  
Final$TotalMut = as.numeric(Final$TotalMut)
summary(Final$TotalMut)

#### FROM APMutSpectrum.Cancer.R:
# VecOfHighlyMutatedSamples.i. nrow.temp. FractionT_C
# 1 2779fa01-ac93-4e80-a997-3385f72172c3         33   0.8181818
# 2                            ICGC_MB34         24   0.3750000
# 3                            ICGC_0441         18   0.7222222

table(Final$Tier2)
FinalBreast = Final[Final$Tier2 == 'Breast',]
FinalBreast$hypoxia_score_buffa = as.numeric(FinalBreast$hypoxia_score_buffa)
cor.test(FinalBreast$TotalMut,FinalBreast$hypoxia_score_buffa, method = 'spearman') # weak positive
plot(FinalBreast$TotalMut,FinalBreast$hypoxia_score_buffa)
FinalBreast[FinalBreast$TotalMut > 30,]$hypoxia_score_buffa # = -4
nrow(FinalBreast[FinalBreast$hypoxia_score_buffa <= -4,])
nrow(FinalBreast) # 28/76 = -0.37 - почти половина, а не нормоксичный хвост как я ожидал...


Final$AhGhfr = as.numeric(Final$AhGhfr)
### out of 828 cancers there are 5 champions with >= 13 mutations. Do they have increased A>G? Are they less hypoxic (more normoxic)?
plot(Final$TotalMut,Final$hypoxia_score_buffa) # 5 champions with >= 13
nrow(Final[Final$TotalMut >= 13,])
summary(Final[Final$TotalMut >= 13,]$AhGhfr)
summary(Final[Final$TotalMut < 13,]$AhGhfr)



