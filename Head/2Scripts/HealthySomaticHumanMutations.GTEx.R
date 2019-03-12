rm(list=ls(all=TRUE))

############ read human ref seq:
RefSeq = read.table("../../Body/1Raw/chrM_refAllele.txt", header = FALSE, sep = '')
names(RefSeq)=c('Position','AncestralAllele')

############ read GTEx mutations:
SomMut = read.table("../../Body/1Raw/TissueSpecificMtDNAMutationsGTEx.txt", header = TRUE, sep = '\t')
SomMut$Position = gsub('_(.*)','',SomMut$Mutation)  # gsub("_(.*)",'',READS$patient)
SomMut$DerivedAllele = gsub('(.*)_','',SomMut$Mutation)  # gsub("_(.*)",'',READS$patient)

Som = merge(SomMut,RefSeq, by = 'Position')

############ 

Som = Som[Som$DerivedAllele != Som$AncestralAllele,]
Som$Substitution = paste(Som$AncestralAllele,Som$DerivedAllele, sep = '_')
table(Som$Substitution)
# A_C A_G A_T C_A C_G C_T G_A G_C G_T T_A T_C T_G 
# 43 473  49  49   6 336 728  63  57  42 693  26

Som$T_C = 0; Som$G_A = 0;
Som$C_T = 0; Som$A_G = 0;
Som$Tv = 0
for (i in 1:nrow(Som))
{
  if (Som$Substitution[i] == 'T_C') {Som$T_C[i] = 1;}
  if (Som$Substitution[i] == 'C_T') {Som$C_T[i] = 1;}
  if (Som$Substitution[i] == 'G_A') {Som$G_A[i] = 1;}
  if (Som$Substitution[i] == 'A_G') {Som$A_G[i] = 1;}
}

table(Som$tissue)
Som$Tissue = as.character(Som$tissue)
for (i in 1:nrow(Som))
{ # i = 1 
  if (length(grep('Brain',as.character(Som$tissue[i]))))       {Som$Tissue[i] = 'CNS';}
  if (length(grep('Colon',as.character(Som$tissue[i]))))       {Som$Tissue[i] = 'Colon/Rectum';}
  if (length(grep('Heart',as.character(Som$tissue[i]))))       {Som$Tissue[i] = 'Heart';}
  if (length(grep('Adipose',as.character(Som$tissue[i]))))     {Som$Tissue[i] = 'Adipose';}
  if (length(grep('Esophagus',as.character(Som$tissue[i]))))   {Som$Tissue[i] = 'Esophagus';}
  if (length(grep('Artery',as.character(Som$tissue[i]))))      {Som$Tissue[i] = 'Artery';}
  if (length(grep('Skin',as.character(Som$tissue[i]))))        {Som$Tissue[i] = 'Skin';}
  if (length(grep('Intestine',as.character(Som$tissue[i]))))   {Som$Tissue[i] = 'Intestine';}
  if (length(grep('Muscle - Skeletal',as.character(Som$tissue[i]))))   {Som$Tissue[i] = 'Muscle';}
  if (length(grep('Kidney - Cortex',as.character(Som$tissue[i]))))   {Som$Tissue[i] = 'Kidney';}
  if (length(grep('Cells',as.character(Som$tissue[i]))))       {Som$Tissue[i] = 'Cells';}
  if (length(grep('Breast',as.character(Som$tissue[i]))))       {Som$Tissue[i] = 'Breast';}
  if (length(grep('Gland',as.character(Som$tissue[i]))))       {Som$Tissue[i] = 'Gland';}
}
Data = data.frame(table(Som$Tissue))
Data = Data[Data$Freq > 0,]

### derive table is here: https://docs.google.com/document/d/1UECub1DIdmuXwPqDLK8WRcZ6uIjEQiVmHXN1_XNVvUU/edit?usp=sharing  
Tissue = c('Bladder','Bone/SoftTissue','Breast','Biliary','Cervix','Lymphoid','Myeloid','Colon/Rectum','Prostate','Esophagus','Stomach','CNS','Head/Neck','Kidney','Liver','Lung','Ovary','Pancreas','Skin','Thyroid','Uterus')  
TurnOverDays = c(200,5373,84.5,200,6,30,30,5,120,11,5.5,10000,16,1000,400,5143,11000,360,147,4138,4)
Turn = data.frame(Tissue,TurnOverDays)
Turn$Type = 'NA'
for (i in 1:nrow(Turn))
{
  if (Turn$TurnOverDays[i] <= 30) {Turn$Type[i] = 'FAST';}
  if (Turn$TurnOverDays[i] > 30 & Turn$TurnOverDays[i] <= 1000) {Turn$Type[i] = 'MIDDLE';}
  if (Turn$TurnOverDays[i] > 1000) {Turn$Type[i] = 'SLOW';}
}

Som = merge(Som,Turn)
Agg = aggregate(list(Som$T_C,Som$G_A), by = list(Som$Type), FUN = mean); names(Agg)=c('TissueType','T_C','G_A')
Agg$TC_GA = Agg$T_C/Agg$G_A

cor.test(Agg$TissueTurnOverDays,Agg$T_C, method = 'spearman')

Agg = aggregate(list(Som$T_C,Som$G_A), by = list(Som$Tissue), FUN = mean); names(Agg)=c('Tissue','T_C','G_A')

Agg1 = merge(Agg,Turn, by = 'Tissue')
Agg1 = Agg1[order(Agg1$TurnOverDays),]
Agg1$TC_GA = Agg1$T_C/Agg1$G_A
cor.test(Agg1$TurnOverDays,Agg1$T_C, method = 'spearman')


FAST = ALL[ALL$TurnOverDays <= 30,]                      
MIDDLE = ALL[ALL$TurnOverDays > 30 & ALL$TurnOverDays <= 1000,]  # 4138/360
SLOW = ALL[ALL$TurnOverDays > 1000,]


