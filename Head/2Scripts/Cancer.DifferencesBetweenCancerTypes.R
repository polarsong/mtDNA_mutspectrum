###################################
###### 26.03.2018: Cancer mut Spectrum - between different cancers 
###################################

#### nucleotide content of the human genome: A T G C: 4993	3871	2159	5357
  rm(list=ls(all=TRUE))
  
  ALL = read.table("../../Body/1Raw/mtDNA_snv_Oct2016.txt", head = TRUE, sep = '\t')  # 7611
  
  pdf("../../Body/4Figures/Cancer.DifferencesBetweenCancerTypes.01.pdf" , height = 15, width = 30)
  
  ### DERIVE NECESSARY TRAITS:
  ALL$TumorVarFreq = ALL$tumor_reads2/(ALL$tumor_reads1 + ALL$tumor_reads2); summary(ALL$TumorVarFreq)  # 0.01000 0.01738 0.04540 0.20268 0.26278 0.99864
  ALL$NormalVarFreq = ALL$normal_reads2/(ALL$normal_reads1 + ALL$normal_reads2); summary(ALL$NormalVarFreq) 
  ALL$Subs = paste(ALL$ref,ALL$var, sep = '_'); table(ALL$Subs)
  VecOfTransitionSubstitutions = c('A_G','G_A','C_T','T_C') # all tr
  VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # ALL Tv
  CancerType = unique(ALL$tissue); length(CancerType)    # 40
  CancerTissue = unique(ALL$Tier2); length(CancerTissue) # 21
  ALL$CancerTissue = ALL$Tier2; unique(ALL$CancerTissue) # 
  ALL$CancerType = ALL$tissue;
  ALL$CancerTypeAndTissues = paste(ALL$CancerType,ALL$CancerTissue, sep = '_')
  table(ALL$CancerTypeAndTissues)
  #BLCA_Bladder BOCA_Bone/SoftTissue          BRCA_Breast         BTCA_Biliary          CESC_Cervix        CLLE_Lymphoid         CMDI_Myeloid    COAD_Colon/Rectum 
  #99                   90                  689                   28                   38                  118                   59                  167 
  #DLBC_Lymphoid        EOPC_Prostate       ESAD_Esophagus         GACA_Stomach              GBM_CNS       HNSC_Head/Neck          KICH_Kidney          KIRC_Kidney 
  #16                  122                  409                   80                   46                  111                  233                  112 
  #KIRP_Kidney         LAML_Myeloid              LGG_CNS           LICA_Liver           LIHC_Liver           LINC_Liver         LIRI_Biliary           LIRI_Liver 
  #191                   19                   30                   18                  238                  102                   62                  915 
  #LUAD_Lung            LUSC_Lung        MALY_Lymphoid            MELA_Skin       ORCA_Head/Neck             OV_Ovary        PACA_Pancreas        PAEN_Pancreas 
  #141                  127                  182                  144                   35                  402                  707                  154 
  #PBCA_CNS        PRAD_Prostate    READ_Colon/Rectum          RECA_Kidney SARC_Bone/SoftTissue            SKCM_Skin         STAD_Stomach         THCA_Thyroid 
  #186                  617                   57                  286                   65                   73                  130                  168 
  #UCEC_Uterus 
  #145 
  
########### LIFETIME RISK OF DEVELOPING CANCER BASED ON TISSUE:
  
  unique(ALL$CancerTissue)
  # Bladder BLCA (Bladder Urothelial Cancer) 0.023 - lifetime risk (https://seer.cancer.gov/statfacts/html/urinb.html)
  # Bone/SoftTissue BOCA (Bone Cancer) 0.01 (https://seer.cancer.gov/statfacts/)  The term “sarcoma” encompasses a broad family of rare cancers that can affect soft tissue or bone throughout the body, and sometimes both. 
  # Breast BRCA (Breast Cancer) 0.124 of women (https://seer.cancer.gov/statfacts/html/breast.html)
  # Biliary (~ Liver and Intrahepatic Bile Duct Cancer) 0.01 (https://seer.cancer.gov/statfacts/html/livibd.html)
  # Cervix (Cervical Cancer) 0.006 (https://seer.cancer.gov/statfacts/html/cervix.html); This cancer forms in tissues of the cervix (the organ connecting the uterus and vagina). It is usually a slow-growing cancer that may not have symptoms but can be found with regular Pap tests (a procedure in which cells are scraped from the cervix and looked at under a microscope). Cervical cancer is almost always caused by human papillomavirus (HPV) infection.
  # Lymphoid 0.002 + 0.021 https://seer.cancer.gov/statfacts/html/hodg.html + https://seer.cancer.gov/statfacts/html/nhl.html 
  # Myeloid 0.008 (https://seer.cancer.gov/statfacts/html/mulmy.html)
  # colon/rectum 0.042 (https://seer.cancer.gov/statfacts/html/colorect.html)
  # Prostate 0.112 of men (https://seer.cancer.gov/statfacts/html/prost.html)
  # Esophagus 0.005 (https://seer.cancer.gov/statfacts/html/esoph.html)
  # Stomach 0.008 (https://seer.cancer.gov/statfacts/html/esoph.html)
  # CNS 0.006 (https://seer.cancer.gov/statfacts/html/brain.html)
  # Head/Neck (~ Soft Tissue including Heart Cancer) 0.003 (https://seer.cancer.gov/statfacts/html/soft.html)
  # Kidney (~ Kidney and Renal Pelvis Cancer) 0.017 (https://seer.cancer.gov/statfacts/html/kidrp.html)
  # Liver (~ Liver and Intrahepatic Bile Duct Cancer)  0.01 (https://seer.cancer.gov/statfacts/html/livibd.html)
  # Lung (Lung and Bronchus Cancer) 0.062 (https://seer.cancer.gov/statfacts/html/lungb.html)
  # Ovary 0.013 of women (https://seer.cancer.gov/statfacts/html/ovary.html)
  # Pancreas 0.016 (https://seer.cancer.gov/statfacts/html/pancreas.html)
  # Skin (~Melanoma of the Skin) 0.023 (https://seer.cancer.gov/statfacts/html/melan.html)
  # Thyroid 0.012 (https://seer.cancer.gov/statfacts/html/thyro.html)
  # Uterus 0.029 of women (https://seer.cancer.gov/statfacts/html/corp.html)
  
  CancerTissue = c('Bladder','Bone/SoftTissue','Breast','Biliary','Cervix','Lymphoid','Myeloid','Colon/Rectum','Prostate','Esophagus','Stomach','CNS','Head/Neck','Kidney','Liver','Lung','Ovary','Pancreas','Skin','Thyroid','Uterus')  
  LifeTimeRisk = c(0.023,0.01,0.124,0.01,0.006,0.023,0.008,0.042,0.112,0.005,0.008,0.006,0.003,0.017,0.01,0.062,0.013,0.016,0.023,0.012,0.029)
  LifeRisk = data.frame(CancerTissue,LifeTimeRisk)
  ALL = merge(ALL, LifeRisk, by = 'CancerTissue')
  
#################### NumOfCellDivPerLife (not for all!!!!)
######### division rate: Cancer number of divisions (from Table S1 one before last column Tomasetti and Vogelstein 2015 science)
  
Cells = data.frame(matrix( c("Thyroid","Head/Neck","Ovary","Esophagus","Pancreas","Skin","Lung","Liver","Colon/Rectum","Lymphoid","CNS","Bone/SoftTissue","Myeloid",                             
                               7,1720,0,1390,80,199,5.6,88,5840,960,0,5,960), ncol = 2))
names(Cells) = c('CancerTissue','NumOfCellDivPerLife')
Cells$NumOfCellDivPerLife = as.numeric(as.character(Cells$NumOfCellDivPerLife))
ALL = merge(ALL,Cells, all.x = TRUE)
  
####################### TURNOVER OF CELLS (number of days)
### derive table is here: https://docs.google.com/document/d/1UECub1DIdmuXwPqDLK8WRcZ6uIjEQiVmHXN1_XNVvUU/edit?usp=sharing  
CancerTissue = c('Bladder','Bone/SoftTissue','Breast','Biliary','Cervix','Lymphoid','Myeloid','Colon/Rectum','Prostate','Esophagus','Stomach','CNS','Head/Neck','Kidney','Liver','Lung','Ovary','Pancreas','Skin','Thyroid','Uterus')  
TurnOverDays = c(200,5373,84.5,200,6,30,30,5,120,11,5.5,10000,16,1000,400,5143,11000,360,147,4138,4)
Turn = data.frame(CancerTissue,TurnOverDays)
Turn = Turn[order(Turn$TurnOverDays),]
ALL = merge(ALL,Turn)

######## glycolisis: from Andrey Yurchenko and https://www.nature.com/articles/s41467-018-07232-8#MOESM5
  
CancerType =c('PRAD','LUNG','COAD','BRCA','KIRC','BLAD','CESC','CHOL','COADREAD','ESCA','GBM','HNSC','THCA','THYM','STAD','SKCM','SARC','READ','PCPG','PAAD','LUSC','LUAD','LIHC','KIRP','KICH','UCEC','MSKCCTvN')
Glycolysis=c(12.38978,21.35984,59.51087,129.9068,99.2046,12.13424,14.89038,49.08608,26.55547,7.960641,3.820559,60.00785,46.08075,20.28906,19.16362,3.914527,10.47275,25.35386,43.00674,0.5165879,64.19461,25.59793,14.61152,77.26818,74.4774,29.7614,3.739063)
OxidativePhosphorylation =c(15.76979,9.893736,14.81407,14.28531,249.1245,2.216577,3.199398,12.84608,7.390275,8.302373,6.407421,13.32523,4.457943,0.9468787,13.56102,2.319901,1.623032,11.08796,40.05531,0.3871189,26.50188,12.20636,13.09163,166.9113,51.49203,14.65809,2.185885)

Glyc = data.frame(CancerType,Glycolysis,OxidativePhosphorylation)
ALL = merge(ALL,Glyc, all.x=TRUE)
  
####### mtDNA copies: https://www.ncbi.nlm.nih.gov/pubmed/26901439 (figure 3 gives a rank of cancer types)
CancerType = c('KIRC','BRCA','BLCA','LIHC','HNSC','ESCA','KIRP','STAD','UCEC','KICH','COAD','THCA','PAAD','PRAD','LUAD')
mtDNACopyNumberRank = c(1:15)
CopyNumber = data.frame(CancerType,mtDNACopyNumberRank)
ALL = merge(ALL,CopyNumber, all.x=TRUE)
  
####### can we derive mtDNA copies from the table?
####### for example compare total coverage in tumor versus total coverage in normal versus cancer types?
  
ALL$CovNor = ALL$normal_reads1 + ALL$normal_reads2
ALL$CovTum = ALL$tumor_reads1 + ALL$tumor_reads2
ALL$CovTumToNor = log2(ALL$CovTum/ALL$CovNor)
boxplot(CovTumToNor ~ CancerType, data = ALL, outline = FALSE, notch = TRUE, ylab = 'log2(TotalCovTum/TotalCovNor)')
abline(h=0, col = 'red')
  
# get medain of log2(TotalCovTum/TotalCovNor) for each cancer
AggMed = aggregate(ALL$CovTumToNor, by = list(ALL$CancerType,ALL$sample), FUN = median)
AggMed = aggregate(AggMed$x, by = list(AggMed$Group.1), FUN = median)
names(AggMed)=c('CancerType','MedianCovTumToNor')
AggMed = AggMed[order(AggMed$MedianCovTumToNor),]
  
ALL = merge(ALL,AggMed, all.x=TRUE)
  
#######################
####### ANALYSES
#######################

ALL_T_C = ALL[ALL$Subs == 'T_C',]; ALL_T_C$T_C = 1
ALL_not_T_C = ALL[!ALL$Subs == 'T_C',]; ALL_not_T_C$T_C = 0
ALL = rbind(ALL_T_C,ALL_not_T_C)

ALL_G_A = ALL[ALL$Subs == 'G_A',]; ALL_G_A$G_A = 1
ALL_not_G_A = ALL[!ALL$Subs == 'G_A',]; ALL_not_G_A$G_A = 0
ALL = rbind(ALL_G_A,ALL_not_G_A)

ALL_Ts = ALL[ALL$Subs %in% VecOfTransitionSubstitutions,]; ALL_Ts$Ts = 1  
ALL_Tv = ALL[ALL$Subs %in% VecOfTransversionSubstitutions,]; ALL_Tv$Ts = 0  
ALL = rbind(ALL_Ts,ALL_Tv)
  
### TEST -1: TurnOver and Ts/Tv
length(unique(ALL$TurnOverDays)) # 19 different rates (21 tissues) = 7 in each group
sort(unique(ALL$TurnOverDays))
Turn
FAST = ALL[ALL$TurnOverDays <= 30,]                      
MIDDLE = ALL[ALL$TurnOverDays > 30 & ALL$TurnOverDays <= 1000,]  
SLOW = ALL[ALL$TurnOverDays > 1000,]

TEMP = FAST
TsTv = nrow(TEMP[TEMP$Subs %in% VecOfTransitionSubstitutions,]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); TsTv # 10.02817
TCTv = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); TCTv # 2.816901
GATv = nrow(TEMP[TEMP$Subs == 'G_A',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); GATv # 5.014085
TCGA = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  == 'G_A',]); TCGA # 0.5617978
TEMP = MIDDLE
TsTv = nrow(TEMP[TEMP$Subs %in% VecOfTransitionSubstitutions,]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); TsTv # 12.56941
TCTv = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); TCTv # 3.827195
GATv = nrow(TEMP[TEMP$Subs == 'G_A',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); GATv # 6.634561
TCGA = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  == 'G_A',]); TCGA # 0.5768574
TEMP = SLOW
TsTv = nrow(TEMP[TEMP$Subs %in% VecOfTransitionSubstitutions,]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); TsTv # 14.49383
TCTv = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); TCTv # 4.555556
GATv = nrow(TEMP[TEMP$Subs == 'G_A',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); GATv # 7.518519
TCGA = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  == 'G_A',]); TCGA # 0.6059113

####### bootstrep to perturb these three numbers (Ts/Tv)? to draw three boxplots!!!!
VecFast = c()
for (i in 1:1000) {
  TEMP = FAST[sample(nrow(FAST),nrow(FAST)/2),]; 
  TsTv = nrow(TEMP[TEMP$Subs %in% VecOfTransitionSubstitutions,]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]) 
  VecFast = c(VecFast,TsTv) }
summary(VecFast)

VecMiddle = c()
for (i in 1:1000) {
  TEMP = MIDDLE[sample(nrow(MIDDLE),nrow(MIDDLE)/2),]; 
  TsTv = nrow(TEMP[TEMP$Subs %in% VecOfTransitionSubstitutions,]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]) 
  VecMiddle = c(VecMiddle,TsTv) }
summary(VecMiddle)

VecSlow = c()
for (i in 1:1000) {
  TEMP = SLOW[sample(nrow(SLOW),nrow(SLOW)/2),]; 
  TsTv = nrow(TEMP[TEMP$Subs %in% VecOfTransitionSubstitutions,]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]) 
  VecSlow = c(VecSlow,TsTv) }
summary(VecSlow)

boxplot(VecFast,VecMiddle,VecSlow, notch = TRUE) # dev.off()

####### logistic regression with all 21 Turnovers and 2 Dummy 

FAST$TurnOverDummyFast = 1; FAST$TurnOverDummySlow = 0; FAST$TurnOverRank = 1
MIDDLE$TurnOverDummyFast = 0; MIDDLE$TurnOverDummySlow = 0; MIDDLE$TurnOverRank = 2
SLOW$TurnOverDummyFast = 0; SLOW$TurnOverDummySlow = 1; SLOW$TurnOverRank = 3
ALL = rbind(FAST, MIDDLE, SLOW)

glm_1a <-glm(ALL$T_C ~ scale(ALL$TurnOverDays) + scale(ALL$TumorVarFreq), family = binomial())  # total number of mutations? total disruption?
summary(glm_1a) # both are significant

glm_1b <-glm(ALL$T_C ~ scale(ALL$TurnOverRank) + scale(ALL$TumorVarFreq), family = binomial())  # total number of mutations? total disruption?
summary(glm_1b) # both are significant

glm_2 <-glm(ALL$G_A ~ ALL$TurnOverDays + ALL$TumorVarFreq, family = binomial())  # total number of mutations? total disruption?
summary(glm_2) # TurnOverDays are not siginficant, VAF is negatively associated

glm_3 <-glm(ALL$Ts ~ ALL$TurnOverDays + ALL$TumorVarFreq, family = binomial())  # total number of mutations? total disruption?
summary(glm_3) # only VAF

glm_4 <-glm(ALL$Ts ~ ALL$TurnOverDummyFast + ALL$TumorVarFreq, family = binomial())  # total number of mutations? total disruption?
summary(glm_4) # both!

glm_4 <-glm(ALL$Ts ~ ALL$TurnOverDummySlow + ALL$TumorVarFreq, family = binomial())  # total number of mutations? total disruption?
summary(glm_4) # only VAF

glm_5a <-glm(ALL$T_C ~ ALL$TurnOverDummyFast + ALL$TumorVarFreq, family = binomial())  # total number of mutations? total disruption?
summary(glm_5a) # both!

glm_5b <-glm(ALL$T_C ~ ALL$TurnOverDummySlow + ALL$TumorVarFreq, family = binomial())  # total number of mutations? total disruption?
summary(glm_5b) # only VAF

glm_6 <-glm(ALL$G_A ~ ALL$TurnOverDummyFast + ALL$TumorVarFreq, family = binomial())  # total number of mutations? total disruption?
summary(glm_6) # both are negative!

glm_6 <-glm(ALL$G_A ~ ALL$TurnOverDummySlow + ALL$TumorVarFreq, family = binomial())  # total number of mutations? total disruption?
summary(glm_6) # only VAF is negative

#### plot glm_1

ALL$prob <- predict(glm_1a, type = 'response') # sqrt(P-EN/pi)
summary(ALL$prob) #  0.2616  0.2632  0.2667  0.2785  0.2918  0.3622
ALL$area = (ALL$prob - min(ALL$prob))/(max(ALL$prob)-min(ALL$prob)) # normilaze data to 0-1 range
summary(ALL$area)
plot(log2(ALL$TurnOverDays) ~ log2(ALL$TumorVarFreq), pch = '',  xlab = 'VAF', ylim = c(1,15), xlim = c(-6.7,0.2), ylab = 'TurnOverDays', main = 'T>C') #  dev.off()
symbols(log2(ALL[ALL$T_C == 0,]$TumorVarFreq), log2(ALL[ALL$T_C == 0,]$TurnOverDays), circles = ALL[ALL$T_C == 0,]$area, add = TRUE, fg = rgb(0.1,0.1,0.1,0.1))
symbols(log2(ALL[ALL$T_C == 1,]$TumorVarFreq), log2(ALL[ALL$T_C == 1,]$TurnOverDays), circles = ALL[ALL$T_C == 1,]$area, add = TRUE, fg = rgb(1,0.0,0.0,0.2))

ALL$prob <- predict(glm_1b, type = 'response') # sqrt(P-EN/pi)
summary(ALL$prob) #  0.2616  0.2632  0.2667  0.2785  0.2918  0.3622
ALL$area = (ALL$prob - min(ALL$prob))/(max(ALL$prob)-min(ALL$prob)) # normilaze data to 0-1 range
summary(ALL$area)
plot(ALL$TurnOverRank ~ log2(ALL$TumorVarFreq), pch = '', xlim = c(-6.7,0.2), xlab = 'VAF', ylab = 'TurnOverDays', main = 'prob (T>C)') #  dev.off()
symbols(log2(ALL[ALL$T_C == 0,]$TumorVarFreq), ALL[ALL$T_C == 0,]$TurnOverRank, circles = ALL[ALL$T_C == 0,]$area, add = TRUE, fg = rgb(0.1,0.1,0.1,0.1))
symbols(log2(ALL[ALL$T_C == 1,]$TumorVarFreq), ALL[ALL$T_C == 1,]$TurnOverRank, circles = ALL[ALL$T_C == 1,]$area, add = TRUE, fg = rgb(1,0.0,0.0,0.2))

cor.test(ALL[ALL$TurnOverRank == 1,]$TumorVarFreq,ALL[ALL$TurnOverRank == 1,]$T_C, method='spearman') # pos
cor.test(ALL[ALL$TurnOverRank == 2,]$TumorVarFreq,ALL[ALL$TurnOverRank == 2,]$T_C, method='spearman') # pos
cor.test(ALL[ALL$TurnOverRank == 3,]$TumorVarFreq,ALL[ALL$TurnOverRank == 3,]$T_C, method='spearman') # nonsign

### TurnOver may correlate with predisposition to switch metabolism, for example number of copies seems is increasing in slow-dividing tissues (they are more aerobic - they need mtDNA).
cor.test(ALL$CovTumToNor,ALL$TurnOverDays,method='spearman') # very significant positive: the slower turnover the more mtDNA in tumor versus normal
boxplot(ALL$CovTumToNor ~ ALL$TurnOverDays)

dev.off()


########################################################
############## OTHER MORE NOISY IDEAS::: glycolysis, mtDNA copies etx
########################################################

####### rank cor
Final = c()
CancerVec = unique(ALL$CancerTissue)
CancerTypeVec = unique(ALL$CancerType)  

for (i in 1:length(CancerVec))
{ # i = 1
  TEMP = ALL[ALL$CancerTissue == CancerVec[i],]
  NRow = nrow(TEMP)
  TurnOver = TEMP$TurnOverDays[1]
  CancerTissue = as.character(TEMP$CancerTissue[1])
  TsTv = nrow(TEMP[TEMP$Subs %in% VecOfTransitionSubstitutions,]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,])  #
  # TsTv = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,])  # 
  OneLine = c(as.character(CancerVec[i]),CancerTissue,TsTv,TurnOver,NRow)
  Final = rbind(Final,OneLine)
}
Final = data.frame(Final)
names(Final)=c('CancerType','CancerTissue','TsTv','TurnOverDays','NRow')
Final$TsTv = as.numeric(as.character(Final$TsTv))
Final$TurnOverDays = as.numeric(as.character(Final$TurnOverDays))
Final$NRow = as.numeric(as.character(Final$NRow))
Final = Final[order(Final$TurnOverDays),]

cor.test(Final$TsTv,Final$TurnOverDays,method = 'spearman')
plot(log2(Final$TsTv),log2(Final$TurnOverDays), pch = '') # dev.off()
text(log2(Final$TsTv),log2(Final$TurnOverDays),Final$CancerTissue) # dev.off()

boxplot(Final[Final$TurnOverDays < 2000,]$TsTv,Final[Final$TurnOverDays > 2000,]$TsTv)
wilcox.test(Final[Final$TurnOverDays < 2000,]$TsTv,Final[Final$TurnOverDays > 2000,]$TsTv)

glm_1 <-glm(ALL$T_C ~ ALL$TurnOverDays + ALL$TumorVarFreq, family = binomial())  # total number of mutations? total disruption?
summary(glm_1)

### TEST 0: lifetyme risk and T_C/G_A
Final = c()  
CancerVec = unique(ALL$CancerTissue)  
for (i in 1:length(CancerVec))
{ # i = 1
  TEMP = ALL[ALL$CancerTissue == CancerVec[i],]
  TC_GA = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs == 'G_A',]) 
  OneLine = c(as.character(CancerVec[i]),TC_GA)
  Final = rbind(Final,OneLine)
}
Final = data.frame(Final)
names(Final)=c('CancerTissue','TC_GA')
Final$TC_GA = as.numeric(as.character(Final$TC_GA))

Risk =  aggregate(ALL$T_C, by = list(ALL$sample,ALL$CancerTissue,ALL$LifeTimeRisk), FUN = mean)
Risk =  aggregate(Risk$x, by = list(Risk$Group.2,Risk$Group.3), FUN = median) 
names(Risk)=c('CancerTissue','LifeTimeRisk','FrT_C')
cor.test(Risk$LifeTimeRisk,Risk$FrT_C, method = 'spearman') # nonsign positive Try T_C/G_A
Risk = merge(Risk, Final, by = 'CancerTissue')
cor.test(Risk$LifeTimeRisk,Risk$TC_GA, method = 'spearman') # nonsign positive Try T_C/G_A
cor.test(Risk$FrT_C,Risk$TC_GA, method = 'spearman') # positive  - the higher T_C, the highet TC_GA.


# more fastly dividing are more aerobic?
### TEST 1: are there correlations between division rate, glycolisis, mtDNA copies etc between cancer types? 
#  I would think that fastly dividing are more glycolitic?
# nonsignificant positive - probably because there are many NAs

DivisionGlycolysis = ALL[grep("CancerType|CancerTissue|CancerTypeAndTissue|NumOfCellDivPerLife|Glycolysis|mtDNACopyNumberRank|MedianCovTumToNor|OxidativePhosphorylation|LifeTimeRisk",colnames(ALL))]
DivisionGlycolysis = unique(DivisionGlycolysis)
cor.test(DivisionGlycolysis$NumOfCellDivPerLife,DivisionGlycolysis$LifeTimeRisk, method='spearman') # nothing, 
cor.test(DivisionGlycolysis$NumOfCellDivPerLife,DivisionGlycolysis$Glycolysis, method='spearman') # positive = the more divisions, the more we can destroy glycolitic pathway
cor.test(DivisionGlycolysis$NumOfCellDivPerLife,DivisionGlycolysis$mtDNACopyNumberRank, method='spearman') # negative - the more divisions, the less mtDNA
cor.test(DivisionGlycolysis$Glycolysis,DivisionGlycolysis$mtDNACopyNumberRank, method='spearman') # negative - the stronger destryoed the glycilitic pathway, the less mtDNA. The best correlation!
cor.test(DivisionGlycolysis$MedianCovTumToNor,DivisionGlycolysis$mtDNACopyNumberRank, method='spearman') # positive - the stronger destryoed the glycilitic pathway, the less mtDNA. The best correlation!
par(mfrow=c(2,2), pch = 2)

plot(DivisionGlycolysis$MedianCovTumToNor,DivisionGlycolysis$mtDNACopyNumberRank, pch = '.')
text(DivisionGlycolysis$MedianCovTumToNor,DivisionGlycolysis$mtDNACopyNumberRank,DivisionGlycolysis$CancerTypeAndTissues)

plot(DivisionGlycolysis$MedianCovTumToNor,DivisionGlycolysis$Glycolysis, pch = '.')
text(DivisionGlycolysis$MedianCovTumToNor,DivisionGlycolysis$Glycolysis,DivisionGlycolysis$CancerTypeAndTissues)

plot(DivisionGlycolysis$NumOfCellDivPerLife,DivisionGlycolysis$Glycolysis, pch = '.')
text(DivisionGlycolysis$NumOfCellDivPerLife,DivisionGlycolysis$Glycolysis,DivisionGlycolysis$CancerTypeAndTissues)

# cor.test(DivisionGlycolysis$Glycolysis,DivisionGlycolysis$mtDNACopyNumberRank, method='spearman') # negative - the stronger destryoed the glycilitic pathway, the less mtDNA
boxplot(DivisionGlycolysis[DivisionGlycolysis$mtDNACopyNumberRank <= 7,]$Glycolysis,DivisionGlycolysis[DivisionGlycolysis$mtDNACopyNumberRank > 7,]$Glycolysis, names = c('lowmtDNA','highmtDNA'), ylab = 'Glycolysis Disturbance')
wilcox.test(DivisionGlycolysis[DivisionGlycolysis$mtDNACopyNumberRank <= 7,]$Glycolysis,DivisionGlycolysis[DivisionGlycolysis$mtDNACopyNumberRank > 7,]$Glycolysis) # 0.1

## HOW TO CLUSTER CANCERS??

par(mfrow=c(1,1))
row.names(DivisionGlycolysis) = DivisionGlycolysis$CancerTypeAndTissues
DIST = dist(DivisionGlycolysis[,4:8])
HCLUST = hclust(DIST, method = "complete")
plot(HCLUST)

#### DIFFERENCE BETWEEN CANCERS:

par(mfrow=c(2,2))

# div rate versus T>C
ALL_T_C = ALL[ALL$Subs == 'T_C',]  # 936
ALL_T_C$T_C = 1
ALL_not_T_C = ALL[ALL$Subs != 'T_C',]  # 2412
ALL_not_T_C$T_C = 0
ALL = rbind(ALL_T_C,ALL_not_T_C)

AGG = aggregate(ALL$T_C, by = list(ALL$CancerTissue, ALL$NumOfCellDivPerLife), FUN = mean)
names(AGG) = c('CancerTissue','NumOfCellDivPerLife','FrT_C')
cor.test(AGG$NumOfCellDivPerLife,AGG$FrT_C, method = 'spearman')
plot(AGG$NumOfCellDivPerLife,AGG$FrT_C, pch = '.')
text(AGG$NumOfCellDivPerLife,AGG$FrT_C,AGG$CancerTissue)

# div rate versus T>C
ALL_T_C = ALL[ALL$Subs == 'G_A',]  # 936
ALL_T_C$T_C = 1
ALL_not_T_C = ALL[ALL$Subs != 'G_A',]  # 2412
ALL_not_T_C$T_C = 0
ALL = rbind(ALL_T_C,ALL_not_T_C)

AGG = aggregate(ALL$T_C, by = list(ALL$CancerTissue, ALL$NumOfCellDivPerLife), FUN = mean)
names(AGG) = c('CancerTissue','NumOfCellDivPerLife','FrT_C')
cor.test(AGG$NumOfCellDivPerLife,AGG$FrT_C, method = 'spearman')
plot(AGG$NumOfCellDivPerLife,AGG$FrT_C, pch = '.')
text(AGG$NumOfCellDivPerLife,AGG$FrT_C,AGG$CancerTissue)

#### RANK ALL TISSUES



#### LOGISTIC REGRESSIONS
######## logistic regression: T>C as a function of CancerGlycolisis and VAF

#tissue =c('PRAD','LUNG','COAD','BRCA','KIRC','BLAD','CESC','CHOL','COADREAD','ESCA','GBM','HNSC','THCA','THYM','STAD','SKCM','SARC','READ','PCPG','PAAD','LUSC','LUAD','LIHC','KIRP','KICH','UCEC','MSKCCTvN')
#Glycolysis=c(12.38978,21.35984,59.51087,129.9068,99.2046,12.13424,14.89038,49.08608,26.55547,7.960641,3.820559,60.00785,46.08075,20.28906,19.16362,3.914527,10.47275,25.35386,43.00674,0.5165879,64.19461,25.59793,14.61152,77.26818,74.4774,29.7614,3.739063)
#OxidativePhosphorylation =c(15.76979,9.893736,14.81407,14.28531,249.1245,2.216577,3.199398,12.84608,7.390275,8.302373,6.407421,13.32523,4.457943,0.9468787,13.56102,2.319901,1.623032,11.08796,40.05531,0.3871189,26.50188,12.20636,13.09163,166.9113,51.49203,14.65809,2.185885)

#Glycolysis = data.frame(tissue,Glycolysis,OxidativePhosphorylation)
#plot(Glycolysis$Glycolysis,Glycolysis$OxidativePhosphorylation)
#ALL = merge(ALL,Glycolysis, by = 'tissue')  # 3348

## T>C: both coefficients are positive
TomacettiLifetimeIncidence = c(0.0041,0.3,0.0052,0.048,1,0.5,0.0003,0.035,0.001938,0.0028,0.00219,0.0138,0.07935,0.0071,0.071,0.0045,0.081,0.00011,0.0203,0.00035,0.00004,0.0000302,0.00022,0.00003,0.000411,0.013589,0.000194,0.0007,0.0037,0.01026,0.000324)
length(TomacettiLifetimeIncidence)
TomacettiStemCellDivPerLifetime = c(960,608,960,5840,5840,5840,1947,1947,1390,47,0,1720,1720,88,88,5.6,5.6,0,199,5,5,5,5,5,0,80,80,2920,463,7,7)
length(TomacettiStemCellDivPerLifetime)
cor.test(TomacettiLifetimeIncidence,TomacettiStemCellDivPerLifetime, method= 'spearman') # positive: 0.0003269, 0.6033422
plot(TomacettiLifetimeIncidence,TomacettiStemCellDivPerLifetime) # dev.off()

ALL_T_C = ALL[ALL$Subs == 'T_C',]  # 936
ALL_T_C$T_C = 1
ALL_not_T_C = ALL[ALL$Subs != 'T_C',]  # 2412
ALL_not_T_C$T_C = 0
ALL = rbind(ALL_T_C,ALL_not_T_C)
par(mfrow=c(1,1), cex = 2)
glm_1 <-glm(ALL$T_C ~ ALL$Glycolysis + ALL$TumorVarFreq, family = binomial())  # total number of mutations? total disruption?
summary(glm_1)
glm_2 <-glm(ALL$T_C ~ ALL$OxidativePhosphorylation  + ALL$TumorVarFreq, family = binomial())
summary(glm_2) # nothing
glm_3 <-glm(ALL$T_C ~ scale(ALL$NumOfCellDivPerLife) + scale(ALL$TumorVarFreq), family = binomial())
summary(glm_3)  # significant!!!

glm_5 <-glm(ALL$T_C ~ scale(ALL$LifeTimeRisk) + scale(ALL$TumorVarFreq), family = binomial())
summary(glm_5)  
glm_5 <-glm(ALL$T_C ~ scale(ALL$LifeTimeRisk), family = binomial())
summary(glm_5)  

## не правильно !!!! чем больше делится рак, тем больше в нем мутаций и в том числе - T>C. Нам надо добавить дамми переменные - типы раков тогда. Или random variable!!!!!!


glm_6 <-glm(ALL$T_C ~ ALL$FrequentCancer + ALL$TumorVarFreq, family = binomial())
summary(glm_6)  


glm_6 <-glm(ALL$T_C ~ scale(ALL$LifeTimeRisk) + scale(ALL$NumOfCellDivPerLife), family = binomial())
summary(glm_6)  


ALL$prob <- predict(glm_1, type = 'response') # sqrt(P-EN/pi)
summary(ALL$prob) #  0.2455  0.2544  0.2717  0.2796  0.2940  0.3713
ALL$area = (ALL$prob - min(ALL$prob))/(max(ALL$prob)-min(ALL$prob)) # normilaze data to 0-1 range
summary(ALL$area)
plot(ALL[ALL$T_C == 0,]$Glycolysis ~ ALL[ALL$T_C == 0,]$TumorVarFreq, pch = '', xlim = c(0,1.05), ylim = c(0,140), xlab = 'VAF', ylab = 'Glykolisis', main = 'T>C') #
symbols(ALL[ALL$T_C == 0,]$TumorVarFreq, ALL[ALL$T_C == 0,]$Glycolysis, circles = ALL[ALL$T_C == 0,]$area, add = TRUE, fg = "gray")
symbols(ALL[ALL$T_C == 1,]$TumorVarFreq, ALL[ALL$T_C == 1,]$Glycolysis, circles = ALL[ALL$T_C == 1,]$area, add = TRUE, fg = "red")

## G>A: both coefficients are negative
ALL_T_C = ALL[ALL$Subs == 'G_A',]  # 936
ALL_T_C$T_C = 1
ALL_not_T_C = ALL[ALL$Subs != 'G_A',]  # 2412
ALL_not_T_C$T_C = 0
ALL = rbind(ALL_T_C,ALL_not_T_C)
par(mfrow=c(1,1), cex = 2)
glm_1 <-glm(ALL$T_C ~ ALL$Glycolysis + ALL$TumorVarFreq, family = binomial())
summary(glm_1)
glm_2 <-glm(ALL$T_C ~ ALL$OxidativePhosphorylation + ALL$TumorVarFreq, family = binomial())
summary(glm_2)
ALL$prob <- predict(glm_1, type = 'response') # sqrt(P-EN/pi)
summary(ALL$prob) #  0.2455  0.2544  0.2717  0.2796  0.2940  0.3713
ALL$area = (ALL$prob - min(ALL$prob))/(max(ALL$prob)-min(ALL$prob)) # normilaze data to 0-1 range
summary(ALL$area)
plot(ALL[ALL$T_C == 0,]$Glycolysis ~ ALL[ALL$T_C == 0,]$TumorVarFreq, pch = '', xlim = c(0,1.05), ylim = c(0,140), xlab = 'VAF', ylab = 'Glykolisis', main = 'G>A') #
symbols(ALL[ALL$T_C == 0,]$TumorVarFreq, ALL[ALL$T_C == 0,]$Glycolysis, circles = ALL[ALL$T_C == 0,]$area, add = TRUE, fg = "gray")
symbols(ALL[ALL$T_C == 1,]$TumorVarFreq, ALL[ALL$T_C == 1,]$Glycolysis, circles = ALL[ALL$T_C == 1,]$area, add = TRUE, fg = "red")

### T>C are more common among early mutations (normal), while G>A are more common among late (cancer specific). Why I don't see it in previous analyses with VAF?
### - it might reflect just the same process. 


LowGlycolysis = DivisionGlycolysis[!is.na(DivisionGlycolysis$Glycolysis) & DivisionGlycolysis$Glycolysis <= median(DivisionGlycolysis[!is.na(DivisionGlycolysis$Glycolysis),]$Glycolysis),]$CancerType
HighGlycolysis = DivisionGlycolysis[!is.na(DivisionGlycolysis$Glycolysis) & DivisionGlycolysis$Glycolysis > median(DivisionGlycolysis[!is.na(DivisionGlycolysis$Glycolysis),]$Glycolysis),]$CancerType

LowMtDNA = DivisionGlycolysis[!is.na(DivisionGlycolysis$mtDNACopyNumberRank) & DivisionGlycolysis$mtDNACopyNumberRank <= median(DivisionGlycolysis[!is.na(DivisionGlycolysis$mtDNACopyNumberRank),]$mtDNACopyNumberRank),]$CancerType
HighMtDNA = DivisionGlycolysis[!is.na(DivisionGlycolysis$mtDNACopyNumberRank) & DivisionGlycolysis$mtDNACopyNumberRank > median(DivisionGlycolysis[!is.na(DivisionGlycolysis$mtDNACopyNumberRank),]$mtDNACopyNumberRank),]$CancerType

LowDiv =  DivisionGlycolysis[!is.na(DivisionGlycolysis$NumOfCellDivPerLife) & DivisionGlycolysis$NumOfCellDivPerLife <= median(DivisionGlycolysis[!is.na(DivisionGlycolysis$NumOfCellDivPerLife),]$NumOfCellDivPerLife),]$CancerType
HighDiv = DivisionGlycolysis[!is.na(DivisionGlycolysis$NumOfCellDivPerLife) & DivisionGlycolysis$NumOfCellDivPerLife >  median(DivisionGlycolysis[!is.na(DivisionGlycolysis$NumOfCellDivPerLife),]$NumOfCellDivPerLife),]$CancerType

# HighDiv + LowMtDNA + LowGlycolysis => mouse
Mouse = c(as.character(HighDiv),as.character(LowMtDNA),as.character(LowGlycolysis))
M = data.frame(table(Mouse))
Mouse.Cancer = M[M$Freq >=2,]$Mouse # HNSC LIHC READ SKCM STAD

Elephant = c(as.character(LowDiv),as.character(HighMtDNA),as.character(HighGlycolysis))
E = data.frame(table(Elephant))
Elephant.Cancer = E[E$Freq >=2,]$Elephant # COAD KICH LUAD LUSC THCA UCEC

Unknown = intersect(Mouse,Elephant) # 

ALL.M = ALL[ALL$CancerType %in% Mouse.Cancer,]  # 609
ALL.E = ALL[ALL$CancerType %in% Elephant.Cancer,] # 981

TsTv.M = nrow(ALL.M[ALL.M$Subs %in% VecOfTransitionSubstitutions,])/nrow(ALL.M[ALL.M$Subs %in% VecOfTransversionSubstitutions,]) # 12
TsTv.E = nrow(ALL.E[ALL.E$Subs %in% VecOfTransitionSubstitutions,])/nrow(ALL.E[ALL.E$Subs %in% VecOfTransversionSubstitutions,]) # 9.66


par(mfcol = c(2,3))

            
for (methods in 1:3)
{ # methods = 1
  if (methods == 1)  {C = ALL; title = paste('ALL: ',nrow(C),sep = '') }
  if (methods == 2)  {C = ALL[ALL$mtDB == 'notinDB',]; title = paste('NEW: ',nrow(C),sep = '') } # new only (5609)
  if (methods == 3)  {C = ALL[ALL$mtDB != 'notinDB',];  title = paste('OLD: ',nrow(C),sep = '')} # old only (2002)
  
  table(C$Subs)
  C$NumberOfMut = 1
  
  VecOfTransitionSubstitutions = c('T_C')  # c('A_G','G_A','C_T','T_C') # all tr
  VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # ALL Tv 
  
  MUT_TR = C[C$Subs %in% VecOfTransitionSubstitutions,]
  MUT_TV = C[C$Subs %in% VecOfTransversionSubstitutions,]
  
  AGG1 = aggregate(MUT_TR$NumberOfMut, by = list(MUT_TR$Tier2), FUN = sum); names(AGG1) = c('CancerTypes','Tr')
  AGG2 = aggregate(MUT_TV$NumberOfMut, by = list(MUT_TV$Tier2), FUN = sum); names(AGG2) = c('CancerTypes','Tv')
  TrTv = merge(AGG1,AGG2, by = 'CancerTypes')
  TrTv$TrTv = TrTv$Tr/TrTv$Tv
  TrTv = TrTv[order(TrTv$TrTv),]
  TrTv
  VecOfTissues = TrTv$CancerTypes 
  # Kidney          Myeloid         Bladder         CNS             Stomach         Esophagus       Colon/Rectum    Lymphoid        Head/Neck       Cervix          Breast         
  # Pancreas        Liver           Biliary         Thyroid         Bone/SoftTissue Ovary           Uterus          Skin            Lung            Prostate
  
  TrTv = merge(TrTv,Cells)
  TrTv = TrTv[order(TrTv$TrTv),]
  str(TrTv)
  
  plot(TrTv$NumOfCellDivPerLife,TrTv$TrTv, pch = '', ylab = 'Ts/Tv', xlab = 'Number Of Divisoins Of Each Stemm Cell Per Lifetime', main = title, xlim = c(-500, 6500)); 
  cor.test(TrTv$TrTv,TrTv$NumOfCellDivPerLife, method = 'spearman') # -0.6804434, p-value = 0.01048
  text(TrTv$NumOfCellDivPerLife,TrTv$TrTv,TrTv$CancerTypes)
  boxplot(TrTv[TrTv$TrTv < median(TrTv$TrTv),]$NumOfCellDivPerLife,TrTv[TrTv$TrTv > median(TrTv$TrTv),]$NumOfCellDivPerLife,names =c('LowTsTv','HighTsTv'), ylab = 'Number Of Divisoins Of Each Stemm Cell Per Lifetime')
}









cor.test(Glyc$AllTsTv,Glyc$Glycolysis, method = 'spearman') # p = 0.014, rho 0.5706914
plot(Glyc$AllTsTv,Glyc$Glycolysis,pch = '.')
text(Glyc$AllTsTv,Glyc$Glycolysis,Glyc$Tissue, cex = )

cor.test(Glyc$LateTsTv,Glyc$Glycolysis, method = 'spearman')  # p = 0.003583, -0.6487607
cor.test(Glyc$EarlyTsTv,Glyc$Glycolysis, method = 'spearman') # p = 0.05017,  -0.4679755 

Glyc = Glyc[order(Glyc$Glycolysis),]
boxplot(Glyc$AllTsTv[1:5],Glyc$AllTsTv[14:18], names = c('oxidative','glycolitic'), ylab = 'Ts/Tv')
wilcox.test(Glyc$AllTsTv[1:5],Glyc$AllTsTv[14:18], alternative = 'greater')

