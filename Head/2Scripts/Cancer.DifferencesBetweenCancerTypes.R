###################################
###### 26.03.2018: Cancer mut Spectrum - between different cancers 
###################################

#### nucleotide content of the human genome: A T G C: 4993	3871	2159	5357
rm(list=ls(all=TRUE))
  
  ALL = read.table("../../Body/2Derived/mtDNA_snv_Oct2016.PatientInfo.txt", head = TRUE, sep = '\t')  # 7611
  
  pdf("../../Body/4Figures/Cancer.DifferencesBetweenCancerTypes.01.pdf" , height = 30, width = 30)
  
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
MIDDLE = ALL[ALL$TurnOverDays > 30 & ALL$TurnOverDays <= 1000,]  # 4138/360
SLOW = ALL[ALL$TurnOverDays > 1000,]

# Ts/Tv is increasing, TC/Tv is increasing, TCOtherTs is increasing, TC/GA is increasing, GAOtherTs is not increasing  
TEMP = FAST
TsTv = nrow(TEMP[TEMP$Subs %in% VecOfTransitionSubstitutions,]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); TsTv # 10.02817
TCTv = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); TCTv # 2.816901
GATv = nrow(TEMP[TEMP$Subs == 'G_A',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); GATv # 5.014085
TCGA = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  == 'G_A',]); TCGA # 0.5617978
TCOtherTs = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  %in% c('C_T','G_A','A_G'),]); TCOtherTs # 0.390625
GAOtherTs = nrow(TEMP[TEMP$Subs == 'G_A',]) / nrow(TEMP[TEMP$Subs  %in% c('C_T','T_C','A_G'),]); GAOtherTs # 1 !!!!!!!!!
TEMP = MIDDLE
TsTv = nrow(TEMP[TEMP$Subs %in% VecOfTransitionSubstitutions,]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); TsTv # 12.56941
TCTv = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); TCTv # 3.827195
GATv = nrow(TEMP[TEMP$Subs == 'G_A',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); GATv # 6.634561
TCGA = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  == 'G_A',]); TCGA # 0.5768574
TCOtherTs = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  %in% c('C_T','G_A','A_G'),]); TCOtherTs # 0.4377835
GAOtherTs = nrow(TEMP[TEMP$Subs == 'G_A',]) / nrow(TEMP[TEMP$Subs  %in% c('C_T','T_C','A_G'),]); GAOtherTs # 1.1179
TEMP = SLOW
TsTv = nrow(TEMP[TEMP$Subs %in% VecOfTransitionSubstitutions,]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); TsTv # 14.49383
TCTv = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); TCTv # 4.555556
GATv = nrow(TEMP[TEMP$Subs == 'G_A',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]); GATv # 7.518519
TCGA = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  == 'G_A',]); TCGA # 0.6059113
TCOtherTs = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  %in% c('C_T','G_A','A_G'),]); TCOtherTs # 0.4583851
GAOtherTs = nrow(TEMP[TEMP$Subs == 'G_A',]) / nrow(TEMP[TEMP$Subs  %in% c('C_T','T_C','A_G'),]); GAOtherTs # 1.077876

####### bootstrep to perturb these numbers to draw three boxplots!!!!
VecFastTsTv = c()
VecFastTCTv = c()
VecFastGATv = c()
VecFastTCOtherTs = c()
VecFastGAOtherTs = c()
VecFastTCGA = c()
for (i in 1:1000) {
  TEMP = FAST[sample(nrow(FAST),nrow(FAST)/2),]; 
  TsTv = nrow(TEMP[TEMP$Subs %in% VecOfTransitionSubstitutions,]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,])
  TCTv = nrow(TEMP[TEMP$Subs  == 'T_C',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]) 
  GATv = nrow(TEMP[TEMP$Subs  == 'G_A',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]) 
  TCOtherTs = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  %in% c('C_T','G_A','A_G'),]) 
  GAOtherTs = nrow(TEMP[TEMP$Subs == 'G_A',]) / nrow(TEMP[TEMP$Subs  %in% c('C_T','T_C','A_G'),]) 
  TCGA = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  == 'G_A',]) 
  VecFastTsTv = c(VecFastTsTv,TsTv)
  VecFastTCTv = c(VecFastTCTv,TCTv) 
  VecFastGATv = c(VecFastGATv,GATv) 
  VecFastTCOtherTs = c(VecFastTCOtherTs,TCOtherTs)
  VecFastGAOtherTs = c(VecFastGAOtherTs,GAOtherTs)
  VecFastTCGA = c(VecFastTCGA,TCGA)
                  }

VecMiddleTsTv = c()
VecMiddleTCTv = c()
VecMiddleGATv = c()
VecMiddleTCOtherTs = c()
VecMiddleGAOtherTs = c()
VecMiddleTCGA = c()
for (i in 1:1000) {
  TEMP = MIDDLE[sample(nrow(MIDDLE),nrow(MIDDLE)/2),]; 
  TsTv = nrow(TEMP[TEMP$Subs %in% VecOfTransitionSubstitutions,]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]) 
  TCTv = nrow(TEMP[TEMP$Subs  == 'T_C',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]) 
  GATv = nrow(TEMP[TEMP$Subs  == 'G_A',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]) 
  TCOtherTs = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  %in% c('C_T','G_A','A_G'),]) 
  GAOtherTs = nrow(TEMP[TEMP$Subs == 'G_A',]) / nrow(TEMP[TEMP$Subs  %in% c('C_T','T_C','A_G'),]) 
  TCGA = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  == 'G_A',]) 
  VecMiddleTsTv = c(VecMiddleTsTv,TsTv) 
  VecMiddleTCTv = c(VecMiddleTCTv,TCTv) 
  VecMiddleGATv = c(VecMiddleGATv,GATv) 
  VecMiddleTCOtherTs = c(VecMiddleTCOtherTs,TCOtherTs)
  VecMiddleGAOtherTs = c(VecMiddleGAOtherTs,GAOtherTs)
  VecMiddleTCGA = c(VecMiddleTCGA,TCGA)
                  }

VecSlowTsTv = c()
VecSlowTCTv = c()
VecSlowGATv = c()
VecSlowTCOtherTs = c()
VecSlowGAOtherTs = c()
VecSlowTCGA = c()
for (i in 1:1000) {
  TEMP = SLOW[sample(nrow(SLOW),nrow(SLOW)/2),]; 
  TsTv = nrow(TEMP[TEMP$Subs %in% VecOfTransitionSubstitutions,]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]) 
  TCTv = nrow(TEMP[TEMP$Subs  == 'T_C',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]) 
  GATv = nrow(TEMP[TEMP$Subs  == 'G_A',]) / nrow(TEMP[TEMP$Subs %in% VecOfTransversionSubstitutions,]) 
  TCOtherTs = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  %in% c('C_T','G_A','A_G'),]) 
  GAOtherTs = nrow(TEMP[TEMP$Subs == 'G_A',]) / nrow(TEMP[TEMP$Subs  %in% c('C_T','T_C','A_G'),]) 
  TCGA = nrow(TEMP[TEMP$Subs == 'T_C',]) / nrow(TEMP[TEMP$Subs  == 'G_A',]) 
  VecSlowTsTv = c(VecSlowTsTv,TsTv) 
  VecSlowTCTv = c(VecSlowTCTv,TCTv) 
  VecSlowGATv = c(VecSlowGATv,GATv) 
  VecSlowTCOtherTs = c(VecSlowTCOtherTs,TCOtherTs)
  VecSlowGAOtherTs = c(VecSlowGAOtherTs,GAOtherTs)
  VecSlowTCGA = c(VecSlowTCGA,TCGA)
  }

par(mfrow=c(2,3))

### TC is increasing not only if normalized by Tv, but also if normalized by all other Ts and if normalized by GA!!!
boxplot(VecFastTsTv,VecMiddleTsTv,VecSlowTsTv, notch = TRUE, names = c('Fast','Middle','Slow'), outline = FALSE, ylab = 'Ts/Tv', col = rgb(0.1,0.1,0.1,0.3)) # dev.off()
boxplot(VecFastTsTv,VecFastTCTv,VecFastGATv,VecMiddleTsTv,VecMiddleTCTv,VecMiddleGATv,VecSlowTsTv,VecSlowTCTv,VecSlowGATv, notch = TRUE, names = c('FastTs','FastTC','FastGA','MiddleTs','MiddleTC','MiddleGA','SlowTs','SlowTC','SlowGA'), outline = FALSE, ylab = 'Ts/Tv', col = c(rgb(0.1,0.1,0.1,0.3),rgb(0.1,1,0.1,0.3),rgb(1,0.1,0.1,0.3))) # dev.off()
boxplot(VecFastTCOtherTs,VecMiddleTCOtherTs,VecSlowTCOtherTs, notch = TRUE, names = c('Fast','Middle','Slow'), outline = FALSE, ylab = 'TC/OtherTs') # dev.off()
boxplot(VecFastGAOtherTs,VecMiddleGAOtherTs,VecSlowGAOtherTs, notch = TRUE, names = c('Fast','Middle','Slow'), outline = FALSE, ylab = 'GA/OtherTs') # dev.off()
boxplot(VecFastTCOtherTs,VecFastGAOtherTs,VecMiddleTCOtherTs,VecMiddleGAOtherTs,VecSlowTCOtherTs,VecSlowGAOtherTs, notch = TRUE, names = c('FastTC','FastGA','MiddleTC','MiddleGA','SlowTC','SlowGA'), outline = FALSE, ylab = 'Common Transition / Other Ts', col = c(rgb(0.1,1,0.1,0.3),rgb(1,0.1,0.1,0.3))) # dev.off()
boxplot(VecFastTCGA,VecMiddleTCGA,VecSlowTCGA, notch = TRUE, names = c('Fast','Middle','Slow'), outline = FALSE, ylab = 'TC/GA', col = rgb(0.1,1,0.1,0.3)) # dev.off()
wilcox.test(VecFastTCGA,VecMiddleTCGA); wilcox.test(VecMiddleTCGA,VecSlowTCGA); 

####### logistic regression with all 21 Turnovers and 2 Dummy 
summary(ALL$TumorVarFreq)
quantile(ALL$TumorVarFreq,0.05)
quantile(ALL$TumorVarFreq,0.1)
quantile(ALL$TumorVarFreq,0.25) # 0.01738204
quantile(ALL$TumorVarFreq,0.5)  # 0.04539986

FAST$TurnOverDummyFast = 1; FAST$TurnOverDummySlow = 0; FAST$TurnOverRank = 1
MIDDLE$TurnOverDummyFast = 0; MIDDLE$TurnOverDummySlow = 0; MIDDLE$TurnOverRank = 2
SLOW$TurnOverDummyFast = 0; SLOW$TurnOverDummySlow = 1; SLOW$TurnOverRank = 3
ALL = rbind(FAST, MIDDLE, SLOW)

cor.test(ALL$mt_copies,ALL$TurnOverDays, method = 'spearman') # very positive - there are many mtDNA copies in slow-dividing cancers!!!!!
cor.test(ALL$mt_copies,ALL$TumorVarFreq, method = 'spearman') # a bit positive 


glm_1 <-glm(ALL$T_C ~ scale(ALL$TurnOverDays) + scale(ALL$TumorVarFreq), family = binomial())  # total number of mutations? total disruption?
summary(glm_1) # MODEL 1

glm_1 <-glm(ALL$T_C ~ scale(ALL$TurnOverDays) + scale(ALL$TumorVarFreq) + scale(ALL$Consensus_age) + scale(ALL$mt_copies), family = binomial())  # total number of mutations? total disruption?
summary(glm_1) # MODEL 1 intermediate

glm_1 <-glm(ALL$T_C ~ scale(ALL$TurnOverDays) + scale(ALL$TumorVarFreq) + scale(ALL$mt_copies), family = binomial())  # total number of mutations? total disruption?
summary(glm_1) # MODEL 1 final

glm_1a <-glm(ALL$T_C ~ scale(ALL$TurnOverDays) + scale(ALL$TurnOverDummyFast) + scale(ALL$TumorVarFreq), family = binomial())  # total number of mutations? total disruption?
summary(glm_1a) # TurnOverDummyFast is more significant than ALL$TurnOverDays => remove ALL$TurnOverDays

glm_1a <-glm(ALL$T_C ~ scale(ALL$TurnOverDummyFast) + scale(ALL$TumorVarFreq) , family = binomial())  # total number of mutations? total disruption?
summary(glm_1a) # MODEL 1A intermediate

glm_1a <-glm(ALL$T_C ~ scale(ALL$TurnOverDummyFast) + scale(ALL$TumorVarFreq) + scale(ALL$Consensus_age) + scale(ALL$mt_copies), family = binomial())  # total number of mutations? total disruption?
summary(glm_1a) # MODEL 1A intermediate

glm_1a <-glm(ALL$T_C ~ scale(ALL$TurnOverDummyFast) + scale(ALL$TumorVarFreq) + scale(ALL$mt_copies), family = binomial())  # total number of mutations? total disruption?
summary(glm_1a) # MODEL 1A final

glm_1b <-glm(ALL[ALL$TumorVarFreq > 0.01738204,]$T_C ~ scale(ALL[ALL$TumorVarFreq > 0.01738204,]$TurnOverDummyFast) + scale(ALL[ALL$TumorVarFreq > 0.01738204,]$TumorVarFreq), family = binomial())  # total number of mutations? total disruption?
summary(glm_1b)  # MODEL 1B

glm_1b <-glm(ALL[ALL$TumorVarFreq > 0.01738204,]$T_C ~ scale(ALL[ALL$TumorVarFreq > 0.01738204,]$TurnOverDummyFast) + scale(ALL[ALL$TumorVarFreq > 0.01738204,]$TumorVarFreq) + scale(ALL[ALL$TumorVarFreq > 0.01738204,]$Consensus_age) + scale(ALL[ALL$TumorVarFreq > 0.01738204,]$mt_copies), family = binomial())  # total number of mutations? total disruption?
summary(glm_1b)  # MODEL 1B Intermediate

glm_1b <-glm(ALL[ALL$TumorVarFreq > 0.01738204,]$T_C ~ scale(ALL[ALL$TumorVarFreq > 0.01738204,]$TurnOverDummyFast) + scale(ALL[ALL$TumorVarFreq > 0.01738204,]$TumorVarFreq) + scale(ALL[ALL$TumorVarFreq > 0.01738204,]$mt_copies), family = binomial())  # total number of mutations? total disruption?
summary(glm_1b)  # MODEL 1B Final mtcopies are not significant anymore!!! because all varaints are high quality!!!

glm_1b <-glm(ALL[ALL$TumorVarFreq > 0.01738204,]$T_C ~ scale(ALL[ALL$TumorVarFreq > 0.01738204,]$TurnOverDummyFast) + scale(ALL[ALL$TumorVarFreq > 0.01738204,]$TumorVarFreq), family = binomial())  # total number of mutations? total disruption?
summary(glm_1b)  # MODEL 1B Final mtcopies are not significant anymore!!! because all varaints are high quality!!!

##### now the same but only with Transitions!!!

ALL1 = ALL[ALL$Subs %in% c('T_C','C_T','G_A','A_G'),] 
glm_2 <-glm(ALL1$T_C ~ scale(ALL1$TurnOverDummyFast) + scale(ALL1$TumorVarFreq), family = binomial())  # total number of mutations? total disruption?
summary(glm_2) # both! THE FIRST BEST

glm_2 <-glm(ALL1$T_C ~ scale(ALL1$TurnOverDummyFast) + scale(ALL1$TumorVarFreq) + scale(ALL1$Consensus_age) + scale(ALL1$mt_copies), family = binomial())  # total number of mutations? total disruption?
summary(glm_2) # both! THE FIRST BEST

glm_2 <-glm(ALL1$T_C ~ scale(ALL1$TurnOverDummyFast) + scale(ALL1$TumorVarFreq)  + scale(ALL1$mt_copies), family = binomial())  # total number of mutations? total disruption?
summary(glm_2) # both! THE FIRST BEST

##### now the same but only with TC GA!!!

ALL1 = ALL1[ALL1$Subs %in% c('T_C','G_A'),] 
glm_5a <-glm(ALL1$T_C ~ scale(ALL1$TurnOverDummyFast) + scale(ALL1$TumorVarFreq) + scale(ALL1$mt_copies), family = binomial())  # total number of mutations? total disruption?
summary(glm_5a) # second only significant, THE FIRST BEST

glm_1a <-glm(ALL1$T_C ~ scale(ALL1$TurnOverDays) + scale(ALL1$TumorVarFreq), family = binomial())  # total number of mutations? total disruption?
summary(glm_1a) # second only significant, THE SECOND BEST
glm_1b <-glm(ALL1$T_C ~ scale(ALL1$TurnOverRank) + scale(ALL1$TumorVarFreq), family = binomial())  # total number of mutations? total disruption?
summary(glm_1b) # second only significant, THE THIRD BEST

#### plot glm_1

par(mfrow=c(1,1), cex = 3)

ALL$prob <- predict(glm_1aFIN, type = 'response') # sqrt(P-EN/pi)
summary(ALL$prob) #  0.2616  0.2632  0.2667  0.2785  0.2918  0.3622
ALL$area = (ALL$prob - min(ALL$prob))/(max(ALL$prob)-min(ALL$prob)) # normilaze data to 0-1 range
summary(ALL$area)
plot(log2(ALL$TurnOverDays) ~ log2(ALL$TumorVarFreq), pch = '',  xlab = 'VAF', ylim = c(1,15), xlim = c(-6.7,0.2), ylab = 'TurnOverDays', main = 'T>C') #  dev.off()
symbols(log2(ALL[ALL$T_C == 0,]$TumorVarFreq), log2(ALL[ALL$T_C == 0,]$TurnOverDays), circles = ALL[ALL$T_C == 0,]$area, add = TRUE, fg = rgb(0.1,0.1,0.1,0.1))
symbols(log2(ALL[ALL$T_C == 1,]$TumorVarFreq), log2(ALL[ALL$T_C == 1,]$TurnOverDays), circles = ALL[ALL$T_C == 1,]$area, add = TRUE, fg = rgb(1,0.0,0.0,0.2))

ALL$prob <- predict(glm_1bFIN, type = 'response') # sqrt(P-EN/pi)
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


