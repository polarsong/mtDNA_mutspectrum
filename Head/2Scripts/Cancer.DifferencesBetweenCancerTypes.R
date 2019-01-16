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
  
  ######### division rate: Cancer number of divisions (from Table S1 one before last column Tomasetti and Vogelstein 2015 science) and T>C fraction (from bioarchive Campbell 2018)
  
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
  summary(LifeRisk$LifeTimeRisk)  # 
  LifeRisk$FrequentCancer = 0
  for (i in 1:nrow(LifeRisk)) 
    {
    if (LifeRisk$LifeTimeRisk[i] >= 0.023) {LifeRisk$FrequentCancer[i] = 1}
    }   
  hist(LifeRisk$LifeTimeRisk, breaks = 30)
  ALL = merge(ALL, LifeRisk, by = 'CancerTissue')
  
  ##################### NumOfCellDivPerLife
  Cells = data.frame(matrix( c("Thyroid","Head/Neck","Ovary","Esophagus","Pancreas","Skin","Lung","Liver","Colon/Rectum","Lymphoid","CNS","Bone/SoftTissue","Myeloid",                             
                               7,1720,0,1390,80,199,5.6,88,5840,960,0,5,960), ncol = 2))
  names(Cells) = c('CancerTissue','NumOfCellDivPerLife')
  Cells$NumOfCellDivPerLife = as.numeric(as.character(Cells$NumOfCellDivPerLife))
  ALL = merge(ALL,Cells, all.x = TRUE)
  
  #### turnover of cells (number of days)
  # http://book.bionumbers.org/how-quickly-do-different-cells-in-the-body-replace-themselves/
  # It has become apparent that adult stem cells not only reside and function in highly regenerative tissues like the bone marrow, intestine and epidermis where they produce a steady supply of differentiated cells to maintain their respective tissues but are also found in
  # tissues of low cell turnover, such as neural, liver, prostate and pancreas. In these tissues, adult stem cells function in maintaining tissue homeostasis by replenishing functional tissue cells lost by apoptosis.
  # Bladder: [200]
    # Urothelial stem cells localize in the basal cell layer and can generate all types of urothelial cells: Under homeostatic conditions, the adult urothelium is quiescent with one of the slowest turnover rates amongst mammalian epithelia (Hicks, 1975; Jost and Potten, 1986; Jost, 1989), however, there is a dramatic upregulation in urothelial proliferation in response to injury, which results in complete restoration of newly differentiated superficial cells within 72 hours # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5354960/
    # The rate of turnover of the cells in the normal transitional epithelium is slow but continuous: it takes more than 200 days for the progeny of a USC to mature and eventually be shed.      https://reproductivesciences.wustl.edu/laboratories/mysorekar-lab/urothelial-stem-cell-regeneration/
  # Bone/SoftTissue: 360 / 0.067 = [5373] Finally, bone stem cells divide every ~15 years (Tomacetti 2015) 15*360 = 5400 days; 
  # Breast tissue: an average turnover of [84.5] days (Tomacetti 2017)
  # Biliary: [400]
      # The biliary tract, (biliary tree or biliary system) refers to the liver, gall bladder and bile ducts, and how they work together to make, store and secrete bile.
      # don't know = took the same as liver
  # Cervix: [6] days http://book.bionumbers.org/how-quickly-do-different-cells-in-the-body-replace-themselves/
  # Lymphoid: [30] Hematopoietic stem cells have ben estimated to divide every ~15 (35), 17 (36), 30 (37, 38), 57 (39) days. Thus, we assume they divide every 30 days.  (Tomasetti 2015, Chronic lymphocytic leukemia)
  # Myeloid:  [30] Hematopoietic stem cells have ben estimated to divide every ~15 (35), 17 (36), 30 (37, 38), 57 (39) days. Thus, we assume they divide every 30 days.  (Tomasetti 2015, Acute myeloid leukemia)
  # Colon/Rectum: [5] Thus, we estimate that the total number of stem cells in the large intestine is ~1.5·107 ·15~ 2·108 and that stem cells divide, on average, every 5 days (Tomasetti 2015, Colorectal adenocarcinoma, Colorectal adenocarcinoma in FAP patients, Colorectal adenocarcinoma in patients with HNPCC)
  # Prostate: 360/3 = [120]: Thus, we estimated that normal prostatic tissue divides at 0.83/4.4=18.9% of that of prostate cancer cells. As the division rate in prostatic cancer cells is every 23.1 days or 365/23.1 =15.8 times a year (48), the normal prostate turnover is 15.8· 0.189 = 2.986, i.e., 3 times a year (Tomacetti 2017)
  # Esophagus: [11] days - The epithelial turnover rate of the healthy human esophagus is ~11 days (50). (Tomacetti 2017)
  # Stomach: [5.5]
    # 2-9 days http://book.bionumbers.org/how-quickly-do-different-cells-in-the-body-replace-themselves/
    # Stomach Adenocarcinoma
  # CNS: 360/0 = [Inf~10000]  virtually no cell division after birth (Tomacetti 2015)
  # Head/Neck: 360/21.5 = [16] days. Stem cells in the oral mucosa divide every 14 to 20 days (Tomacetti 2015)
  # Kidney: [low for example 1000]
      # The kidney is a complex organ that is derived from interactions between two distinct embryonic appendages and has more than 24 distinct mature cells. Therefore, the adult mammalian kidney is likely to have multiple stem cells and niches. Here, we will discuss some of the most promising kidney stem cell candidates and the evidence supporting their role in kidney regeneration. As isolation of stem cells using selective culture conditions has several limitations, such as contamination by blood cells and dedifferentiation of mature cells under culture conditions, we have not discussed this method in the present review. We have also not discussed the role of exogenous administration and paracrine effects of stem cells in kidney regeneration.
      # While the magnitude of renal cell turnover is lower than other organs, homeostatic mechanism(s) are still needed to maintain kidney functionality. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3665577/
      # Organ complexity and low turnover of renal cells has made stem cell identification difficult and lead to the investigation of multiple possible populations. https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-017-0339-7
  # Liver: 360/0.9125 ~ [400] days. Hepatic stem cells. Cell turnover has been estimated to be very low, with estimates ranging from 300 to 500 days (81).(Tomacetti 2015, Hepatocellular carcinoma)
  # Lung: 360/0.07 ~ [5143] days. Finally, the stem cell division rate for lung tissue (specifically the alveoli) has been estimated to be 7% per year. (Tomacetti 2015, Lung adenocarcinoma) 
  # Ovary: 360/0 = [Inf ~ 11000] There does not appear to be any oocyte divisions after week 15 of embryogenesis. (Tomacetti 2015, Ovarian germ cell cancer) 
  # Pancreas: [360] days. These stem cells are estimated to divide about once a year (Tomacetti 2015, Pancreatic ductal adenocarcinoma) 
  # Skin: [147] days. Melanocytes have been estimated to divide every ~147 days (120), i.e. 2.48 times per year. (Tomacetti 2015, Melanoma)
  # Thyroid: 360/0.087 = [4138] days. Cell turnover has been estimated to be very low, with stem cell divisions occurring once every 8.5 – 14.4 years (Tomacetti 2015, Thyroid papillary and follicular carcinoma)
  # Uterus: [4]
    # UCEC - Uterine Corpus Endometrial Carcinoma, 
    # Human endometrium is the mucosal lining of the uterus, which undergoes more than 400 cycles of regeneration, differentiation, and shedding during a woman’s reproductive life [14]. Each month, 5–10 mm of new endometrial mucosa grows during the first 4–10 days of the menstrual cycle from the residual basalis layer (0.5–1 mm thick) to generate a new functionalis layer into which the embryo subsequently implants (Fig. 1A). Within 48 hours of endometrial shedding, regeneration begins with rapid repair/re-epithelialization of the endometrial surface to cover the exposed basalis surface
    # Stem cell division is modeled (Fig. 4) as minimal before menarche (once per week), rapid with menstruation (average of one division per day), and minimal after menopause (once per week). Ages at menarche or menopause are variable among women but are modeled as occurring at 12 and 52 years of age. https://www.pnas.org/content/102/49/17739
    # So it is 7 days for a half of a life and 1 day for another half => 4 days on average (KP)
  
  CancerTissue = c('Bladder','Bone/SoftTissue','Breast','Biliary','Cervix','Lymphoid','Myeloid','Colon/Rectum','Prostate','Esophagus','Stomach','CNS','Head/Neck','Kidney','Liver','Lung','Ovary','Pancreas','Skin','Thyroid','Uterus')  
  TurnOverDays = c(200,5373,84.5,200,6,30,30,5,120,11,5.5,10000,16,1000,400,5143,11000,360,147,4138,4)
  Turn = data.frame(CancerTissue,TurnOverDays)
  
  
  
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
  

  # div rate versus T>C
  ALL_T_C = ALL[ALL$Subs == 'T_C',]  # 936
  ALL_T_C$T_C = 1
  ALL_not_T_C = ALL[!ALL$Subs == 'T_C',]  # 2412
  ALL_not_T_C$T_C = 0
  ALL = rbind(ALL_T_C,ALL_not_T_C)
  
### TEST 0: lifetyme risk and Ts/Tv
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
dev.off()






















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
dev.off()








cor.test(Glyc$AllTsTv,Glyc$Glycolysis, method = 'spearman') # p = 0.014, rho 0.5706914
plot(Glyc$AllTsTv,Glyc$Glycolysis,pch = '.')
text(Glyc$AllTsTv,Glyc$Glycolysis,Glyc$Tissue, cex = )

cor.test(Glyc$LateTsTv,Glyc$Glycolysis, method = 'spearman')  # p = 0.003583, -0.6487607
cor.test(Glyc$EarlyTsTv,Glyc$Glycolysis, method = 'spearman') # p = 0.05017,  -0.4679755 

Glyc = Glyc[order(Glyc$Glycolysis),]
boxplot(Glyc$AllTsTv[1:5],Glyc$AllTsTv[14:18], names = c('oxidative','glycolitic'), ylab = 'Ts/Tv')
wilcox.test(Glyc$AllTsTv[1:5],Glyc$AllTsTv[14:18], alternative = 'greater')

