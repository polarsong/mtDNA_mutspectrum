###################################
###### 26.03.2018: Cancer mut Spectrum - between different cancers 
###################################

#### nucleotide content of the human genome: A T G C: 4993	3871	2159	5357
rm(list=ls(all=TRUE))

ALL = read.table("../../Body/1Raw/mtDNA_snv_Oct2016.txt", head = TRUE, sep = '\t')  # 7611

### DERIVE NECESSARY TRAITS:
ALL$TumorVarFreq = ALL$tumor_reads2/(ALL$tumor_reads1 + ALL$tumor_reads2); summary(ALL$TumorVarFreq)  # 0.01000 0.01738 0.04540 0.20268 0.26278 0.99864
ALL$NormalVarFreq = ALL$normal_reads2/(ALL$normal_reads1 + ALL$normal_reads2); summary(ALL$NormalVarFreq) 
ALL$Subs = paste(ALL$ref,ALL$var, sep = '_'); table(ALL$Subs)
VecOfTransitionSubstitutions = c('A_G','G_A','C_T','T_C') # all tr
VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # ALL Tv
CancerType = unique(ALL$tissue); length(CancerType)    # 40
CancerTissue = unique(ALL$Tier2); length(CancerTissue) # 21
ALL$CancerTissue = ALL$Tier2
ALL$CancerType = ALL$tissue
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

Cells = data.frame(matrix( c("Thyroid","Head/Neck","Ovary","Esophagus","Pancreas","Skin","Lung","Liver","Colon/Rectum","Lymphoid","CNS","Bone/SoftTissue","Myeloid",                             
                             7,1720,0,1390,80,199,5.6,88,5840,960,0,5,960), ncol = 2))
names(Cells) = c('CancerTissue','NumOfCellDivPerLife')
Cells$NumOfCellDivPerLife = as.numeric(as.character(Cells$NumOfCellDivPerLife))
ALL = merge(ALL,Cells, all.x = TRUE)

######## glycolisis: from Andrey Yurchenko and https://www.nature.com/articles/s41467-018-07232-8#MOESM5

CancerType =c('PRAD','LUNG','COAD','BRCA','KIRC','BLAD','CESC','CHOL','COADREAD','ESCA','GBM','HNSC','THCA','THYM','STAD','SKCM','SARC','READ','PCPG','PAAD','LUSC','LUAD','LIHC','KIRP','KICH','UCEC','MSKCCTvN')
Glycolysis=c(12.38978,21.35984,59.51087,129.9068,99.2046,12.13424,14.89038,49.08608,26.55547,7.960641,3.820559,60.00785,46.08075,20.28906,19.16362,3.914527,10.47275,25.35386,43.00674,0.5165879,64.19461,25.59793,14.61152,77.26818,74.4774,29.7614,3.739063)
OxPhos = 
Glyc = data.frame(CancerType,Glycolysis)
ALL = merge(ALL,Glyc, all.x=TRUE)

####### mtDNA copies: https://www.ncbi.nlm.nih.gov/pubmed/26901439 (figure 3 gives a rank of cancer types)
CancerType = c('KIRC','BRCA','BLCA','LIHC','HNSC','ESCA','KIRP','STAD','UCEC','KICH','COAD','THCA','PAAD','PRAD','LUAD')
mtDNACopyNumberRank = c(1:15)
CopyNumber = data.frame(CancerType,mtDNACopyNumberRank)
ALL = merge(ALL,CopyNumber, all.x=TRUE)

#### TEST 1: is there correlation between division rate and glycolisis? I would think that fastly dividing are more glycolitic?
# nonsignificant positive - probably because there are many NAs

DivisionGlycolysis = ALL[grep("CancerType|CancerTissue|CancerTypeAndTissue|NumOfCellDivPerLife|Glycolysis|mtDNACopyNumberRank",colnames(ALL))]
DivisionGlycolysis = unique(DivisionGlycolysis)
cor.test(DivisionGlycolysis$NumOfCellDivPerLife,DivisionGlycolysis$Glycolysis, method='spearman') # positive = the more divisions, the more we can destroy glycolitic pathway
cor.test(DivisionGlycolysis$NumOfCellDivPerLife,DivisionGlycolysis$mtDNACopyNumberRank, method='spearman') # negative - the more divisions, the less mtDNA
cor.test(DivisionGlycolysis$Glycolysis,DivisionGlycolysis$mtDNACopyNumberRank, method='spearman') # negative - the stronger destryoed the glycilitic pathway, the less mtDNA. The best correlation!
plot(DivisionGlycolysis$NumOfCellDivPerLife,DivisionGlycolysis$Glycolysis, pch = '.')
text(DivisionGlycolysis$NumOfCellDivPerLife,DivisionGlycolysis$Glycolysis,DivisionGlycolysis$CancerTypeAndTissues)

# cor.test(DivisionGlycolysis$Glycolysis,DivisionGlycolysis$mtDNACopyNumberRank, method='spearman') # negative - the stronger destryoed the glycilitic pathway, the less mtDNA
boxplot(DivisionGlycolysis[DivisionGlycolysis$mtDNACopyNumberRank <= 7,]$Glycolysis,DivisionGlycolysis[DivisionGlycolysis$mtDNACopyNumberRank > 7,]$Glycolysis, names = c('lowmtDNA','highmtDNA'), ylab = 'Glycolysis Disturbance')
wilcox.test(DivisionGlycolysis[DivisionGlycolysis$mtDNACopyNumberRank <= 7,]$Glycolysis,DivisionGlycolysis[DivisionGlycolysis$mtDNACopyNumberRank > 7,]$Glycolysis) # 0.1

## I CAN merge subset of cancers with low mtDNA or high glucose disturbance or high division rate

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


#### pdf out
pdf("../../Body/4Figures/Cancer.DifferencesBetweenCancerTypes.01.pdf" , height = 20, width = 15)
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

### do it better = sort each variant by VAF and calculate Ts/Tv based for example on the last 20 substitutions.

######## logistic regression: T>C as a function of CancerGlycolisis and VAF

tissue =c('PRAD','LUNG','COAD','BRCA','KIRC','BLAD','CESC','CHOL','COADREAD','ESCA','GBM','HNSC','THCA','THYM','STAD','SKCM','SARC','READ','PCPG','PAAD','LUSC','LUAD','LIHC','KIRP','KICH','UCEC','MSKCCTvN')
Glycolysis=c(12.38978,21.35984,59.51087,129.9068,99.2046,12.13424,14.89038,49.08608,26.55547,7.960641,3.820559,60.00785,46.08075,20.28906,19.16362,3.914527,10.47275,25.35386,43.00674,0.5165879,64.19461,25.59793,14.61152,77.26818,74.4774,29.7614,3.739063)
OxidativePhosphorylation =c(15.76979,9.893736,14.81407,14.28531,249.1245,2.216577,3.199398,12.84608,7.390275,8.302373,6.407421,13.32523,4.457943,0.9468787,13.56102,2.319901,1.623032,11.08796,40.05531,0.3871189,26.50188,12.20636,13.09163,166.9113,51.49203,14.65809,2.185885)

Glycolysis = data.frame(tissue,Glycolysis,OxidativePhosphorylation)
plot(Glycolysis$Glycolysis,Glycolysis$OxidativePhosphorylation)
ALL = merge(ALL,Glycolysis, by = 'tissue')  # 3348

## T>C: both coefficients are positive
ALL_T_C = ALL[ALL$Subs == 'T_C',]  # 936
ALL_T_C$T_C = 1
ALL_not_T_C = ALL[ALL$Subs != 'T_C',]  # 2412
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

#plot(RL[RL$StatusRL == 'LC',]$lnBM_g,RL[RL$StatusRL == 'LC',]$lnDN_DS, pch = '', xlim = c(2,23), ylim = c(-5.5,-1.4), xlab = 'ln(Body mass(gramm))', ylab = 'ln(Kn/Ks)')
#text(RL[RL$StatusRL == 'LC',]$lnBM_g,RL[RL$StatusRL == 'LC',]$lnDN_DS, RL[RL$StatusRL == 'LC',]$NUMBER, xlim = c(0,23), ylim = c(-5.5,-1.4), xlab = 'ln(Body mass(gramm))', ylab = 'ln(Kn/Ks)', cex = 0.7, col = 'black')
#par(new = TRUE)
#plot(RL[RL$StatusRL == 'EN',]$lnBM_g,RL[RL$StatusRL == 'EN',]$lnDN_DS, pch = '', xlim = c(2,23), ylim = c(-5.5,-1.4), xlab = 'ln(Body mass(gramm))', ylab = 'ln(Kn/Ks)')
#text(RL[RL$StatusRL == 'EN',]$lnBM_g,RL[RL$StatusRL == 'EN',]$lnDN_DS, RL[RL$StatusRL == 'EN',]$NUMBER, xlim = c(0,23), ylim = c(-5.5,-1.4), xlab = 'ln(Body mass(gramm))', ylab = 'ln(Kn/Ks)', cex = 0.7, col = 'red')
#legend(20,-1.4,legend = head(RL$NUM_SP,76), cex = 0.5)
#legend(22,-1.4,legend = tail(RL$NUM_SP,75), cex = 0.5)
