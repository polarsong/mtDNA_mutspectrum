### derive a table for patient-specific analysis

rm(list=ls(all=TRUE))

########################################
######### 1: DERIVE TABLE Som
########################################

############ A: read human ref seq (note that it is hg19 - not cambridge reference!):
RefSeq = read.table("../../Body/1Raw/chrM_refAllele.txt", header = FALSE, sep = '')
names(RefSeq)=c('Position','AncestralAllele')

############ B: read GTEx mutations and merge them with Ref Seq:
SomMut = read.table("../../Body/1Raw/TissueSpecificMtDNAMutationsGTEx.txt", header = TRUE, sep = '\t')
SomMut$Position = gsub('_(.*)','',SomMut$Mutation)  # gsub("_(.*)",'',READS$patient)
SomMut$DerivedAllele = gsub('(.*)_','',SomMut$Mutation)  # gsub("_(.*)",'',READS$patient)

nrow(SomMut) # 2565
Som = merge(SomMut,RefSeq, by = 'Position')
nrow(SomMut) # 2565

Som = Som[Som$DerivedAllele != Som$AncestralAllele,] # correct - in all lines Derived != Ancestral
nrow(SomMut) # 2565
Som$Substitution = paste(Som$AncestralAllele,Som$DerivedAllele, sep = '_')
table(Som$Substitution)
# A_C A_G A_T C_A C_G C_T G_A G_C G_T T_A T_C T_G 
# 43 473  49  49   6 336 728  63  57  42 693  26

########## C: read subjectID and merge patient info:
SubjId = read.table("../../Body/1Raw/GtexTaiEvgenii/pheno.csv", header = TRUE, sep = ',')
SubjId = SubjId[,c(1,6,7,8)] # keep only SRR and patient specific important info (sex, age, Sibjid)
names(SubjId)[names(SubjId) == "run"] <- "SRR" # rename run into SRR

Som = merge(Som,SubjId, by = 'SRR')

########## D: read GTEx expression data:
Expr = read.table("../../Body/1Raw/TissueSpecificMtDNACoverageGTEx.txt", header = TRUE, sep = '\t')
ExprTissues = data.frame(table(Expr$Tissue))
ExprTissues = ExprTissues[order(ExprTissues$Var1),]

# merge by SRRID (lost 2565 - 2523 observations, but got expression level (percentMito) for each sample and alternative name of the Tissue (more simple))
nrow(Som) # 2565
Som = merge(Som,Expr, by.x = 'SRR', by.y = 'SRRID', all.x = TRUE) 
nrow(Som) # 2565  - but not all of them have expression level defined  (2523 out of 2565)

TissueLongToTissueShort = unique(data.frame(Som$tissue,Som$Tissue))
# now, based on TissueLongToTissueShort and common sence I  manually create tissue vocabular:

TissuesVocab = c(
'Heart - Left Ventricle','Heart',
'Heart - Atrial Appendage','Heart',
'Cells - Transformed fibroblasts','Skin',
'Esophagus - Mucosa', 'Esophagus',
'Colon - Sigmoid', 'Colon',
'Artery - Aorta', 'Blood Vessel',
'Breast - Mammary Tissue', 'Breast',
'Cells - EBV-transformed lymphocytes', 'Blood',
'Liver', 'Liver',
'Whole Blood', 'Blood',
'Muscle - Skeletal', 'Muscle',
'Colon - Transverse', 'Colon',
'Adrenal Gland', 'Adrenal Gland',
'Vagina', 'Vagina',
'Stomach', 'Stomach',
'Pancreas', 'Pancreas',
'Nerve - Tibial', 'Nerve',
'Uterus', 'Uterus',
'Adipose - Subcutaneous', 'Adipose Tissue',
'Adipose - Visceral (Omentum)', 'Adipose Tissue', # my decision
'Skin - Not Sun Exposed (Suprapubic)', 'Skin',
'Skin - Sun Exposed (Lower leg)', 'Skin',
'Artery - Coronary', 'Blood Vessel',
'Artery - Tibial', 'Blood Vessel',
'Prostate', 'Prostate',
'Spleen', 'Spleen',
'Esophagus - Muscularis', 'Esophagus',
'Testis', 'Testis',
'Kidney - Cortex', 'Kidney',
'Lung', 'Lung',
'Ovary', 'Ovary',
'Small Intestine - Terminal Ileum', 'Small Intestine',
'Esophagus - Gastroesophageal Junction','Esophagus',
'Minor Salivary Gland', 'Salivary Gland',
'Thyroid', 'Thyroid',
'Pituitary', 'Pituitary',
'Brain - Putamen (basal ganglia)','Brain',
'Brain - Spinal cord (cervical c-1)','Brain',
'Brain - Caudate (basal ganglia)','Brain',
'Brain - Cerebellar Hemisphere','Brain',
'Brain - Substantia nigra','Brain',
'Brain - Anterior cingulate cortex (BA24)','Brain',
'Brain - Amygdala','Brain',
'Brain - Hypothalamus','Brain',
'Brain - Hippocampus','Brain',
'Brain - Cortex','Brain',
'Brain - Nucleus accumbens (basal ganglia)','Brain',
'Brain - Cerebellum','Brain',
'Brain - Frontal Cortex (BA9)','Brain'
)
TISS = data.frame(TissuesVocab[c(TRUE, FALSE)],TissuesVocab[c(FALSE, TRUE)]); names(TISS)=c('tissue','TissueShortName')
setdiff(as.character(unique(Som$tissue)),as.character(TISS$TissueLongName))
Som = merge(Som, TISS) # add TissueShortName to each 'tissue'

########### E: Cell lifespans from the table here: https://docs.google.com/document/d/1UECub1DIdmuXwPqDLK8WRcZ6uIjEQiVmHXN1_XNVvUU/edit?usp=sharing 
VecOfTissues = c('Adipose Tissue','Adrenal Gland','Bladder','Blood','Blood Vessel','Bone Marrow','Brain','Breast','Cervix Uteri','Colon','Esophagus','Heart','Kidney','Liver','Lung','Muscle','Nerve','Ovary','Pancreas','Prostate','Salivary Gland','Skin','Small Intestine','Spleen','Stomach','Testis','Thyroid','Uterus','Vagina','Pituitary')
VecOfTurnOvers = c(  2448,  455,  49,  30,  67.5,  3.2,  24637.5,  65.75,  5.7,  4.25,  10.5,  25300,  270,  363.5,  126,  5510,  24637.5,  30000,  315,  120,  60,  64,  7.05,  7.8,  1.4,  63.5,  3679.5,  13,  3.9, 25000)
All = data.frame(VecOfTissues,VecOfTurnOvers); names(All) = c('TissueShortName','TurnOverRate')
Som = merge(Som,All, by = 'TissueShortName') 

########### F: write derive for future analyses:

write.table(Som,"../../Body/2Derived/HealthySomaticHumanMutations.GTEx.PatientSpecificDerive.txt", row.names = FALSE, quote = FALSE, sep = '\t')
