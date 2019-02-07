#####  COUNT NUCLEOTIDE CONTENT CAREFULLY (BODY/2Derived/polarizedbr_data => external + More Shallow => codons, 4fold nucl, FrA,T,G,C) => barplot?
#####  normalized average MutSpec (pie charts or 12 boxplots for each class)?

rm(list=ls(all=TRUE))

VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG')
length(unique(VecOfSynFourFoldDegenerateSites)) # 32

unzip("../../Body/2Derived/POLARIZEDBR_DATA.zip", exdir= "../../Body/2Derived/")
List = list.files("../../Body/2Derived/POLARIZEDBR_DATA/")

for (i in 1:10) # length(List))
{ # i = 1
  infile = paste("../../Body/2Derived/POLARIZEDBR_DATA/",as.character(List[i]),sep='') 
  Species = gsub('\\..*','',as.character(List[i]))
  Gene = gsub(Species,'',as.character(List[i])); Gene = gsub('.POLARISED.txt','',Gene); Gene = gsub('\\.POLARISED.txt','',Gene); Gene = gsub('\\.','',Gene); 
  GeneSpecies = read.table(infile, header = TRUE)
  GeneSpecies = GeneSpecies[GeneSpecies$BranchPosition == 'External',]
  ExternalSeqsTogether = paste(GeneSpecies$MoreShallowNodeSeq,collapse = '')
  ExternalSeqsTogether = unlist(strsplit(ExternalSeqsTogether,'')) # 5700/3
  CodonsVec = c(); StartNuc = 1
  if (length(ExternalSeqsTogether)/3 == round(length(ExternalSeqsTogether)/3))  # if divide by 3 without the rest
  {
  for (j in 1:(length(ExternalSeqsTogether)/3))
    {
    CodonsVec = c(CodonsVec,paste(ExternalSeqsTogether[StartNuc : (StartNuc+2)],collapse = ''))
    StartNuc = StartNuc+3
    }
  AllCodons = length(CodonsVec)        # 1021
  CodonsVecNeutral = CodonsVec[CodonsVec %in% VecOfSynFourFoldDegenerateSites]
  NeutralCodons = length(CodonsVecNeutral) # 1900
  data.frame(table(CodonsVecNeutral))
  
  CodonsVecNeutral = gsub("CTA|GTA|TCA|CCA|ACA|GCA|CGA|GGA",'A',CodonsVecNeutral)
  CodonsVecNeutral = gsub("CTT|GTT|TCT|CCT|ACT|GCT|CGT|GGT",'T',CodonsVecNeutral)
  CodonsVecNeutral = gsub("CTG|GTG|TCG|CCG|ACG|GCG|CGG|GGG",'G',CodonsVecNeutral)
  CodonsVecNeutral = gsub("CTC|GTC|TCC|CCC|ACC|GCC|CGC|GGC",'C',CodonsVecNeutral)
  
  Line=c(Species,Gene,length(CodonsVecNeutral[CodonsVecNeutral == 'A']),length(CodonsVecNeutral[CodonsVecNeutral == 'T']),length(CodonsVecNeutral[CodonsVecNeutral == 'G']),length(CodonsVecNeutral[CodonsVecNeutral == 'C']), AllCodons, NeutralCodons)
  if (i == 1) {Final = Line}
  if (i >  1) {Final = rbind(Final,Line)}
  }
}

Final = as.data.frame(Final); names(Final)=c('Species','Gene','A','T','G','C',"NumberOfFourFoldDegenCodons",'NumberOfAllCodons')






MUT = read.table("../../Body/3Results/Mutational_spectra_in_Chordata_ML.txt", header = TRUE)
length(unique(MUT$Species)) 

Taxa = read.table("../../Body/1Raw/TaxaFromKostya.Names.stat", sep = '\t',header = FALSE) 
Taxa$Species = gsub(";.*",'',Taxa$V1); Taxa$Species = gsub(" ",'_',Taxa$Species);  
Taxa$Class = gsub(";Chordata;.*",'',Taxa$V1); Taxa$Class = gsub(".*;",'',Taxa$Class); table(Taxa$Class)
Taxa$Class = gsub('Actinopteri','Actinopterygii',Taxa$Class)
Taxa$Class = gsub("Testudines|Squamata|Crocodylia|Sphenodontia",'Reptilia',Taxa$Class)
length(unique(Taxa$Species)) # 1708
table(Taxa$Class)

SpeciesWithoutTaxonomy = setdiff(unique(MUT$Species),unique(Taxa$Species)); length(SpeciesWithoutTaxonomy)
write.table(SpeciesWithoutTaxonomy, "../../Body/2Derived/NeedTaxa.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

##### FILTER 1: to take only normal substitutions and filter out species with too high fraction (> 5%) of unnormal substitutions
VecOfNormalSubstitutions <- c('A_C','C_A',
                              'A_G','G_A',
                              'C_G','G_C',
                              'C_T','T_C',
                              'G_T','T_G',
                              'T_A','A_T')
nrow(MUT)
table(MUT$Subs)   # MANY CRAPPY SUBSTITUTIONS!!!!!!!!!!!!!!!!!! WHY?????????????????????
SP = data.frame(table(MUT$Species)); names(SP) = c('Species','NumberOfAllSubst')
SPN = data.frame(table(MUT[MUT$Subs %in% VecOfNormalSubstitutions,]$Species)); names(SPN) = c('Species','NumberOfNormalSubst')
SP = merge(SP,SPN); SP$FractionOfNormal = SP$NumberOfNormalSubst/SP$NumberOfAllSubst
hist(SP$FractionOfNormal)
summary(SP$FractionOfNormal) # how many to delete? ask to have more than 95% of substitutions as normal
SpeciesToDelete = SP[SP$FractionOfNormal <=0.95,]$Species; length(SpeciesToDelete) # 279 - delete
MUT = MUT[!MUT$Species %in% SpeciesToDelete,]
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]

##### FILTER 2: Synonymous Substitutions
nrow(MUT) # 461215
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]; nrow(MUT) # 395157

##### FILTER 3: fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG')
length(unique(VecOfSynFourFoldDegenerateSites)) # 32
nrow(MUT) # 395157
MUT = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]; nrow(MUT) # 209120

##### FILTER: Gene
table(MUT$Gene)
# ATP6   ATP8   COX1   COX2   COX3   CytB    ND1    ND2    ND3    ND4   ND4L 
# 9752    625  13233   4264   4873 109049  14314  33809   2635  14511   2055  # there are quite many genes - do I need to analyse all of them focus on CYTB? Better - all and from time to time to check that results based on CYTB only are robust!!!

MUT = merge(MUT,Taxa, all.x = TRUE)  ##### NOT ALL SPECIES HAVE TAXONOMY!!!!

##### derive class specific spectra without normalization

Equil = c()
par(mfrow = c(1,5))
VecOfClasses = c('Actinopterygii','Amphibia','Reptilia','Mammalia','Aves')
for (i in 1:length(VecOfClasses))
{ # i = 1
  Temp = MUT[MUT$Class == VecOfClasses[i],]
  NumberOfSpecies = length(unique(Temp$Species))
  title = paste(VecOfClasses[i],', N = ',NumberOfSpecies, sep = '')
  Temp$number = 1
  Agg = aggregate(Temp$number, by = list(Temp$Subs), FUN = sum); names(Agg) = c('Subs','Number')
  Agg$Number = Agg$Number/sum(Agg$Number)
  sum(Agg$Number) # 1 - 100%
  pie(Agg$Number, labels = Agg$Subs, main = title, col=rainbow(length(Agg$Subs)))
  
  ToGFromG = sum(Agg[Agg$Subs %in% c('A_G','T_G','C_G'),]$Number) /  sum(Agg[Agg$Subs %in% c('G_A','G_T','G_C'),]$Number)
  ToTFromT = sum(Agg[Agg$Subs %in% c('A_T','G_T','C_T'),]$Number) /  sum(Agg[Agg$Subs %in% c('T_A','T_G','T_C'),]$Number)
  ToAFromA = sum(Agg[Agg$Subs %in% c('T_A','G_A','C_A'),]$Number) /  sum(Agg[Agg$Subs %in% c('A_T','A_G','A_C'),]$Number)
  ToCFromC =  sum(Agg[Agg$Subs %in% c('T_C','G_C','A_C'),]$Number) / sum(Agg[Agg$Subs %in% c('C_T','C_G','C_A'),]$Number)
  
  Equil = rbind(Equil,c(i,VecOfClasses[i],ToGFromG,ToTFromT,ToAFromA,ToCFromC))
}

Equil = data.frame(Equil); names(Equil) = c('i','Classes','G','T','C','A')
# i        Classes                 G                 T                C                 A
# 1 Actinopterygii 0.594741697416974  0.88695652173913 1.64319645994664  1.04326864396436
# 2       Amphibia   0.5125284738041  1.04096989966555 1.76627712854758 0.917777777777778
# 3       Reptilia 0.532318741450068 0.869721115537849 1.69920494699647   1.0426957687153
# 4       Mammalia 0.538739462063428 0.899331180811808 1.77489626556017  1.02578796561605
# 5           Aves 0.569868001820665 0.766223612197029  1.8578811369509  1.08779093507554

### equilibrium for each species:
Equilibrium = c()
MUT$number = 1
AGG = aggregate(MUT$number, by = list(MUT$Subs,MUT$Species,MUT$Class), FUN = sum); names(AGG) = c('Subs','Species','Class','Number')
VecOfSpecies = unique(MUT$Species)
for (i in 1:length(VecOfSpecies))
  { # i = 1
  Agg = AGG[AGG$Species == VecOfSpecies[i],]
  if (sum(Agg$Number) >= 15)
    {
    ToGFromG = log2(sum(Agg[Agg$Subs %in% c('A_G','T_G','C_G'),]$Number) /  sum(Agg[Agg$Subs %in% c('G_A','G_T','G_C'),]$Number))
    ToTFromT = log2(sum(Agg[Agg$Subs %in% c('A_T','G_T','C_T'),]$Number) /  sum(Agg[Agg$Subs %in% c('T_A','T_G','T_C'),]$Number))
    ToAFromA = log2(sum(Agg[Agg$Subs %in% c('T_A','G_A','C_A'),]$Number) /  sum(Agg[Agg$Subs %in% c('A_T','A_G','A_C'),]$Number))
    ToCFromC =  log2(sum(Agg[Agg$Subs %in% c('T_C','G_C','A_C'),]$Number) / sum(Agg[Agg$Subs %in% c('C_T','C_G','C_A'),]$Number))
    
    Equilibrium = rbind(Equilibrium,c(as.character(VecOfSpecies[i]),Agg$Class[1],ToGFromG,ToTFromT,ToAFromA,ToCFromC))
    }
  }
Equilibrium = data.frame(Equilibrium); names(Equilibrium)=c('Species','Class','G','T','A','C')
Equilibrium[,3:6] <- sapply(Equilibrium[,3:6], function(x) as.numeric(as.character(x)))
par(mfrow=c(1,1))
Equilibrium = Equilibrium[Equilibrium$Class %in% VecOfClasses,]
boxplot(
Equilibrium[Equilibrium$Class == 'Actinopterygii',]$G,Equilibrium[Equilibrium$Class == 'Actinopterygii',]$T, Equilibrium[Equilibrium$Class == 'Actinopterygii',]$C, Equilibrium[Equilibrium$Class == 'Actinopterygii',]$A, 
Equilibrium[Equilibrium$Class == 'Amphibia',]$G,Equilibrium[Equilibrium$Class == 'Amphibia',]$T, Equilibrium[Equilibrium$Class == 'Amphibia',]$C, Equilibrium[Equilibrium$Class == 'Amphibia',]$A, 
Equilibrium[Equilibrium$Class == 'Reptilia',]$G,Equilibrium[Equilibrium$Class == 'Reptilia',]$T, Equilibrium[Equilibrium$Class == 'Reptilia',]$C, Equilibrium[Equilibrium$Class == 'Reptilia',]$A, 
Equilibrium[Equilibrium$Class == 'Mammalia',]$G,Equilibrium[Equilibrium$Class == 'Mammalia',]$T, Equilibrium[Equilibrium$Class == 'Mammalia',]$C, Equilibrium[Equilibrium$Class == 'Mammalia',]$A, 
Equilibrium[Equilibrium$Class == 'Aves',]$G,Equilibrium[Equilibrium$Class == 'Aves',]$T, Equilibrium[Equilibrium$Class == 'Aves',]$C, Equilibrium[Equilibrium$Class == 'Aves',]$A, 
outline = FALSE, notch = TRUE, col = c(ColG,ColT,ColC,ColA), names = rep(c('G','T','C','A'),5))
abline(h = 0, col = 'red')

dev.off()

# TO DO:
#####  ADD TAXA => MORE SPECIES
#####  COUNT NUCLEOTIDE CONTENT CAREFULLY (BODY/2Derived/polarizedbr_data => external + More Shallow => codons, 4fold nucl, FrA,T,G,C) => barplot?
#####  notmalized average MutSpec (pie charts or 12 boxplots for each class)?

####################
#### PCA FOR ALL CHORDATA
####################

rm(list=ls(all=TRUE))

### SETTINGS:
ParamQuantileNumberOfSubst = 0.25  # used 0.25 to generate FISH.txt file
ParamNumberOfZeroes = 12 # do not use this filter
ParamScaleOrNot = FALSE # TRUE # FALSE

library(dplyr)  # dplyr
library(purrr)  # install.packages("purrr") # as_vector()

user = 'Kostya'
if (user == 'Kostya') {setwd('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/')}

###### LOAD MUTATION SPECTRUM 
# METHOD = 'PARSIMONY'  
METHOD = 'MAXLIKELIHOOD'

if (METHOD == 'MAXLIKELIHOOD') {MUT = read.table('2_DERIVED/Table_fixed.txt', header = TRUE)}
if (METHOD == 'PARSIMONY')     {MUT = read.table('2_DERIVED/Table_fixed_parsymony.txt', header = TRUE)}


##### LOAD LONGEVITY OF MAMMALS
GenerTime = read.table('1_RAW/GenerationLenghtforAllMammals/GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
GenerTime$Species = gsub(' ','_',GenerTime$Scientific_name)
GenerTime = GenerTime[,grepl("Species|GenerationLength_d", names(GenerTime))]

##### LOAD ECOLOGY FOR ALL CHORDATA
AnAge = read.table('1_RAW/anage_data.txt', header = TRUE, sep = '\t')
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') # names(AnAge)
AnAge = AnAge[,grepl("Species|Female.maturity..days.|Adult.weight..g.|Maximum.longevity..yrs.|Metabolic.rate..W.|Temperature..K.",names(AnAge))]

##### LOAD TEMPERATURE AND MBR
Temp = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/EnvironmentalTemperature/TemperatureDataSet.txt', header = TRUE, sep = '\t')
Temp$Species = gsub(" ", '_', Temp$Species)
cor.test(Temp[Temp$Tax == 'fishes',]$T..oC.,Temp[Temp$Tax == 'fishes',]$Basal.metabolic.rate..watts., method = 'spearman') # nothing
cor.test(Temp[Temp$Tax == 'fishes',]$T..oC.,Temp[Temp$Tax == 'fishes',]$Body.mass..gram., method = 'spearman') # nothing
cor.test(Temp[Temp$Tax == 'fishes',]$Body.mass..gram.,Temp[Temp$Tax == 'fishes',]$Basal.metabolic.rate..watts., method = 'spearman') # positive!! and strong

##### FILTER: take only normal substitutions and subset species with high fraction (0.95) of normal substitutions
VecOfNormalSubstitutions <- c('A_C','C_A',
                            'A_G','G_A',
                            'C_G','G_C',
                            'C_T','T_C',
                            'G_T','T_G',
                            'T_A','A_T')
nrow(MUT)
table(MUT$Subs)   # MANY CRAPPY SUBSTITUTIONS!!!!!!!!!!!!!!!!!! WHY?????????????????????
MUT %>% count(Species, sort = TRUE) -> cs; hist(cs$n, breaks = 50)
MUT[MUT$Subs %in% VecOfNormalSubstitutions,] %>% count(Species, sort = TRUE) -> csN
CS = merge(cs,csN, by = 'Species', all = FALSE)
CS$FractionOfNormal = CS$n.y/CS$n.x; summary(CS$FractionOfNormal)
SpeciesWithHighFractionOfNomralSubstitutions = as.character(CS[CS$FractionOfNormal > 0.95,]$Species); 
length(SpeciesWithHighFractionOfNomralSubstitutions)
length(unique(SpeciesWithHighFractionOfNomralSubstitutions)); head(SpeciesWithHighFractionOfNomralSubstitutions)
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions & MUT$Species %in% SpeciesWithHighFractionOfNomralSubstitutions,]; nrow(MUT)


##### FILTER: species with many substitutions (more than the lower quartile)
MUT %>% count(Species, sort = TRUE) -> cs
hist(cs$n)
summary(cs$n)
quantile(cs$n, .05) # 4 substitutions
quantile(cs$n, .1)  # 6 substitutions
quantile(cs$n, .15) # 9 substitutions
quantile(cs$n, .25) # 15 substitutions
quantile(cs$n, .3)  # 20 substitutions
ListOfSpeciesWithManySubst <- filter(cs, n >= quantile(n, ParamQuantileNumberOfSubst)) %>% select(Species) %>% as_vector()
MUT <- MUT %>% filter(Species %in% ListOfSpeciesWithManySubst)

##### FILTER: class
MUT = merge(MUT, Taxa, by = 'Species')  # some species disappears!!!! they are not in taxa file!
MUT = MUT[MUT$Class == 'Actinopteri',]  # Actinopteri   Mammalia

###### DERIVE MUTATIONAL SPECTRUM:
### NORMALIZATION of the 'NumberOfSynMutPerSpecies' by ancestral nucleotide count in the third position of four-fold synonymous substitutions:
SpeciesInMut = unique(MUT$Species); length(SpeciesInMut) # 1172
NUC = read.table('2_DERIVED/ATGC_counts_in_SYN_codons_with_full_gene.txt', header = TRUE)
NUC$Gene = gsub("(.*)\\.",'',NUC$Species)
NUC$Species = gsub("\\.(.*)",'',NUC$Species)
SpeciesInNuc = unique(NUC$Species); length(SpeciesInNuc)
MUT = merge(MUT,NUC, by = c("Species","Gene"))  # compare CountA.x and CountA.y  - they should be identical.
nrow(MUT) # 
setdiff(SpeciesInMut,unique(MUT$Species))  # 

EXTRACT = function(x) {first = unlist(strsplit(as.character(x),'_'))[1]; return(first);}; MUT$AncestralNuc = apply(as.matrix(MUT$Subs), 1, EXTRACT)
MUT$NumberOfSynMutPerSpecies = 1
MUT_A = MUT[MUT$AncestralNuc == 'A',]; MUT_T = MUT[MUT$AncestralNuc == 'T',]; MUT_G = MUT[MUT$AncestralNuc == 'G',]; MUT_C = MUT[MUT$AncestralNuc == 'C',] # 64145+123587+97657+128195=413584 
MUT_A$NumberOfSynMutPerSpecies = MUT_A$NumberOfSynMutPerSpecies/MUT_A$CountA_Syn;
MUT_T$NumberOfSynMutPerSpecies = MUT_T$NumberOfSynMutPerSpecies/MUT_T$CountT_Syn;
MUT_G$NumberOfSynMutPerSpecies = MUT_G$NumberOfSynMutPerSpecies/MUT_G$CountG_Syn;
MUT_C$NumberOfSynMutPerSpecies = MUT_C$NumberOfSynMutPerSpecies/MUT_C$CountC_Syn;
MUT = rbind(MUT_A,MUT_T,MUT_G,MUT_C)

### COUNT THE TOTAL NUMBER OF NORMALIZED MUTATIONS PER SPECIES
# create a dataset with total number of mutations and with all 12 types of substitutions
AggTotalMutSpectrumPerSpecies = aggregate(MUT$NumberOfSynMutPerSpecies, by = list(MUT$Species), FUN = sum); 
names(AggTotalMutSpectrumPerSpecies) = c('Species','TotalMutRate')
MutTypes = data.frame(VecOfNormalSubstitutions); names(MutTypes)=c('MutType')
nrow(AggTotalMutSpectrumPerSpecies)
AggTotalMutSpectrumPerSpecies = merge(AggTotalMutSpectrumPerSpecies,MutTypes)
nrow(AggTotalMutSpectrumPerSpecies) 

# count total number of all types of substitutions and merge with dataset above
AggTotalMutSpectrumPerSpeciesPerMutType = aggregate(MUT$NumberOfSynMutPerSpecies, by = list(MUT$Species,MUT$Subs), FUN = sum); 
names(AggTotalMutSpectrumPerSpeciesPerMutType) = c('Species','MutType','MutTypeRate')
nrow(AggTotalMutSpectrumPerSpeciesPerMutType) # 
ALL = merge(AggTotalMutSpectrumPerSpecies,AggTotalMutSpectrumPerSpeciesPerMutType, by = c('Species','MutType'), all.x=TRUE)
ALL[is.na(ALL)]<-0   # some classes of substitutions are absent from some animals => they are transformed to zeroes. 

#### SHOULD WE DELETE SPECIES WITH TOO MANY ZEROES?
ALL[ALL$MutTypeRate ==0,] %>% count(Species, sort = TRUE) -> cs; hist(cs$n, breaks = 50)
SpeciesWithTooManyZeroMutTypes = unique(cs[cs$n >= ParamNumberOfZeroes,]$Species); length(SpeciesWithTooManyZeroMutTypes) # 240 species to delete
ALL = ALL[!ALL$Species %in% SpeciesWithTooManyZeroMutTypes,]
ALL$Fraction = ALL$MutTypeRate/ALL$TotalMutRate; summary(ALL$Fraction)

###### create matrix for PCA:
AT = ALL[ALL$MutType == 'A_T',]; AT = AT[c(1,5)]; names(AT) = c('Species','AT'); AT = AT[order(AT$Species),]
AG = ALL[ALL$MutType == 'A_G',]; AG = AG[c(1,5)]; names(AG) = c('Species','AG'); AG = AG[order(AG$Species),]
AC = ALL[ALL$MutType == 'A_C',]; AC = AC[c(1,5)]; names(AC) = c('Species','AC'); AC = AC[order(AC$Species),]
TA = ALL[ALL$MutType == 'T_A',]; TA = TA[c(1,5)]; names(TA) = c('Species','TA'); TA = TA[order(TA$Species),]
TG = ALL[ALL$MutType == 'T_G',]; TG = TG[c(1,5)]; names(TG) = c('Species','TG'); TG = TG[order(TG$Species),]
TC = ALL[ALL$MutType == 'T_C',]; TC = TC[c(1,5)]; names(TC) = c('Species','TC'); TC = TC[order(TC$Species),]
CA = ALL[ALL$MutType == 'C_A',]; CA = CA[c(1,5)]; names(CA) = c('Species','CA'); CA = CA[order(CA$Species),]
CG = ALL[ALL$MutType == 'C_G',]; CG = CG[c(1,5)]; names(CG) = c('Species','CG'); CG = CG[order(CG$Species),]
CT = ALL[ALL$MutType == 'C_T',]; CT = CT[c(1,5)]; names(CT) = c('Species','CT'); CT = CT[order(CT$Species),]
GA = ALL[ALL$MutType == 'G_A',]; GA = GA[c(1,5)]; names(GA) = c('Species','GA'); GA = GA[order(GA$Species),]
GC = ALL[ALL$MutType == 'G_C',]; GC = GC[c(1,5)]; names(GC) = c('Species','GC'); GC = GC[order(GC$Species),]
GT = ALL[ALL$MutType == 'G_T',]; GT = GT[c(1,5)]; names(GT) = c('Species','GT'); GT = GT[order(GT$Species),]

###### PCA 
MATRIX = cbind(AT,AG[,2],AC[,2],TA[,2],TG[,2],TC[,2],CA[,2],CG[,2],CT[,2],GA[,2],GC[,2],GT[,2]); names(MATRIX) = c('Species','AT','AG','AC','TA','TG','TC','CA','CG','CT','GA','GC','GT')
row.names(MATRIX)=MATRIX$Species
head(MATRIX)
matrix = MATRIX[,c(2:13)]
PCA = prcomp(matrix, center = TRUE, scale = ParamScaleOrNot) #FALSE) # I don't scale because we analyze the same units (fraction from MutSpec) 
print(PCA)  
summary(PCA)
MATRIX$Pca1 = PCA$x[,1]
MATRIX$Pca2 = PCA$x[,2]
MATRIX$Pca3 = PCA$x[,3]
MATRIX$Pca4 = PCA$x[,4]
MATRIX$Pca5 = PCA$x[,5]
MATRIX$Pca6 = PCA$x[,6]
MATRIX$Pca7 = PCA$x[,7]
MATRIX$Pca8 = PCA$x[,8]
MATRIX$Pca9 = PCA$x[,9]
MATRIX$Pca10 = PCA$x[,10]
MATRIX$Pca11 = PCA$x[,11]
MATRIX$Pca12 = PCA$x[,12]
MATRIX = merge(MATRIX,Taxa, by = 'Species', all.x = TRUE) # some animals are missing
MATRIX = merge(MATRIX,GenerTime, by = 'Species', all.x = TRUE) # some animals are missing
MATRIX = merge(MATRIX,AnAge, by = 'Species', all.x = TRUE)
MATRIX = merge(MATRIX,Temp, by = 'Species', all.x = TRUE)
nrow(MATRIX[is.na(MATRIX$Species),])
plot(MATRIX$Pca1,MATRIX$Pca2)
plot(MATRIX$Pca2,MATRIX$Pca3)

table(MATRIX$Class)

###### FIRST COMPONENT (only Mammals and Actinopteri are two big taxa which we have to analyze - they show that GA positively correlate with BMR:
boxplot(MATRIX[MATRIX$Class == 'Actinopteri',]$Pca1,MATRIX[MATRIX$Class == 'Amphibia',]$Pca1,MATRIX[MATRIX$Class == 'Squamata' | MATRIX$Class == 'Testudines' |  MATRIX$Class == 'Crocodylia',]$Pca1,MATRIX[MATRIX$Class == 'Mammalia',]$Pca1,MATRIX[MATRIX$Class == 'Aves',]$Pca1, notch = TRUE, outline = FALSE, varwidth = TRUE)
boxplot(MATRIX[MATRIX$Class == 'Actinopteri',]$GA,MATRIX[MATRIX$Class == 'Amphibia',]$GA,MATRIX[MATRIX$Class == 'Squamata' | MATRIX$Class == 'Testudines' |  MATRIX$Class == 'Crocodylia',]$GA,MATRIX[MATRIX$Class == 'Mammalia',]$GA,MATRIX[MATRIX$Class == 'Aves',]$GA, notch = TRUE, outline = FALSE, varwidth = TRUE)

# G>A in fishes positively correlate with BMR (probably BMR represent body temperature)
MATRIX[MATRIX$Class == 'Actinopteri',]$Species
MATRIX[MATRIX$Species == 'Thunnus_thynnus',]

nrow(MATRIX[MATRIX$Class == 'Actinopteri',]) # 470
setwd('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED')
FISH = MATRIX[!is.na(MATRIX$Class) & MATRIX$Class == 'Actinopteri',]
FISH = FISH[order(FISH$GA),]
write.table(FISH, file = 'FISHES.txt')

nrow(MATRIX[MATRIX$Class == 'Actinopteri' & !is.na(MATRIX$Basal.metabolic.rate..watts.),])
MATRIX[MATRIX$Class == 'Actinopteri' & !is.na(MATRIX$Basal.metabolic.rate..watts.),]$Basal.metabolic.rate..watts.
MATRIX[MATRIX$Class == 'Actinopteri' & !is.na(MATRIX$Basal.metabolic.rate..watts.),]$Species
MATRIX[MATRIX$Class == 'Actinopteri' & !is.na(MATRIX$Basal.metabolic.rate..watts.),]$GA

cor.test(MATRIX[MATRIX$Class == 'Actinopteri',]$GA,MATRIX[MATRIX$Class == 'Actinopteri',]$Basal.metabolic.rate..watts., method = 'spearman') # (if min = 4 subst (quantile = 0.05); rho = 0.33; p = 0.023, N = 49)
plot(MATRIX[MATRIX$Class == 'Actinopteri',]$GA,MATRIX[MATRIX$Class == 'Actinopteri',]$Basal.metabolic.rate..watts.)
cor.test(MATRIX[MATRIX$Class == 'Actinopteri',]$GA,MATRIX[MATRIX$Class == 'Actinopteri',]$T..oC., method = 'spearman') # a bit negative with temperature

# G>A in mammals negatively correlate with BMR (probably BMR negatively correlate with temperature in mammals)
cor.test(MATRIX[MATRIX$Class == 'Mammalia',]$Adult.weight..g.,(MATRIX[MATRIX$Class == 'Mammalia',]$Metabolic.rate..W.), method = 'spearman') # super positive
cor.test(MATRIX[MATRIX$Class == 'Mammalia',]$Adult.weight..g.,MATRIX[MATRIX$Class == 'Mammalia',]$Temperature..K., method = 'spearman') # positive
cor.test(MATRIX[MATRIX$Class == 'Mammalia',]$Metabolic.rate..W.,MATRIX[MATRIX$Class == 'Mammalia',]$Temperature..K., method = 'spearman') # positive
cor.test(MATRIX[MATRIX$Class == 'Mammalia',]$GA,(MATRIX[MATRIX$Class == 'Mammalia',]$Metabolic.rate..W.), method = 'spearman') # nothing
cor.test(MATRIX[MATRIX$Class == 'Mammalia',]$GA,MATRIX[MATRIX$Class == 'Mammalia',]$Adult.weight..g., method = 'spearman') # nothing
cor.test(MATRIX[MATRIX$Class == 'Mammalia',]$GA,MATRIX[MATRIX$Class == 'Mammalia',]$Temperature..K., method = 'spearman') # nothing
cor.test(MATRIX[MATRIX$Class == 'Mammalia',]$GA,(MATRIX[MATRIX$Class == 'Mammalia',]$Metabolic.rate..W./MATRIX[MATRIX$Class == 'Mammalia',]$Adult.weight..g.), method = 'spearman') # positive
names(MATRIX)

plot(MATRIX[MATRIX$Class == 'Mammalia',]$GA,log2(MATRIX[MATRIX$Class == 'Mammalia',]$Metabolic.rate..W./MATRIX[MATRIX$Class == 'Mammalia',]$Adult.weight..g.)) # positive
cor.test(MATRIX[MATRIX$Class == 'Mammalia',]$Pca1,(MATRIX[MATRIX$Class == 'Mammalia',]$Metabolic.rate..W./MATRIX[MATRIX$Class == 'Mammalia',]$Adult.weight..g.), method = 'spearman') # positive
plot(MATRIX[MATRIX$Class == 'Mammalia',]$Pca1,log2(MATRIX[MATRIX$Class == 'Mammalia',]$Metabolic.rate..W./MATRIX[MATRIX$Class == 'Mammalia',]$Adult.weight..g.)) # positive

###### SECOND COMPONENT: 
cor.test(MATRIX$Pca2,MATRIX$GenerationLength_d)
plot(MATRIX$Pca2,log2(MATRIX$GenerationLength_d))
cor.test(MATRIX$TC,MATRIX$GenerationLength_d)

cor.test(MATRIX[MATRIX$Class == 'Actinopteri',]$TC,MATRIX[MATRIX$Class == 'Actinopteri',]$Female.maturity..days., method = 'spearman')
cor.test(MATRIX[MATRIX$Class == 'Actinopteri',]$Pca2,MATRIX[MATRIX$Class == 'Actinopteri',]$Female.maturity..days., method = 'spearman')
plot(MATRIX[MATRIX$Class == 'Actinopteri',]$Pca2,MATRIX[MATRIX$Class == 'Actinopteri',]$Female.maturity..days.)



  









####################
#### FIND ECOLOGICAL INTERPRETATIONS OF PCs
####################

###### EXTREMES
MATRIX[MATRIX$Pca1 < quantile(MATRIX$Pca1,0.05),]$Species 
MATRIX[MATRIX$Pca1 > quantile(MATRIX$Pca1,0.95),]$Species 
MATRIX[MATRIX$Pca2 < quantile(MATRIX$Pca2,0.05),]$Species 
MATRIX[MATRIX$Pca2 > quantile(MATRIX$Pca2,0.95),]$Species 
MATRIX[MATRIX$Pca3 < quantile(MATRIX$Pca3,0.05),]$Species 
MATRIX[MATRIX$Pca3 > quantile(MATRIX$Pca3,0.95),]$Species   

###### GENERATION LENGTH
GenerTime = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/GenerationLenghtforAllMammals/GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
GenerTime$Species = gsub(' ','_',GenerTime$Scientific_name)
GenerTime = GenerTime[,c(11,13)]

Test = merge(MATRIX,GenerTime) 
cor.test(Test$Pca1,Test$GenerationLength_d, method = 'spearman') # no
cor.test(Test$Pca2,Test$GenerationLength_d, method = 'spearman') # super negative (-AG) -0.34
cor.test(Test$Pca3,Test$GenerationLength_d, method = 'spearman') # positive (+TC)       +0.26
cor.test(Test$Pca4,Test$GenerationLength_d, method = 'spearman') # no
cor.test(Test$Pca5,Test$GenerationLength_d, method = 'spearman') # no

####### AnAge 
AnAge = read.table('1_RAW/anage_data.txt', header = TRUE, sep = '\t')
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') 
AnAge = AnAge[AnAge$Class == 'Mammalia',];  
names(AnAge)

### Adult.weight..g.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Adult.weight..g."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca3,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Maximum.longevity..yrs.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Maximum.longevity..yrs."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca3,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # a bit positive (but unlikely survive multiple test correction)

### Female.maturity..days.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Female.maturity..days."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca3,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Weaning..days.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Weaning..days."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca3,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Birth.weight..g.  
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Birth.weight..g."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca3,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # a bit positive (unlikely will pass multiple test correction)

### Litters.Clutches.per.year
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Litters.Clutches.per.year"))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca3,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Litter.Clutch.size
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Litter.Clutch.size"))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca3,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Weaning.weight..g.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Weaning.weight..g."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca3,Test[,19], method = 'spearman') # no
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Metabolic.rate..W. 
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Metabolic.rate..W."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)  
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # no
cor.test(Test$Pca3,Test[,19], method = 'spearman') # no
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Temperature..K.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Temperature..K."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # no
cor.test(Test$Pca3,Test[,19], method = 'spearman') # no
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Growth.rate..1.days.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Growth.rate..1.days."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # no
cor.test(Test$Pca3,Test[,19], method = 'spearman') # no
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

###### HETEROTERMS
Hib = read.table('1_RAW/hibernation/tornew5.csv', header = TRUE, sep = ',')
Hib$Species = gsub(" ",'_',Hib$Species)
HibOnlySpecies = Hib[Hib$Type == 'HIB',]$Species; length(HibOnlySpecies)
DtOnlySpecies = Hib[Hib$Type == 'DT',]$Species; length(DtOnlySpecies)
HibAndDtSpecies = Hib[Hib$Type == 'HIB' || Hib$Type == 'DT',]$Species; length(HibAndDtSpecies)
boxplot(MATRIX[MATRIX$Species %in% HibOnlySpecies,]$Pca1,MATRIX[!MATRIX$Species %in% HibOnlySpecies,]$Pca1, notch = TRUE, names = c('HibOnlySpecies','Other'), ylab = 'PC1'); # !!!! HibSpecies have a bit lower Pc1
boxplot(MATRIX[MATRIX$Species %in% DtOnlySpecies,]$Pca1,MATRIX[!MATRIX$Species %in% DtOnlySpecies,]$Pca1, notch = TRUE, names = c('HibOnlySpecies','Other'), ylab = 'PC1');   # nothing
# control for body mass? they are small! # other tests - decreased BMR as compared to body mass - compare strong outliers (cold and hot)...; animals, which live more than should according to bosy mass (naked mole rat)

###### MARSUPIALS, BATS... AGAIN NEED NORMAL TAXONOMY !!!!!!!!
Marsupials = c(AnAge[AnAge$Order == 'Diprotodontia' | AnAge$Order == 'Didelphimorphia' | AnAge$Order == 'Dasyuromorphia',]$Species)
Placental = c(AnAge[AnAge$Order == 'Artiodactyla' | AnAge$Order == 'Cetacea' | AnAge$Order == 'Carnivora',]$Species)
boxplot(MATRIX[MATRIX$Species %in% Marsupials,]$Pca1,MATRIX[MATRIX$Species %in% Placental,]$Pca1, notch = TRUE); # !!!! HibSpecies have a bit lower Pc1

###### FIGURES:

MATRIX = merge(MATRIX,GenerTime)
MATRIX = MATRIX[order(MATRIX$GenerationLength_d),]
MATRIX$Col = c(rep('green',150),rep('gray',187),rep('red',150))
summary(PCA)
print(PCA)

pdf('4_FIGURES/PCA.pdf', width = 14, height = 14)
par(mfcol=c(2,3))
summary(PCA)
#plot(PCA)
plot(MATRIX$Pca1,MATRIX$Pca2), col = MATRIX$Col)
plot(MATRIX$Pca2,MATRIX$Pca3, col = MATRIX$Col)
# plot(PCA$x[,1],MATRIX$GenerationLength_d); cor.test(PCA$x[,1],MATRIX$GenerationLength_d, method = 'spearman') # nothing  - First mutagen signature! Body mass normalized BMR!
plot(MATRIX$Pca2,log2(MATRIX$GenerationLength_d)); cor.test(MATRIX$Pca2,MATRIX$GenerationLength_d, method = 'spearman')
plot(MATRIX$Pca3,log2(MATRIX$GenerationLength_d)); cor.test(MATRIX$Pca2,MATRIX$GenerationLength_d, method = 'spearman') 
biplot(PCA, col = c('grey','black'), cex = 0.5)
biplot(PCA, choices=c(2,3), col = c('grey','black'), cex = 0.5) #  biplot(princomp(USArrests),choices=c(1,3))
dev.off()

####################
#### FIND DNA POLYMERAZE SIGNATURE (last PCs)!!!!!
####################
# LOGIC IS THE NEXT: first - third components are driven by temperature and ecology. 
# other components are not driven by ecology, physiology... sow we need to subtract effect of the first three PC and get naked signature of DNA polymeraze!!!
# how to do it carefully?!

pdf('4_FIGURES/DnaPolymerazeSignature.pdf', width = 14, height = 14)
biplot(PCA, choices=c(4,5), col = c('grey','black'), cex = 0.5) #  biplot(princomp(USArrests),choices=c(1,3))
dev.off()

PCA$x # PC's
PCA$sdev # the eigenvalues (res$sdev) giving information on the magnitude of each PC, 
PCA$rotation # and the loadings (res$rotation).

########### exercise with PCA and going back to data with all PCs - REMOVE something
par(mfrow=c(2,2))
PcaScale = prcomp(matrix, center = TRUE, scale = TRUE) 
PcaScale$scale  # sd of each column from the original matrix: sd(matrix[,1]) 
PcaScale$center # mean of each column from the original matrix: mean(matrix[,1])
head(PcaScale$x)
start = 1; end = 12
BACK <- PcaScale$x[,start:end] %*% t(PcaScale$rotation[,start:end]) # reconstruct everything taking into account some PCs
BACK <- scale(BACK, center = FALSE , scale=1/PcaScale$scale) # divide by sd 
BACK <- scale(BACK, center = -1 * PcaScale$center, scale=FALSE) # 
gamma = data.frame(BACK); 
barplot(c(mean(gamma$AT),mean(gamma$AG),mean(gamma$AC),mean(gamma$TA),mean(gamma$TG),mean(gamma$TC),mean(gamma$CA),mean(gamma$CG),mean(gamma$CT),mean(gamma$GA),mean(gamma$GC),mean(gamma$GT))  , names = c('AT','AG','AC','TA','TG','TC','CA','CG','CT','GA','GC','GT'), main = 'EXO')

PcaScale = prcomp(matrix, center = TRUE, scale = TRUE) 
PcaScale$scale  # sd of each column from the original matrix: sd(matrix[,1]) 
PcaScale$center # mean of each column from the original matrix: mean(matrix[,1])
start = 1; end = 1
BACK <- PcaScale$x[,start:end] %*% t(PcaScale$rotation[,start:end])
BACK <- scale(BACK, center = FALSE , scale=1/PcaScale$scale)
BACK <- scale(BACK, center = -1 * PcaScale$center, scale=FALSE)
gamma = data.frame(BACK)
barplot(c(mean(gamma$AT),mean(gamma$AG),mean(gamma$AC),mean(gamma$TA),mean(gamma$TG),mean(gamma$TC),mean(gamma$CA),mean(gamma$CG),mean(gamma$CT),mean(gamma$GA),mean(gamma$GC),mean(gamma$GT))  , names = c('AT','AG','AC','TA','TG','TC','CA','CG','CT','GA','GC','GT'), main = 'EXO')

start = 4; end = 12
BACK <- PcaScale$x[,start:end] %*% t(PcaScale$rotation[,start:end])
BACK <- scale(BACK, center = FALSE , scale=1/PcaScale$scale)
BACK <- scale(BACK, center = -1 * PcaScale$center, scale=FALSE)
gamma = data.frame(BACK)
barplot(c(mean(gamma$AT),mean(gamma$AG),mean(gamma$AC),mean(gamma$TA),mean(gamma$TG),mean(gamma$TC),mean(gamma$CA),mean(gamma$CG),mean(gamma$CT),mean(gamma$GA),mean(gamma$GC),mean(gamma$GT))  , names = c('AT','AG','AC','TA','TG','TC','CA','CG','CT','GA','GC','GT'), main = 'EXO')
