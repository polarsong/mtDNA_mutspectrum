rm(list=ls(all=TRUE))

if (!require(caper)) install.packages("caper")
if (!require(geiger)) install.packages("geiger")
if (!require(ggpubr)) install.packages("ggpubr")
library(caper)
library(geiger)
library("ggpubr")
###########Taxonomy###################################################################
Taxa = read.table("../../Body/1Raw/TaxaFromKostya.Names.stat", sep = '\t',header = FALSE) 
Taxa$Species = gsub(";.*",'',Taxa$V1); 
for (i in (1:nrow(Taxa)))  {Taxa$Species[i] = paste(unlist(strsplit(as.character(Taxa$Species[i]),split = ' '))[1],unlist(strsplit(as.character(Taxa$Species[i]),split = ' '))[2], sep = '_')}
Taxa$Class = gsub(";Chordata;.*",'',Taxa$V1); Taxa$Class = gsub(".*;",'',Taxa$Class); table(Taxa$Class)
Taxa$Class = gsub('Actinopteri','Actinopterygii',Taxa$Class)
Taxa$Class = gsub("Testudines|Squamata|Crocodylia|Sphenodontia",'Reptilia',Taxa$Class)
length(unique(Taxa$Species)) # 1708
table(Taxa$Class)
Taxa = Taxa[,-1]

TaxaMore = read.table("../../Body/1Raw/TaxaFromKostya.2NeedTaxa.tax.txt", sep = '\t',header = FALSE) 
TaxaMore$Species = ''
for (i in (1:nrow(TaxaMore)))  
{TaxaMore$Species[i] = paste(unlist(strsplit(as.character(TaxaMore$V1[i]),split = ' '))[1],unlist(strsplit(as.character(TaxaMore$V1[i]),split = ' '))[2], sep = '_')}
TaxaMore$Class = gsub("; Chordata;.*",'',TaxaMore$V2); 
TaxaMore$Class = gsub(".*; ",'',TaxaMore$Class); 
TaxaMore$Class = gsub('Actinopteri','Actinopterygii',TaxaMore$Class)
TaxaMore$Class = gsub("Testudines|Squamata|Crocodylia",'Reptilia',TaxaMore$Class)
table(TaxaMore$Class)
TaxaMore = TaxaMore[,-c(1,2)]

Taxa = rbind(Taxa,TaxaMore); Taxa = unique(Taxa)
#########################################################################################

MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)
MUTFROMK = merge(MUT, Taxa)
table(MUTFROMK$Class)
MUTFROMK = merge(MUT, Taxa)
MUTFROMK = MUTFROMK[MUTFROMK$T_C > 0,]
MUTFROMK = MUTFROMK[MUTFROMK$A_G > 0,]
MUTFROMK$TCdivAG = MUTFROMK$T_C / MUTFROMK$A_G
######reading ecology table from us + Kuptsov 
kuptsovtable = read.table("../../Body/2Derived/EcologyMammalianTable01_KuptsovA_ver2_Full.txt", sep='\t', header=TRUE)

MUTALL = merge(MUTFROMK, kuptsovtable)
MUTALL$Temper = as.numeric(gsub(",", ".", MUTALL$Temperature.C._White2003.2006.other.close.species))
MUTALL$GenerationLength_d = as.numeric(gsub(",", ".", MUTALL$GenerationLength_d))
summary(lm(formula = TCdivAG ~ Temper + GenerationLength_d, data = MUTALL))
summary(lm(formula = TCdivAG ~ scale(Temper) + scale(GenerationLength_d), data = MUTALL))


MUTALL$MarsMono = MUTALL$Mars + MUTALL$Mono
table(MUTALL$MarsMono)
formediantemperature = MUTALL[!is.na(MUTALL$Temper),]$Temper
coldspeciesnames = MUTALL[MUTALL$Temper <= mean(formediantemperature) & !is.na(MUTALL$Temper),]$Species
MUTALL$colddummy = 0
MUTALL[MUTALL$Species %in% coldspeciesnames,]$colddummy = 1
MUTALL$allcolddummy = MUTALL$Hib.unconfirmedHib + MUTALL$Daily.unconfirmedDaily + MUTALL$MarsMono + MUTALL$colddummy
table(MUTALL$allcolddummy)
MUTALL[MUTALL$allcolddummy > 0,]$allcolddummy = 1

summary(lm(formula = TCdivAG ~ allcolddummy + GenerationLength_d, data = MUTALL))
summary(lm(formula = TCdivAG ~ scale(GenerationLength_d), data = MUTALL))

unique(MUTALL$Order..or.infraorder.for.Cetacea.) 
MUTALL = MUTALL[MUTALL$Class != "Chiroptera",]
MUTALL = MUTALL[MUTALL$MarsMono != 1,]

tree = read.tree('../../Body/1Raw/mtalign.aln.treefile.rooted')

row.names(MUTALL) = MUTALL$Species

tree_pruned = treedata(tree, MUTALL, sort=T, warnings=T)$phy 

data<-as.data.frame(treedata(tree_pruned, MUTALL, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)

data$TCdivAG = as.numeric(as.character(data$TCdivAG))
data$Temper = as.numeric(as.character(data$Temper))
data$GenerationLength_d = as.numeric(as.character(data$GenerationLength_d))
data$allcolddummy = as.numeric(as.character(data$allcolddummy))


data_comp <- comparative.data(tree_pruned, data[, c('Species', 'TCdivAG',
                                                    'GenerationLength_d', 'Temper', 'allcolddummy')], Species, vcv=TRUE)

model = pgls(TCdivAG ~ scale(Temper) + scale(GenerationLength_d), data_comp, lambda="ML")
summary(model)

model = pgls(TCdivAG ~ allcolddummy + scale(GenerationLength_d), data_comp, lambda="ML")
summary(model)


model = pgls(TCdivAG ~ scale(GenerationLength_d), data_comp, lambda="ML")
summary(model)
