#### read AllGenesCodonUsageNoOverlap.txt

rm(list=ls(all=TRUE))

### unzip, keep in R memory and delete from folder
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip") 
ALL = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t') 
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")} 

### delete ND6
nrow(ALL)
ALL = ALL[ALL$Gene != 'ND6',]
nrow(ALL)

#### aggregate (summ up) neutral nucleotides for each species
SynNuc = aggregate(list(ALL$NeutralA,ALL$NeutralT,ALL$NeutralG,ALL$NeutralC), by = list(ALL$Species,ALL$Class,ALL$Taxonomy), FUN = sum)
names(SynNuc)=c('Species','Class','Taxonomy','NeutralA','NeutralT','NeutralG','NeutralC')
table(SynNuc$Class)

#### estimate fraction of each nucleotide 
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 

### read Bird phenotypes from Valya
PHE = read.table("../../Body/1Raw/DataFromValya/ALL_PHENOTYPES.txt", header = FALSE, sep = ',') 
names(PHE) = c('Species','PhenotypeCode')

### merge SynNuc with PHE by 'Species'
SynNucPhe = merge(SynNuc,PHE, by = 'Species', all.y = TRUE)

VectorOfSpeceisWithourRefSeqData = SynNucPhe[is.na(SynNucPhe$Class),]$Species
length(VectorOfSpeceisWithourRefSeqData)

SynNucPhe = SynNucPhe[!is.na(SynNucPhe$Class),]
table(SynNucPhe$Class)
table(SynNucPhe$PhenotypeCode)
#HI - high altitude
#DI - divers
#FM - long-distance migration
#NF - flightless
#WI - wintering
#AI - outstanding flight abilities
#TROPIC - tropic birds for reference

### first descriptive analyses

# a priory hyp: (low metaboilism => low Ah>Gh)
boxplot(SynNucPhe$FrT~SynNucPhe$PhenotypeCode, notch = TRUE) # !!!! flightless has excess of Ah (Tl)
boxplot(SynNucPhe$FrC~SynNucPhe$PhenotypeCode, notch = TRUE) # !!!! flightless has deficit of Gh (Cl)

boxplot(SynNucPhe$FrA~SynNucPhe$PhenotypeCode, notch = TRUE)
boxplot(SynNucPhe$FrG~SynNucPhe$PhenotypeCode, notch = TRUE) 

### nice plot with all phenotype groups
str(SynNucPhe)
plot(SynNucPhe$FrT,SynNucPhe$FrC, col = labels(as.factor(SynNucPhe$PhenotypeCode)))
legend("topleft", col = labels(as.factor(SynNucPhe$PhenotypeCode)), legend = levels(as.factor(SynNucPhe$PhenotypeCode)))

### TO THINK 
## AnAge and check if there is correlattion of MetRate (temp) with different phenotypes?
## MutSpecFromPolymorphisms (the same) !!!!




