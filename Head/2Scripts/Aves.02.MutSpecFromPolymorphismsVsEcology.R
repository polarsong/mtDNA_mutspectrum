#### read AllGenesCodonUsageNoOverlap.txt

rm(list=ls(all=TRUE))

MutSpec = read.table("../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt", header = TRUE) 

### read Bird phenotypes from Valya
PHE = read.table("../../Body/1Raw/DataFromValya/ALL_PHENOTYPES.txt", header = FALSE, sep = ',') 
names(PHE) = c('Species','PhenotypeCode')

### merge MutSpec with PHE by 'Species'
MutSpecPhe = merge(MutSpec,PHE, by = 'Species')
table(MutSpecPhe$PhenotypeCode)
#HI - high altitude
#DI - divers
#FM - long-distance migration
#NF - flightless
#WI - wintering
#AI - outstanding flight abilities
#TROPIC - tropic birds for reference

### first descriptive analyses

# a priory hyp: (low metabolism => low Ah>Gh)
boxplot(MutSpecPhe$T_C~MutSpecPhe$PhenotypeCode, notch = TRUE) # !!!! flightless has excess of Ah (Tl)


### nice plot with all phenotype groups
str(MutSpecPhe)
#plot(MutSpecPhe$FrT,MutSpecPhe$FrC, col = labels(as.factor(MutSpecPhe$PhenotypeCode)))
#legend("topleft", col = labels(as.factor(SynNucPhe$PhenotypeCode)), legend = levels(as.factor(SynNucPhe$PhenotypeCode)))

### TO THINK:
## correlate AnAge with MutSpec
## AnAge and check if there is correlattion of MetRate (temp) with different phenotypes?
## MutSpecFromPolymorphisms (the same) !!!!

### rank and think
MutSpecPhe = MutSpecPhe[order(MutSpecPhe$T_C),]

## without phenotypes: correlate NucContent (01....R) with Substitutions (02...R) (N = 41) 





