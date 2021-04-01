#### read AllGenesCodonUsageNoOverlap.txt

rm(list=ls(all=TRUE))

MutSpec = read.table("../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt", header = TRUE) 

### read Bird phenotypes from AnAge
anage = read.table('../../Body/1Raw/anage_data.txt', header = TRUE, sep = '\t')
anage$Species = paste(anage$Genus, anage$Species, sep = '_')
table(anage$Class)
anage = anage[anage$Class == 'Aves',]

### merge MutSpec with PHE by 'Species'
MutSpecPhe = merge(MutSpec,anage, by = 'Species')

### descriptive analyses: go column by column and if there are many notNA run cor test:
names(MutSpecPhe)
for (i in 22:42)
{ # i = 40
  if (length(MutSpecPhe[!is.na(MutSpecPhe[,i]),][,i]) > 10)
  {
    print(names(MutSpecPhe)[i])
    print(length(MutSpecPhe[!is.na(MutSpecPhe[,i]),][,i]))
    print(cor.test(MutSpecPhe$T_C,MutSpecPhe[,i], method = 'spearman'))
  }
}

### the higher the clutch size, the higher the Ah>Gh. little relative eggs? R strategy? 
temp = MutSpecPhe[!is.na(MutSpecPhe$Litter.Clutch.size),]
table(temp$Litter.Clutch.size)
summary(temp$Litter.Clutch.size)
boxplot(temp[temp$Litter.Clutch.size<=4,]$T_C,temp[temp$Litter.Clutch.size>4,]$T_C, notch = TRUE)


