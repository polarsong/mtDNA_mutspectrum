#### read AllGenesCodonUsageNoOverlap.txt

rm(list=ls(all=TRUE))

MutSpec = read.table("../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt", header = TRUE)

### add taxonomy
Taxa = read.table("../../Body/1Raw/TaxaWithClasses.txt")

## merge MutSpec with Taxa and keep only Aves
nrow(MutSpec)
MutSpec = merge(MutSpec,Taxa, by = 'Species')
nrow(MutSpec) # lose a bit...
table(MutSpec$Class)
MutSpec = MutSpec[MutSpec$Class == 'Aves',]
nrow(MutSpec)

### order and save
MutSpec = MutSpec[order(MutSpec$T_C),]
write.table(MutSpec, file = '../../Body/3Results/Aves.02B.MutSpecFromPolymorphismsVsAnAge.AllAvesForMeditation.csv', sep = ',', quote = FALSE, row.names = FALSE)

### read Bird phenotypes from AnAge
anage = read.table('../../Body/1Raw/anage_data.txt', header = TRUE, sep = '\t')
anage$Species = paste(anage$Genus, anage$Species, sep = '_')
table(anage$Class)
anage = anage[anage$Class == 'Aves',]

### merge MutSpec with PHE by 'Species'
MutSpecPhe = merge(MutSpec,anage, by = 'Species')

### descriptive analyses: go column by column and if there are many notNA run cor test:
names(MutSpecPhe)
for (i in c(23:34,39:43))
{ # i = 30
  if (length(MutSpecPhe[!is.na(MutSpecPhe[,i]),][,i]) > 10)
  { 
    if (as.numeric(cor.test(MutSpecPhe$T_C,MutSpecPhe[,i], method = 'spearman')[3]) < 0.1)
    {
    print(names(MutSpecPhe)[i])
    print(length(MutSpecPhe[!is.na(MutSpecPhe[,i]),][,i]))
    print(cor.test(MutSpecPhe$T_C,MutSpecPhe[,i], method = 'spearman'))
    }
  }
}

### the higher the clutch size, the higher the Ah>Gh. little relative eggs? R strategy? 
temp = MutSpecPhe[!is.na(MutSpecPhe$Litter.Clutch.size),]
table(temp$Litter.Clutch.size)
summary(temp$Litter.Clutch.size)
boxplot(temp[temp$Litter.Clutch.size<=4,]$T_C,temp[temp$Litter.Clutch.size>4,]$T_C, notch = TRUE)


