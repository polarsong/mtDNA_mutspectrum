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

##### TUNE 1Raw/Aves and clutches tuned.txt AND SAVE IT TO 2Derived/AvesClutchSizes.txt
Clutch = read.table('../../Body/1Raw/Aves and clutches tuned.txt', header = FALSE, sep = '\n')
for (i in 1:100) {
Clutch$V1 = gsub('  ','\t',Clutch$V1);
Clutch$V1 = gsub(' \t','\t',Clutch$V1);
Clutch$V1 = gsub('\t ','\t',Clutch$V1);
Clutch$V1 = gsub('\t\t','\t',Clutch$V1); 
}
Clutch$V1 = gsub(' ','_',Clutch$V1); 
write.table(Clutch,"../../Body/2Derived/AvesClutchSizes.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
Clutch = read.table("../../Body/2Derived/AvesClutchSizes.txt", header = TRUE, sep = "\t")

### merge MutSpec with PHE by 'Species'
MutSpecPhe = merge(MutSpec,Clutch, by = 'Species')
nrow(MutSpecPhe) # 152

### descriptive analyses
summary(MutSpecPhe$Clutch)
cor.test(MutSpecPhe$T_C,MutSpecPhe$Clutch, method = 'spearman')
