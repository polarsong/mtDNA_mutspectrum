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

### read Bird phenotypes from AnAge
anage = read.table('../../Body/1Raw/anage_data.txt', header = TRUE, sep = '\t')
anage$Species = paste(anage$Genus, anage$Species, sep = '_')
table(anage$Class)
anage = anage[anage$Class == 'Aves',]

### merge MutSpec with PHE by 'Species'
SynNucPhe = merge(SynNuc,anage, by = 'Species')

### descriptive analyses: go column by column and if there are many notNA run cor test:
names(SynNucPhe)
for (i in 20:40)
{ # i = 40
  if (length(SynNucPhe[!is.na(SynNucPhe[,i]),][,i]) > 10)
  {
    print(names(SynNucPhe)[i])
    print(length(SynNucPhe[!is.na(SynNucPhe[,i]),][,i]))
    print(cor.test(SynNucPhe$FrT,SynNucPhe[,i], method = 'spearman'))
  }
}

### main results: if Ah>Gh is higher in species with high clutch size, Ah can be:
### 1) less in species with high clutch size (r strategy?)
### 2) less in species with fast growth rate (the more th clutch size the faster should be growth rate ~ r strategy?)
### 3) higher in species with long maturation period (K strategy?)

"Growth.rate..1.days.
[1] 76

	Spearman's rank correlation rho

data:  SynNucPhe$FrT and SynNucPhe[, i]
S = 92314, p-value = 0.02225
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
-0.2619764 

[1] Litter.Clutch.size
[1] 72

	Spearman's rank correlation rho

data:  SynNucPhe$FrT and SynNucPhe[, i]
S = 76484, p-value = 0.05223
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
-0.2297272 

[1] Female.maturity..days.
[1] 111

	Spearman's rank correlation rho

data:  SynNucPhe$FrT and SynNucPhe[, i]
S = 178156, p-value = 0.02133
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho 
0.2183384 "

