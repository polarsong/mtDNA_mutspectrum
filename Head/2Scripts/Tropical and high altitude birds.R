rm(list=ls(all=TRUE))
library('ggplot2')
#merging some phenotypes and life history traits to find some good correlations
#reading phenotypes
pheno = read.table('../../Body/1Raw/DataFromValya/ALL_PHENOTYPES.txt')
#reading anage
anage = read.table('../../Body/1Raw/anage_data.txt',header = TRUE, sep = '\t')
names(pheno) = c('Species', 'Phenotype')
#getting only tropical birds
only_tropical = pheno[pheno$Phenotype == ",\"TROPIC\"",]
#getting only high_altitude birds
only_altitude = pheno[pheno$Phenotype == ",\"HI\"",]
#reading codon usage
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip") 
codus = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t') 
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")} codus = codus[codus$Gene != 'ND6',]
codus = codus[codus$Gene != 'ND6',]
#### aggregate (summ up) neutral nucleotides for each species
SynNuc = aggregate(list(codus$NeutralA,codus$NeutralT,codus$NeutralG,codus$NeutralC), by = list(codus$Species,codus$Class,codus$Taxonomy), FUN = sum)
names(SynNuc)=c('Species','Class','Taxonomy','NeutralA','NeutralT','NeutralG','NeutralC')
table(SynNuc$Class)

#### estimate fraction of each nucleotide 
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
#merging high altitude birds and tropical birds with SynNuc
tropical_and_fr = merge(SynNuc, only_tropical, by = 'Species')
altitude_and_fr = merge(SynNuc, only_altitude, by = 'Species')
random_rows = sample(1:81, 18, replace = FALSE)
random_birds_from_tropical = data.frame()
for (i in random_rows)
{
  random_birds_from_tropical = rbind(random_birds_from_tropical, tropical_and_fr[i,])
}
#graphics
boxplot(altitude_and_fr$FrA, random_birds_from_tropical$FrA)
boxplot(altitude_and_fr$FrG, random_birds_from_tropical$FrG)
boxplot(altitude_and_fr$FrC, random_birds_from_tropical$FrC)
boxplot(altitude_and_fr$FrT, random_birds_from_tropical$FrT)
