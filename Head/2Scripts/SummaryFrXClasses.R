rm(list=ls(all=TRUE))


SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')

#take only 12 genes
SynNuc = SynNuc[SynNuc$Gene !='ND6', ]

# count sum of all nucleotides btw Species and Classes
AGG = aggregate(list(af$NeutralA,af$NeutralT,af$NeutralG,af$NeutralC), by = list(af$Species,af$Class, af$Taxonomy), FUN = sum)
names(AGG) = c('Species','Class','Taxonomy','NeutralA','NeutralT','NeutralG','NeutralC')

## how many classes we have
unique(AGG$Class)

# count fractions of each nucleotide, annot LIGHT CHAIN
AGG$FrA = AGG$NeutralA / (AGG$NeutralA + AGG$NeutralT + AGG$NeutralG + AGG$NeutralC)
AGG$FrT = AGG$NeutralT / (AGG$NeutralA + AGG$NeutralT + AGG$NeutralG + AGG$NeutralC) 
AGG$FrG = AGG$NeutralG / (AGG$NeutralA + AGG$NeutralT + AGG$NeutralG + AGG$NeutralC) 
AGG$FrC = AGG$NeutralC / (AGG$NeutralA + AGG$NeutralT + AGG$NeutralG + AGG$NeutralC)

### summary for FrA btw classes

summary(AGG[AGG$Class == 'Actinopterygii', ]$FrA)
summary(AGG[AGG$Class == 'Amphibia',]$FrA)
summary(AGG[AGG$Class == 'Aves', ]$FrA)
summary(AGG[AGG$Class == 'Reptilia',]$FrA)
summary(AGG[AGG$Class == 'Mammalia',]$FrA)
summary(AGG[AGG$Class == 'AncientFish',]$FrA)


### summary for FrT btw classes
summary(AGG[AGG$Class == 'Actinopterygii', ]$FrT)
summary(AGG[AGG$Class == 'Amphibia',]$FrT)
summary(AGG[AGG$Class == 'Aves', ]$FrT)
summary(AGG[AGG$Class == 'Reptilia',]$FrT)
summary(AGG[AGG$Class == 'Mammalia',]$FrT)
summary(AGG[AGG$Class == 'AncientFish',]$FrT)


### summary for FrC btw classes
summary(AGG[AGG$Class == 'Actinopterygii', ]$FrC)
summary(AGG[AGG$Class == 'Amphibia',]$FrC)
summary(AGG[AGG$Class == 'Aves', ]$FrC)
summary(AGG[AGG$Class == 'Reptilia',]$FrC)
summary(AGG[AGG$Class == 'Mammalia',]$FrC)
summary(AGG[AGG$Class == 'AncientFish',]$FrC)

### summary for FrG btw classes
summary(AGG[AGG$Class == 'Actinopterygii', ]$FrG)
summary(AGG[AGG$Class == 'Amphibia',]$FrG)
summary(AGG[AGG$Class == 'Aves', ]$FrG)
summary(AGG[AGG$Class == 'Reptilia',]$FrG)
summary(AGG[AGG$Class == 'Mammalia',]$FrG)
summary(AGG[AGG$Class == 'AncientFish',]$FrC)



