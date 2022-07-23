rm(list=ls(all=TRUE))

SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')

#take only 12 genes
SynNuc = SynNuc[SynNuc$Gene !='ND6', ]

# count sum of all nucleotides btw Species and Classes
AGG = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species,SynNuc$Class, SynNuc$Taxonomy), FUN = sum)
names(AGG) = c('Species','Class','Taxonomy','NeutralA','NeutralT','NeutralG','NeutralC')

## how many classes we have
unique(AGG$Class)
table(AGG$Class)

# count fractions of each nucleotide, change light chain notation to heavy one:

AGG$FrTh = AGG$NeutralA / (AGG$NeutralA + AGG$NeutralT + AGG$NeutralG + AGG$NeutralC)
AGG$FrAh = AGG$NeutralT / (AGG$NeutralA + AGG$NeutralT + AGG$NeutralG + AGG$NeutralC) 
AGG$FrCh = AGG$NeutralG / (AGG$NeutralA + AGG$NeutralT + AGG$NeutralG + AGG$NeutralC) 
AGG$FrGh = AGG$NeutralC / (AGG$NeutralA + AGG$NeutralT + AGG$NeutralG + AGG$NeutralC)

median(AGG[AGG$Class == 'Mammalia',]$FrTh)
median(AGG[AGG$Class != 'Mammalia',]$FrTh)

median(AGG[AGG$Class == 'Mammalia',]$FrAh)
median(AGG[AGG$Class != 'Mammalia',]$FrAh)

median(AGG[AGG$Class == 'Mammalia',]$FrGh)
median(AGG[AGG$Class != 'Mammalia',]$FrGh)

median(AGG[AGG$Class == 'Mammalia',]$FrCh)
median(AGG[AGG$Class != 'Mammalia',]$FrCh)

