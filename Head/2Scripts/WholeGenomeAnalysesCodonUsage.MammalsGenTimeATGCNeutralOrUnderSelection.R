rm(list=ls(all=TRUE))

###### PURPOSE: derive a table with 649 unique mammalian species with A T G C per whole mit genome, A T G C in neutral sites (with removed overlaps between genes), A T G C in selected sites (whole genome minus neutral) and generation length.

########### 1 read species x gene nucleotide content table:
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip", exdir = "../../Body/3Results/") # AllGenesCodonUsageNoOverlap.txt
NeutralATGC = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")

########### 1A filter out all classes except mammals (788 species x 13 genes)
table(NeutralATGC$Class)
NeutralATGC = NeutralATGC[NeutralATGC$Class == 'Mammalia',]
nrow(NeutralATGC)/13 # 788

########### 1B summ up all neutral A T G C per species
NeutralAtgcSum = aggregate(list(NeutralATGC$NeutralA,NeutralATGC$NeutralT,NeutralATGC$NeutralG,NeutralATGC$NeutralC), by = list(NeutralATGC$Species), FUN = sum)
names(NeutralAtgcSum) = c('species','NeutralA','NeutralT','NeutralG','NeutralC')

########### 2 read A T G C in the whole genome and taxonomy
unzip("../../Body/2Derived/full_table.zip", exdir = "../../Body/2Derived/")
WholeGenomeATGC = read.table("../../Body/2Derived/full_table.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/2Derived/full_table.txt")) file.remove("../../Body/2Derived/full_table.txt")
names(WholeGenomeATGC) = c('species','taxonomy','WholeGenomeA','WholeGenomeC','WholeGenomeG','WholeGenomeT','WholeGenomeX','genome')

########## 3 merge NeutralAtgcSum with WholeGenomeATGC by species
ATGC = merge(NeutralAtgcSum,WholeGenomeATGC, by ='species')
nrow(ATGC) # 788

########## 4 read Generation Length
GL = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GL$species = gsub(' ','_',GL$Scientific_name)
GL = GL[,names(GL) %in% c('species','GenerationLength_d')] 
nrow(GL); length(unique(GL$species)); # 5424, 5426 => two species have two different generation lengthes => take average
GL = aggregate(GL$GenerationLength_d, by = list(GL$species), FUN = mean); names(GL)=c('species','GenerationLength_d')

########## 5 merge ATGC with GL => ATGC
nrow(ATGC); length(unique(ATGC$species))
AtgcGl = merge(ATGC,GL, by = 'species')
nrow(AtgcGl); length(unique(AtgcGl$species)); # 649

########## 6 derive A T G C under selection as whole genome minus neutral

AtgcGl$UnderSelectionA = AtgcGl$WholeGenomeA - AtgcGl$NeutralA
AtgcGl$UnderSelectionT = AtgcGl$WholeGenomeT - AtgcGl$NeutralT
AtgcGl$UnderSelectionG = AtgcGl$WholeGenomeG - AtgcGl$NeutralG
AtgcGl$UnderSelectionC = AtgcGl$WholeGenomeC - AtgcGl$NeutralC

########## 7 write table 
write.table(AtgcGl, "../../Body/3Results/WholeGenomeAnalysesCodonUsage.MammalsGenTimeATGCNeutralOrUnderSelection.R.txt")

