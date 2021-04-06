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

### estimate GhAhSkew
# from our BioRxiv with mammals:
# An excess of GH and deficit of AH in long-lived species determines the GHAH nucleotide skew which approximates the level of asymmetry in the distribution of these two nucleotides and is calculated as 
# (GH-AH)/(GH+AH). As expected we observed positive correlation between GHAH skew and generation length of mammalian species
# in our case we have light chain notation => 
SynNuc$GhAhSkew = (SynNuc$NeutralC-SynNuc$NeutralT)/(SynNuc$NeutralC+SynNuc$NeutralT)
summary(SynNuc$GhAhSkew)

#### keep only Aves and order by GhAhSkew
table(SynNuc$Class) # 432 Aves
SynNuc = SynNuc[SynNuc$Class =='Aves',]
SynNuc=SynNuc[order(SynNuc$GhAhSkew),]

### derive Dummy Variable - Passeriformes
SynNuc$Passeriformes = 0
for (i in 1:nrow(SynNuc))
{ # i = 1
  SynNuc$Passeriformes[i] = as.numeric(grepl('Passeriformes',SynNuc$Taxonomy[i]))
}
table(SynNuc$Passeriformes)
wilcox.test(SynNuc[SynNuc$Passeriformes == 0,]$GhAhSkew,SynNuc[SynNuc$Passeriformes == 1,]$GhAhSkew)
boxplot(SynNuc[SynNuc$Passeriformes == 0,]$GhAhSkew,SynNuc[SynNuc$Passeriformes == 1,]$GhAhSkew, notch = TRUE, outline = FALSE)

#### read Bird phenotypes
Clutch = read.table("../../Body/2Derived/AvesClutchSizes.txt", header = TRUE, sep = "\t")

### merge MutSpec with PHE by 'Species'
MutSpecPhe = merge(SynNuc,Clutch, by = 'Species')
nrow(MutSpecPhe) # 274

### descriptive analyses
summary(MutSpecPhe$Clutch)
cor.test(MutSpecPhe$GhAhSkew,MutSpecPhe$Clutch, method = 'spearman') # a bit positive
cor.test(MutSpecPhe$FrT,MutSpecPhe$Clutch, method = 'spearman') # a bit negative
cor.test(MutSpecPhe$FrA,MutSpecPhe$Clutch, method = 'spearman') # nothing
