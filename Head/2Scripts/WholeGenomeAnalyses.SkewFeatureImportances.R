################################
################################

rm(list=ls(all=TRUE))

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")){file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

names(SynNuc)

### make ND6 complementary:
NotND6 = SynNuc[SynNuc$Gene != 'ND6',]
ND6 = SynNuc[SynNuc$Gene == 'ND6',]
A = ND6$NeutralT
T = ND6$NeutralA
G = ND6$NeutralC
C = ND6$NeutralG
ND6$NeutralA = A
ND6$NeutralT = T
ND6$NeutralG = G
ND6$NeutralC = C
SynNuc = rbind(NotND6,ND6)

VecOfTaxa = unique(SynNuc$Class)
table(SynNuc$Class)/13

########################################## GENOME WIDE SKEW FOR EACH SPECIES

AGG = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species,SynNuc$Class), FUN = sum)
names(AGG) = c('Species','Class','NeutralA','NeutralT','NeutralG','NeutralC')

## all six different skews
AGG$CTSkew = (AGG$NeutralC - AGG$NeutralT)/(AGG$NeutralC + AGG$NeutralT); summary(AGG$CTSkew) # GA on heavy
AGG$CGSkew = (AGG$NeutralC - AGG$NeutralG)/(AGG$NeutralC + AGG$NeutralG); summary(AGG$CGSkew) # 
AGG$CASkew = (AGG$NeutralC - AGG$NeutralA)/(AGG$NeutralC + AGG$NeutralA); summary(AGG$CASkew) # 
AGG$TGSkew = (AGG$NeutralT - AGG$NeutralG)/(AGG$NeutralT + AGG$NeutralG); summary(AGG$TGSkew) # 
AGG$TASkew = (AGG$NeutralT - AGG$NeutralA)/(AGG$NeutralT + AGG$NeutralA); summary(AGG$TASkew) # 
AGG$GASkew = (AGG$NeutralG - AGG$NeutralA)/(AGG$NeutralG + AGG$NeutralA); summary(AGG$GASkew) # 

#AGG$TCSkew = (AGG$NeutralT - AGG$NeutralC)/(AGG$NeutralT + AGG$NeutralC); summary(AGG$CTSkew) # AG on heavy. Added it for fun, just to be sure, that it is opposite to CT (GA on heavy) 

GT = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

M = merge(AGG,GT, by ='Species')
ShortLived = unique(M[M$GenerationLength_d <= median(M$GenerationLength_d),]$Species); length(ShortLived)
LongLived = unique(M[M$GenerationLength_d  > median(M$GenerationLength_d),]$Species);  length(LongLived)
MShort = M[M$Species %in% ShortLived,]; MShort$GT = 0; MLong = M[M$Species %in% LongLived,]; MLong$GT = 1;
M = rbind(MShort,MLong)

M_clean = na.omit(M)

br <- Boruta::Boruta(M_clean$GT~., data=M[, c('CTSkew', 'CGSkew', 'CASkew', 'TGSkew', 'TASkew', 'GASkew')], maxRuns=1000)
dev.off()

plot(br)