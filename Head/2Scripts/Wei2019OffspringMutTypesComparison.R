rm(list=ls())

data = read.table("../../Body/1Raw/Wei2019HumanDuos.txt", sep='\t', header = TRUE)
str(data)
length(unique(data$PairID)) # 928
table(data$MotherOffspring) # Mother: 1043 Offspring: 893
table(data$HeteroplasmiCategory) # de novo: 477; inherited: 416; lost: 614; tansmitted: 429  # inherited - transmitted&&&

### Add mutation type column
data$Allele = as.character(data$Allele)
data$mutType = paste(substr(data$Allele, 1, 1), substr(data$Allele, nchar(data$Allele), nchar(data$Allele)), sep="")

### we need to work with offspring who has several de novo variants (to compare their VAFs)
Off = data[which(data$MotherOffspring == 'Offspring'),]
table(Off$HeteroplasmiCategory) # de novo: 477; inherited: 416  
OffDeNovo = Off[Off$HeteroplasmiCategory == 'de novo',]  # 477

### check that this is indeed light chain notation (default as expected) and check - how reasonable the MutSpec in general
TotalMutSpecLightChain = as.data.frame(table(OffDeNovo$mutType)); TotalMutSpecLightChain = TotalMutSpecLightChain[order(-TotalMutSpecLightChain$Freq),]
TotalMutSpecLightChain
# Var1 Freq
# 10   TC  159 # our A>G. it is the most common here (not the second most common as in cancers and in mammals)! unusual a bit, but from other side - they are oocytes!!!!!
# 7    GA  129 # C>T
# 2    AG  114 # T>C
# 6    CT   52 # 
# 1    AC    7
# 4    CA    4
# 11   TG    4
# 3    AT    3
# 5    CG    2
# 9    TA    2
# 8    GC    1

# invert transitions to heavy notations (to A>G) and collapse all transversions to Tv:
OffDeNovo$mutType = gsub('TC','A>G',OffDeNovo$mutType);  OffDeNovo$mutType = gsub('GA','C>T',OffDeNovo$mutType); OffDeNovo$mutType = gsub('AG','T>C',OffDeNovo$mutType); OffDeNovo$mutType = gsub('CT','G>A',OffDeNovo$mutType); 
OffDeNovo$mutType = gsub("AC|CA|GT|TG|AT|TA|CG|GC", 'Tv',OffDeNovo$mutType);
TotalMutSpecHeavyChain = as.data.frame(table(OffDeNovo$mutType)); TotalMutSpecHeavyChain = TotalMutSpecHeavyChain[order(-TotalMutSpecHeavyChain$Freq),]
TotalMutSpecHeavyChain
# Var1 Freq
#1  A>G  159
#2  C>T  129
#4  T>C  114
#3  G>A   52
#5   Tv   23

#### average heteroplasmy per substitution type: G>A, A>G, Tv, T>C, C>T
Agg1 = aggregate(OffDeNovo$HeteroplasmicFraction, by = list(OffDeNovo$mutType), FUN = mean)
Agg1

boxplot(OffDeNovo$HeteroplasmicFraction ~ OffDeNovo$mutType, varwidth = TRUE, notch = TRUE)


#### OLD 


#### get paired data: average heteroplasmy for A>G and average heteroplasmy for all others
# offspring$mutType1 = offspring$mutType
offspring$mutType1 = gsub("A>G","X",offspring$mutType); offspring$mutType1 = gsub("[^[X]+","others",offspring$mutType1);   # not X one or more times

Agg = aggregate(offspring$HeteroplasmicFraction, by = list(offspring$PairID,offspring$mutType1), FUN = mean)
Agg1 = Agg[Agg$Group.2 == 'X',]; Agg1 = Agg1[,c(1,3)]; names(Agg1) = c('pairID','HeteroplasmyA2G')
Agg2 = Agg[Agg$Group.2 == 'others',]; Agg2 = Agg2[,c(1,3)]; names(Agg2) = c('pairID','HeteroplasmyOthers')
Agg = merge(Agg1,Agg2)
Agg$heteroplasmyA2G2others = log2(Agg$HeteroplasmyA2G / Agg$HeteroplasmyOthers)
summary(Agg$heteroplasmyA2G2others)
wilcox.test(Agg$heteroplasmyA2G2others, mu = 0, alternative = 'less')
hist(Agg$heteroplasmyA2G2others, breaks = 50)

#### if we repeat the same test but focus on more and more rare variants => do we see the same trend? yes, we see, but how we can prove that this is correct analysis? we use heteroplasmy as both - function of the time and neutrality... not very good
offspring = offspring[offspring$HeteroplasmicFraction < 10,]

Agg = aggregate(offspring$HeteroplasmicFraction, by = list(offspring$PairID,offspring$mutType1), FUN = mean)
Agg1 = Agg[Agg$Group.2 == 'X',]; Agg1 = Agg1[,c(1,3)]; names(Agg1) = c('pairID','HeteroplasmyA2G')
Agg2 = Agg[Agg$Group.2 == 'others',]; Agg2 = Agg2[,c(1,3)]; names(Agg2) = c('pairID','HeteroplasmyOthers')
Agg = merge(Agg1,Agg2)
Agg$heteroplasmyA2G2others = log2(Agg$HeteroplasmyA2G / Agg$HeteroplasmyOthers)
summary(Agg$heteroplasmyA2G2others)
wilcox.test(Agg$heteroplasmyA2G2others, mu = 0, alternative = 'less')
hist(Agg$heteroplasmyA2G2others, breaks = 50)



#### number of de novo mutations per sample (too many for some individuals is not good) - looks good
offspring$MutLoad = 1; 
Agg3 = aggregate(offspring$MutLoad, by = list(offspring$PairID), FUN = sum)
table(Agg3$x)
#  1   2   3   4 
# 454 152  33   9 


# BELOW IS OLD CODE BY KRISTINA & ALYA


  
offspringMultiple = offspring[offspring$PairID %in% offspring$PairID[duplicated(offspring$PairID)],][order(offspring$PairID),]
length(unique(offspringMultiple$PairID))

# Only A>G and T>C
offspringMultipleAGTC = offspringMultiple[which(offspringMultiple$mutType == "A>G" | offspringMultiple$mutType == "T>C"),]

# Only G>A and C>T
offspringMultipleGACT = offspringMultiple[which(offspringMultiple$mutType == "G>A" | offspringMultiple$mutType == "C>T"),]


### ANALYSES ###
# 1.  All mutation types
# 1.1 Scatterplot. X: PairID, Y:heteroplasmy. Colour by type
offspringMultiple %>%
ggplot(aes(x=PairID, y=HeteroplasmicFraction, color=mutType)) + 
         geom_point()

# 1.2 Boxplot by mutation types
ggplot(aes(x=mutType, y=HeteroplasmicFraction, fill=mutType)) +
  geom_boxplot()

# 2.  Only A>G and T>C
# 2.1 Scatterplot. Only A>G and T>C
offspringMultipleAGTC %>%
ggplot(aes(x=PairID, y=HeteroplasmicFraction, color=mutType)) + 
  geom_point()

# 2.2 Boxplot by mutation types. Only A>G and T>C
offspringMultipleAGTC %>%
ggplot(aes(x=mutType, y=HeteroplasmicFraction, fill=mutType)) +
  geom_boxplot()

# 3.  Only G>A and C>T
# 3.1 Scatterplot. Only G>A and C>T
offspringMultipleGACT %>%
  ggplot(aes(x=PairID, y=HeteroplasmicFraction, color=mutType)) + 
  geom_point()

# 3.2 Boxplot by mutation types.
offspringMultipleGACT %>%
ggplot(aes(x=mutType, y=HeteroplasmicFraction, fill=mutType)) +
  geom_boxplot()

# 4. Regression. Only A>G and T>C
a <- lm(offspringMultiple$HeteroplasmicFraction ~ offspringMultiple$mutType); summary(a)

