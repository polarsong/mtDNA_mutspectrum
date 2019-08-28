# Clear all vars
rm(list=ls())

# Load libraries  - V ZADNITCY 
#require(stringr)
#require(ggplot2)
#require(magrittr)

### DERIVE ###
# Load the raw data
# data = read.table("~/Desktop/Wei2019HumanDuos.txt", sep='\t', header=TRUE) # KAKOI V ZADNITCY DESKTOP!!!!! SEE LINE BELOW
data = read.table("../../Body/1Raw/Wei2019HumanDuos.txt", sep='\t', header = TRUE)

# Convert contents of allele type column into character values
str(data)
data$Allele = as.character(data$Allele)

# Add mutation type column
data$mutType = paste(substr(data$Allele, 1, 1), substr(data$Allele, nchar(data$Allele), nchar(data$Allele)), sep="")

### Only offspring with several de novo mutations
offspring = data[which(data$MotherOffspring == 'Offspring'),]

# check frequencies => indeed this is light chain notation
TotalMutSpec = as.data.frame(table(offspring$mutType)); TotalMutSpec = TotalMutSpec[order(-TotalMutSpec$Freq),]
#    Var1 Freq
#10  T>C  311
#7   G>A  210
#2   A>G  194
#6   C>T  121
#1   A>C   24

# invert transitions to heavy notations (to A>G):
offspring$mutType = gsub('TC','A>G',offspring$mutType);  offspring$mutType = gsub('GA','C>T',offspring$mutType); offspring$mutType = gsub('AG','T>C',offspring$mutType); offspring$mutType = gsub('CT','G>A',offspring$mutType); 
offspring$mutType = gsub("AC|CA|GT|TG|AT|TA|CG|GC", 'Tv',offspring$mutType);
TotalMutSpec = as.data.frame(table(offspring$mutType)); TotalMutSpec = TotalMutSpec[order(-TotalMutSpec$Freq),]
TotalMutSpec
#     Var1 Freq
#1    A>G  311
#2    C>T  210
#5    T>C  194
#3    G>A  121
#4others   57

#### average heteroplasmy per substitution type: G>A, A>G, Tv, T>C, C>T
Agg1 = aggregate(offspring$HeteroplasmicFraction, by = list(offspring$mutType), FUN = mean)
Agg1
#1     A>G 21.11029
#2     C>T 15.00333
#3     G>A 25.02893
#4      Tv 20.28246
#5     T>C 16.22423

boxplot(offspring$HeteroplasmicFraction ~ offspring$mutType, varwidth = TRUE, notch = TRUE)

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

