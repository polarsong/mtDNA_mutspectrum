# Clear all vars
rm(list=ls())

# Load libraries
require(stringr)
require(ggplot2)
require(magrittr)


### DERIVE ###
# Load the raw data
data = read.table("~/Desktop/Wei2019HumanDuos.txt", sep='\t', header=TRUE)

# Convert contents of allele type column into character values
data$Allele = as.character(data$Allele)

# Add mutation type column
data$mutType = paste(substr(data$Allele, 1, 1), substr(data$Allele, nchar(data$Allele), nchar(data$Allele)), sep=">")

# Only offspring with several mutations
offspring = data[which(data$MotherOffspring == 'Offspring'),]
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

