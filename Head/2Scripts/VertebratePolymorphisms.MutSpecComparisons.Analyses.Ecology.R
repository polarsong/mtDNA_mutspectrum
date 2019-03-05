rm(list=ls(all=TRUE))

MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)

GT = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', sep = '\t', header = TRUE)
GT$Species = gsub(' ','_',GT$Scientific_name)

#########
#### Generation length in mammals and mut spec
#########

MamGt = merge(GT,MUT, by = 'Species')

cor.test(MamGt$GenerationLength_d,MamGt$T_C, method = 'spearman')
