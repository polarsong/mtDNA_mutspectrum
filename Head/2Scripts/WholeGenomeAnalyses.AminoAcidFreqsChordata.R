rm(list=ls(all=TRUE))

library(dplyr)
library(Biostrings)

unzip('../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip', 
          exdir = '../../Body/3Results/')

codons = read.table('../../Body/3Results/AllGenesCodonUsageNoOverlap.txt', header = TRUE, sep='\t')

codons$AminoNoOverlap = as.character(codons$AminoNoOverlap)
codons$Species = as.character(codons$Species)
codons$Gene = as.character(codons$Gene)

aa = codons %>%
  select(Species, AminoNoOverlap, Taxonomy, Class, Gene)


# as.data.frame(t(alphabetFrequency(AAString(aa$AminoNoOverlap[1000]))))
# as.data.frame(t(alphabetFrequency(AAString(aa$AminoNoOverlap[100]))))

df <- data.frame(matrix(ncol = 33, nrow = 0))

for(i in 1:nrow(aa)){
  # i = 1
  seq = AAString(aa$AminoNoOverlap[i])
  temp_df = as.data.frame(t(c(aa$Species[i], aa$Gene[i], t(alphabetFrequency(seq)))))
  df = rbind(df, temp_df)
}

names(df) = c('Species', 'Gene', names(as.data.frame(t(alphabetFrequency(seq)))))
three_letter = c('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val')

names(df)[3:22] = three_letter

final = distinct(inner_join(df, codons[, c('Species', 'Class')]))

write.csv(final, '~/Alina/Kostya/AminoAcidFreqsChordata.csv', quote = FALSE, row.names = FALSE)
