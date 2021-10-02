rm(list=ls(all=TRUE))

library(ggplot2)
library(dplyr)

### take data from Valerian's diffur
expected_cold = data.frame(0.1645046270, 0.2378578262, 0.5271436896, 0.07049385836)
names(expected_cold) =c('FrA','FrG','FrT','FrC')

expected_warm = data.frame(0.1518187869, 0.3236389788, 0.4862288867, 0.03831334710)
names(expected_warm) =c('FrA','FrG','FrT','FrC')

### reverse to put in graph
expected_all = cbind(expected_cold, expected_warm)
expected_all = t(as.matrix(expected_all))

## write types and colnames
expected_all = cbind(expected_all,c('cold_fish','cold_fish','cold_fish','cold_fish','warm_fish','warm_fish','warm_fish','warm_fish'))
expected_all = cbind(expected_all, c('FrA','FrG','FrT','FrC','FrA','FrG','FrT','FrC'))
colnames(expected_all) = c('expected', 'type_of_fish','mutation')

### reading whole genomes database
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
SynNuc = SynNuc[SynNuc$Gene != 'ND6',]

####### obtaining neutral nucleotide fractions in whole genomes
SynNuc = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species), FUN = sum)
names(SynNuc) = c('Species','NeutralA','NeutralT','NeutralG','NeutralC')
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)

SynNuc = SynNuc[,c(1,6,7,8,9)]

##Temperature for fishes
muttemp = read.table('../../Body/2Derived/Supplementary_table_2.txt', header =TRUE) ## ../../

muttemp = muttemp[is.na(muttemp$Temperature) != T,] # delete bad rows

nuc_and_temp = merge(muttemp, SynNuc, by = 'Species') ## take fishes that have codons and temperature

quant = quantile(nuc_and_temp$Temperature, probs = c(0.1, 0.9))

cold = nuc_and_temp[nuc_and_temp$Temperature <= quant[[1]],]

warm = nuc_and_temp[nuc_and_temp$Temperature >= quant[[2]],]

cold_fish = (apply(as.matrix(cold[,17:20]), 2, mean))
warm_fish = (apply(as.matrix(warm[,17:20]), 2, mean))

observed = cbind(t(as.matrix(cold_fish)), t(as.matrix(warm_fish)))
observed = t(as.matrix(observed))

observed = cbind(observed, c('FrT', 'FrA', 'FrC', 'FrG', 'FrT', 'FrA', 'FrC', 'FrG')) ### CHANGE TO HEAVY CHAIN!!!
observed = cbind(observed, c('cold_fish','cold_fish','cold_fish','cold_fish','warm_fish','warm_fish','warm_fish','warm_fish'))
colnames(observed) = c('observed','mutation','type_of_fish')

expectedvsoberved = merge(expected_all, observed, by = c('mutation','type_of_fish'))
expectedvsoberved$observed =as.numeric(as.character(expectedvsoberved$observed))
expectedvsoberved$expected =as.numeric(as.character(expectedvsoberved$expected))

pdf("../../Body/4Figures/5c.ExpectedvsObserved.pdf")

exvsobs = ggplot(data = expectedvsoberved, aes(x = observed, y = expected, group=type_of_fish, col = type_of_fish))+
  geom_point(size = 3.5)+
  geom_abline(col = 'gray3', linetype="longdash", size = 0.6)+
  theme_bw()+
  scale_color_manual(name="Type of Fish", labels = c('Cold fish', 'Warm fish'), values = c('cold_fish'='deepskyblue4', 'warm_fish' = 'firebrick3'))+
  geom_line(aes(group = mutation), col = 'black', size = 0.7)+
  geom_text(aes(label=mutation),hjust=-0.40, vjust=-0.40)+
  xlim(0, 0.55)+
  ylim(0,0.55)+
  labs(x = 'Observed Nucleotide Content',y = 'Expected Nucleotide Content')
  

print(exvsobs)

dev.off()

