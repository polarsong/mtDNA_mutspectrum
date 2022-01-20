rm(list=ls(all=TRUE))

library(dplyr)
library(ggplot2)

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

##Generation Time of Mammals
gen_len = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', header =TRUE, sep = '\t') ## ../../
gen_len = gen_len[ ,c(2,11)] 
gen_len$Scientific_name = gsub(' ','_',gen_len$Scientific_name)


## Read polymorphism

poly = read.table('../../Body/3DerivedVertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE) 

polygen = merge(poly, gen_len, by.x = 'Species', by.y = 'Scientific_name')

nuc_and_gen = merge(polygen, SynNuc, by = 'Species') ## take fishes that have codons and temperature

#### take 10% of long and short mammals
quant = quantile(nuc_and_gen$GenerationLength_d, probs = c(0.1, 0.9))

short = nuc_and_gen[nuc_and_gen$GenerationLength_d <= quant[[1]],]

long = nuc_and_gen[nuc_and_gen$GenerationLength_d >= quant[[2]],]


short_mam = apply(as.matrix(short[,2:13]), 2, mean)
long_mam = apply(as.matrix(long[,2:13]), 2, mean)

names(short_mam) = c('T_A','T_C','T_G','A_T','A_C','A_G','C_T','C_A','C_G','G_T','G_A','G_C')
names(long_mam) = c('T_A','T_C','T_G','A_T','A_C','A_G','C_T','C_A','C_G','G_T','G_A','G_C')

### check mean_mut == 1 in both
sum(short_mam[1:12], na.rm = T)
sum(long_mam[1:12], na.rm=T)

##### INITIALIZE PARAMETERS:

Final = data.frame()
GenomeLength = 10000  
SimulationLengthNumberOfGenerations = 1000000

##### A: INITIALIZE GENOME

for (MutSpecProb in c('short-matu mam','long-matu mam'))
{
  ### choose initial nucleotide frequencies: equal to 25% if InitGenome == 1 or random if InitGenome > 1
  for (InitGenome in 1:10)
  {  #InitGenome = 2
    if (InitGenome == 1) {frA = frG = frC = frT = 0.25} 
    if (InitGenome >  1) 
    {
      frA = runif(1); frG = runif(1);
      frC = runif(1); frT = runif(1);  
      Summa = frA+frT+frG+frC
      frA = frA/Summa; frG = frG/Summa; frC = frC/Summa; frT = frT/Summa; 
      frA+frG+frC+frT # should be 1
    }
    
    ### make a genome
    genome = sample(c(rep('A',round(frA*GenomeLength)),rep('T',round(frT*GenomeLength)),rep('G',round(frG*GenomeLength)),rep('C',round(frC*GenomeLength))))
    length(genome) # should be equal GenomeLength
    
    ##### B: DEFINE MUTATIONAL SPECTRUM
    # Take MutSpec from the two values: short_mam and worm_fish
    
    if (MutSpecProb == 'short-matu mam')
    { #long maturated mammals
      VecMutSpec = c(  
        "T","A",short_mam[['T_A']],
        "T","C",short_mam[['T_C']], 
        "T","G",short_mam[['T_G']],
        "A","T",short_mam[['A_T']],
        "A","C",short_mam[['A_C']],
        "A","G",short_mam[['A_G']], 
        "C","T",short_mam[['C_T']],
        "C","A",short_mam[['C_A']],
        "C","G",short_mam[['C_G']],
        "G","T",short_mam[['G_T']],
        "G","A",short_mam[['G_A']],
        "G","C",short_mam[['G_C']])} 
    
    if (MutSpecProb == 'long-matu mam')
    { #long_mam
      VecMutSpec = c(  
        "T","A",long_mam[['T_A']],
        "T","C",long_mam[['T_C']], 
        "T","G",long_mam[['T_G']],
        "A","T",long_mam[['A_T']],
        "A","C",long_mam[['A_C']],
        "A","G",long_mam[['A_G']], 
        "C","T",long_mam[['C_T']],
        "C","A",long_mam[['C_A']],
        "C","G",long_mam[['C_G']],
        "G","T",long_mam[['G_T']],
        "G","A",long_mam[['G_A']],
        "G","C",long_mam[['G_C']])} 
    
    MutSpec = data.frame(matrix(VecMutSpec, ncol = 3, nrow = 12, byrow = TRUE))
    names(MutSpec) = c('From','To','Prob')
    MutSpec$Prob = as.numeric(as.character(MutSpec$Prob))
    sum(MutSpec$Prob)
    
    ExpectedFrA = sum(MutSpec[MutSpec$To == 'A',]$Prob)/sum(MutSpec[MutSpec$From == 'A',]$Prob)
    ExpectedFrT = sum(MutSpec[MutSpec$To == 'T',]$Prob)/sum(MutSpec[MutSpec$From == 'T',]$Prob)
    ExpectedFrG = sum(MutSpec[MutSpec$To == 'G',]$Prob)/sum(MutSpec[MutSpec$From == 'G',]$Prob)
    ExpectedFrC = sum(MutSpec[MutSpec$To == 'C',]$Prob)/sum(MutSpec[MutSpec$From == 'C',]$Prob)
    Summa = ExpectedFrA + ExpectedFrT + ExpectedFrG + ExpectedFrC
    ExpectedFrA = ExpectedFrA/Summa
    ExpectedFrT = ExpectedFrT/Summa 
    ExpectedFrG = ExpectedFrG/Summa
    ExpectedFrC = ExpectedFrC/Summa
    MutSpecProb
    ExpectedFrA 
    ExpectedFrT 
    ExpectedFrG 
    ExpectedFrC 
    
    sum(MutSpec$Prob) 
    
    for (i in 1:nrow(MutSpec))
    { # i = 1
      if (i == 1) {MutSpec$RulletFrom[i] = 0}
      MutSpec$RulletTo[i] = sum(MutSpec[seq(1:i),]$Prob)
      if (i  > 1) {MutSpec$RulletFrom[i] = MutSpec$RulletTo[i-1]}
    }
    
    ##### C: MUTATE AND SAVE NUCLEOTIDE CONTENT EVERY 1000 GENERATIONS
    for (gener in 1:SimulationLengthNumberOfGenerations)
    {
      ### 3: choose a random position in a genome
      RandomPos = sample(1:length(genome), 1)  
      NucInRandomPos = genome[RandomPos];
      
      ### 2: choose a random mutation
      Rullet = runif(1)
      PotentialSubstitution = MutSpec[MutSpec$RulletFrom <= Rullet & MutSpec$RulletTo > Rullet,]
      PotentialSubstitution$To = as.character(PotentialSubstitution$To)
      ### 3: mutation happens
      if (nrow(PotentialSubstitution) == 1)
      {
        if (NucInRandomPos == PotentialSubstitution$From)  {genome[RandomPos] = PotentialSubstitution$To}
      }
      
      ### 4: print out every 1000 generations
      if ((gener %% 1000) == 0)
      {
        Res = data.frame(table(genome))
        Res = data.frame(t(Res[order(Res$genome),]))
        Res = Res[2,]
        Res$Gener = gener
        Res$InitGenome = InitGenome
        Res$MutSpecProb = MutSpecProb
        names(Res)= c('A','C','G','T','Gener','InitGenome','MutSpecProb')
        Final = rbind(Final,Res)
      }
    }}}

row.names(Final) = 1:nrow(Final)
Final$FrA = as.numeric(as.character(Final$A))/length(genome)
Final$FrT = as.numeric(as.character(Final$T))/length(genome)
Final$FrG = as.numeric(as.character(Final$G))/length(genome)
Final$FrC = as.numeric(as.character(Final$C))/length(genome)

write.table(Final, "/../../Body/3Results/EquilibriumMammals.txt")

expected = Final %>%
  group_by(MutSpecProb) %>% 
  mutate(mFrA = median(FrA)) %>%
  mutate(mFrT = median(FrT)) %>% 
  mutate(mFrG = median(FrG)) %>%
  mutate(mFrC = median(FrC)) %>% 
  select(MutSpecProb, mFrA, mFrT, mFrG, mFrC) %>% 
  unique()

### take data from simulation 
expected_short = data.frame(0.1643164, 0.2374237, 0.5838584, 0.01540154)
names(expected_short) =c('FrA','FrG','FrT','FrC')

expected_long = data.frame(0.1117112, 0.2582258, 0.5892589, 0.04000400)
names(expected_long) =c('FrA','FrG','FrT','FrC')

### reverse to put in graph
expected_all = cbind(expected_short, expected_long)
expected_all = t(as.matrix(expected_all))

## write types and colnames
expected_all = cbind(expected_all,c('short_mam','short_mam','short_mam','short_mam','long_mam','long_mam','long_mam','long_mam'))
expected_all = cbind(expected_all, c('FrA','FrG','FrT','FrC','FrA','FrG','FrT','FrC'))
colnames(expected_all) = c('expected', 'type_of_mam','mutation')


short_mam = (apply(as.matrix(short[,19:22]), 2, mean))
long_mam = (apply(as.matrix(long[,19:22]), 2, mean))
observed = cbind(t(as.matrix(short_mam)), t(as.matrix(long_mam)))
observed = t(as.matrix(observed))

observed = cbind(observed, c('FrT', 'FrA', 'FrC', 'FrG', 'FrT', 'FrA', 'FrC', 'FrG')) ### CHANGE TO HEAVY CHAIN!!!
observed = cbind(observed, c('short_mam','short_mam','short_mam','short_mam','long_mam','long_mam','long_mam','long_mam'))
colnames(observed) = c('observed','mutation','type_of_mam')

expvsobs = merge(expected_all, observed, by = c('mutation','type_of_mam'))
expvsobs$observed =as.numeric(as.character(expvsobs$observed))
expvsobs$expected =as.numeric(as.character(expvsobs$expected))

ggplot(data = expvsobs, aes(x = observed, y = expected, group=type_of_mam, col = type_of_mam))+
  geom_point(size = 3.5)+
  geom_abline(col = 'gray3', linetype="longdash", size = 0.6)+
  theme_bw()+
  scale_color_manual(name="Type of Mammalia", labels = c('Short Maturated Mammalia', 'Long Maturated Mammalia'), values = c('short_mam'='deepskyblue4', 'long_mam' = 'firebrick3'))+
  geom_line(aes(group = mutation), col = 'black', size = 0.7)+
  geom_text(aes(label=mutation),hjust=-0.40, vjust=-0.40)+
  xlim(0, 0.60)+
  ylim(0,0.60)+
  labs(x = 'Observed Nucleotide Content',y = 'Expected Nucleotide Content')

