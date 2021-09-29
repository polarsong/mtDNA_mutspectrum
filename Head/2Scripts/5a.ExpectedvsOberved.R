rm(list = ls(all=TRUE))

library(caper)
library(geiger)

### reading whole genomes database 
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}
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

#### take 10% of cold and warm fishes
quant = quantile(nuc_and_temp$Temperature, probs = c(0.1, 0.9))

cold = nuc_and_temp[nuc_and_temp$Temperature <= quant[[1]],]

warm = nuc_and_temp[nuc_and_temp$Temperature >= quant[[2]],]


cold_fish = apply(as.matrix(cold[,2:13]), 2, mean)
warm_fish = apply(as.matrix(warm[,2:13]), 2, mean)

names(cold_fish) = c('T_A','T_C','T_G','A_T','A_C','A_G','C_T','C_A','C_G','G_T','G_A','G_C')
names(warm_fish) = c('T_A','T_C','T_G','A_T','A_C','A_G','C_T','C_A','C_G','G_T','G_A','G_C')

### check mean_mut == 1 in both
sum(cold_fish[1:12], na.rm = T)
sum(warm_fish[1:12], na.rm=T)

##### INITIALIZE PARAMETERS:

Final = data.frame()
GenomeLength = 10000  
SimulationLengthNumberOfGenerations = 1000000

##### A: INITIALIZE GENOME

for (MutSpecProb in c('cold-water fish','warm-water fish'))
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
    # Take MutSpec from the two values: cold_fish and worm_fish
    
    if (MutSpecProb == 'cold-water fish')
    { #cold_fish
      VecMutSpec = c(  
      "T","A",cold_fish[['T_A']],
      "T","C",cold_fish[['T_C']], 
      "T","G",cold_fish[['T_G']],
      "A","T",cold_fish[['A_T']],
      "A","C",cold_fish[['A_C']],
      "A","G",cold_fish[['A_G']], 
      "C","T",cold_fish[['C_T']],
      "C","A",cold_fish[['C_A']],
      "C","G",cold_fish[['C_G']],
      "G","T",cold_fish[['G_T']],
      "G","A",cold_fish[['G_A']],
      "G","C",cold_fish[['G_C']])} 
    
    if (MutSpecProb == 'warm-water fish')
    { #warm_fish
      VecMutSpec = c(  
        "T","A",warm_fish[['T_A']],
        "T","C",warm_fish[['T_C']], 
        "T","G",warm_fish[['T_G']],
        "A","T",warm_fish[['A_T']],
        "A","C",warm_fish[['A_C']],
        "A","G",warm_fish[['A_G']], 
        "C","T",warm_fish[['C_T']],
        "C","A",warm_fish[['C_A']],
        "C","G",warm_fish[['C_G']],
        "G","T",warm_fish[['G_T']],
        "G","A",warm_fish[['G_A']],
        "G","C",warm_fish[['G_C']])} 
    
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

write.table(Final, "../../Body/Body/2Derived/5A.FromMutSpecToNucContent.R.FinalTable.txt")



