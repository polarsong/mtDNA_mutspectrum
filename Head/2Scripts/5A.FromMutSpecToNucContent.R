rm(list=ls(all=TRUE))

##### INITIALIZE PARAMETERS:

Final = data.frame()
GenomeLength = 10000  
SimulationLengthNumberOfGenerations = 1000000

##### A: INITIALIZE GENOME

for (MutSpecProb in c('average fish','cold-water fish','warm-water fish'))
{
### choose initial nucleotide frequencies: equal to 25% if InitGenome == 1 or random if InitGenome > 1
for (InitGenome in 1:10)
{  # InitGenome = 2
  if (InitGenome == 1) {frA = frG = frC = frT = 0.25} # frA = 0.145342343
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
# Take MutSpec from this file:
# https://github.com/polarsong/mtDNA_mutspectrum/blob/TemperatureVSVertabrates/Body/2Derived/ActinopterMutSpec.txt
# and modify +/- 2.5%  for A>G and T>C for cold- and warm-blooded (5% difference is taken visually from figure 1B)

if (MutSpecProb == 'average fish')
{ VecMutSpec = c(  
  "T","A",0.009953,
  "T","C",0.07020, # add 2.5% for cold-blooded (0.07020 + 0.025 = 0.0952) and extract 2.5% for warm-blooded (0.07020 - 0.025 = 0.0452)
  "T","G",0.013630,
  "A","T",0.02742,
  "A","C",0.008847,
  "A","G",0.1422, # extract 2.5% for cold-blooded (0.1422 - 0.025 = 0.1172) and add 2.5 for warm-blooded (0.1422 + 0.025 = 0.1672)
  "C","T",0.536696, # decreased it by 0.000004 to make a sum == 1 (0.5367 - 0.000004 = 0.536696)
  "C","A",0.025616,
  "C","G",0.05079,
  "G","T",0.026663,
  "G","A",0.07799,
  "G","C",0.009995)} 


if (MutSpecProb == 'cold-water fish')
{ VecMutSpec = c(  
"T","A",0.009953,
"T","C",0.0952, # add 2.5% for cold-blooded (0.07020 + 0.025 = 0.0952) and extract 2.5% for warm-blooded (0.07020 - 0.025 = 0.0452)
"T","G",0.013630,
"A","T",0.02742,
"A","C",0.008847,
"A","G",0.1172, # extract 2.5% for cold-blooded (0.1422 - 0.025 = 0.1172) and add 2.5 for warm-blooded (0.1422 + 0.025 = 0.1672)
"C","T",0.536696, # decreased it by 0.000004 to make a sum == 1 (0.5367 - 0.000004 = 0.536696)
"C","A",0.025616,
"C","G",0.05079,
"G","T",0.026663,
"G","A",0.07799,
"G","C",0.009995)} 

if (MutSpecProb == 'warm-water fish')
{ VecMutSpec = c(  
"T","A",0.009953,
"T","C",0.0452,   # add 2.5% for cold-blooded (0.07020 + 0.025 = 0.0952) and extract 2.5% for warm-blooded (0.07020 - 0.025 = 0.0452)
"T","G",0.013630,
"A","T",0.02742,
"A","C",0.008847,
"A","G",0.1672,   # extract 2.5% for cold-blooded (0.1422 - 0.025 = 0.1172) and add 2.5 for warm-blooded (0.1422 + 0.025 = 0.1672)
"C","T",0.536696, # decreased it by 0.000004 to make a sum == 1 (0.5367 - 0.000004 = 0.536696)
"C","A",0.025616,
"C","G",0.05079,
"G","T",0.026663,
"G","A",0.07799,
"G","C",0.009995)} 
  
MutSpec = data.frame(matrix(VecMutSpec, ncol = 3, nrow = 12, byrow = TRUE))
names(MutSpec) = c('From','To','Prob')
MutSpec$Prob = as.numeric(MutSpec$Prob)
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
ExpectedFrA # 0.09887561
ExpectedFrT # 0.6645893
ExpectedFrG # 0.21168
ExpectedFrC # 0.02485505

sum(MutSpec$Prob)
for (i in 1:nrow(MutSpec))
{ # i = 1
if (i == 1) {MutSpec$RulletFrom[i] = 0}
MutSpec$RulletTo[i] = sum(MutSpec[seq(1:i),]$Prob)
if (i  > 1) MutSpec$RulletFrom[i] = MutSpec$RulletTo[i-1]
}

##### C: MUTATE AND SAVE NUCLEOTIDE CONTENT EVERY 100 GENERATIONS
for (gener in 1:SimulationLengthNumberOfGenerations)
{
### 3: choose a random position in a genome
 RandomPos = sample(1:length(genome), 1)  
 NucInRandomPos = genome[RandomPos];
 
### 2: choose a random mutation
 Rullet = runif(1)
 PotentialSubstitution = MutSpec[MutSpec$RulletFrom <= Rullet & MutSpec$RulletTo > Rullet,]
 
### 3: mutation happens
 if (nrow(PotentialSubstitution) == 1)
 {
  if (NucInRandomPos == PotentialSubstitution$From)  {genome[RandomPos] = PotentialSubstitution$To}
 }

### 4: print out every 100 generations
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

### 5: DERIVE FRACTIONS AND SAVE

Final$FrA = as.numeric(Final$A)/length(genome)
Final$FrT = as.numeric(Final$T)/length(genome)
Final$FrG = as.numeric(Final$G)/length(genome)
Final$FrC = as.numeric(Final$C)/length(genome)

write.table(Final, "../../Body/2Derived/5A.FromMutSpecToNucContent.R.FinalTable.txt")

