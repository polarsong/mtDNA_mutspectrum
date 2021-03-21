rm(list=ls(all=TRUE))

##### INITIALIZE PARAMETERS:

Final = data.frame()
GenomeLength = 10000  
SimulationLengthNumberOfGenerations = 1000000

##### A: INITIALIZE GENOME

for (MutSpecProb in c('cold-water fish','warm-water fish'))
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
# figure 1B: cold-water fishes: A>G = 0.1; T>C = 0.1
# figure 1B: warm-water fishes: A>G = 0.15; T>C = 0.05
if (MutSpecProb == 'cold-water fish')
{ VecMutSpec = c(
'C','T',0.55, 
'A','G',0.1,
'G','A',0.08,
'T','C',0.1,
'T','G',0.02125,
'T','A',0.02125,
'G','T',0.02125,
'G','C',0.02125,
'C','G',0.02125,
'C','A',0.02125,
'A','T',0.02125,
'A','C',0.02125)}

if (MutSpecProb == 'warm-water fish')
{ VecMutSpec = c(
  'C','T',0.55, 
  'A','G',0.15,
  'G','A',0.08,
  'T','C',0.05,
  'T','G',0.02125,
  'T','A',0.02125,
  'G','T',0.02125,
  'G','C',0.02125,
  'C','G',0.02125,
  'C','A',0.02125,
  'A','T',0.02125,
  'A','C',0.02125)}

MutSpec = data.frame(matrix(VecMutSpec, ncol = 3, nrow = 12, byrow = TRUE))
names(MutSpec) = c('From','To','Prob')
MutSpec$Prob = as.numeric(MutSpec$Prob)
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

### 6: PLOT
# real data for cold-water fishes: A ~ 0.24, C ~ 0.08, G ~ 0.34, T ~ 0.34 # total = 1
# real data for warm-water fishes: A ~ 0.20, C ~ 0.06, G ~ 0.36, T ~ 0.38 # total = 1
# summary(Final$InitGenome)
pdf("../../Body/4Figures/5A.FromMutSpecToNucContent.R.pdf",  width = 20, height = 10)
par(mfrow=c(1,2))
for (InitGenome in 1:10)
{
plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$FrA, ylim = c(0,0.6), pch = '.', col = rgb(1,0,0,0.3), main = '', xlab = '', ylab = ''); par(new=TRUE)
plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$FrT, ylim = c(0,0.6), pch = '.', col = rgb(0,1,0,0.3), main = '', xlab = '', ylab = ''); par(new=TRUE)
plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$FrG, ylim = c(0,0.6), pch = '.', col = rgb(0,0,1,0.3), main = '', xlab = '', ylab = ''); par(new=TRUE)
plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$FrC, ylim = c(0,0.6), pch = '.', col = rgb(0,1,1,0.3), main = '', xlab = '', ylab = ''); par(new=TRUE)
}
legend(1,0.6, legend = c('A','T','G','C'), col = c(rgb(1,0,0,0.3), rgb(0,1,0,0.3), rgb(0,0,1,0.3), rgb(0,1,1,0.3)), pch = 16)

for (InitGenome in 1:10)
{
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$FrA, ylim = c(0,0.6), pch = '.', col = rgb(1,0,0,0.6), main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$FrT, ylim = c(0,0.6), pch = '.', col = rgb(0,1,0,0.6), main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$FrG, ylim = c(0,0.6), pch = '.', col = rgb(0,0,1,0.6), main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$FrC, ylim = c(0,0.6), pch = '.', col = rgb(0,1,1,0.6), main = '', xlab = '', ylab = ''); par(new=TRUE)
}
legend(1,0.6, legend = c('A','T','G','C'), col = c(rgb(1,0,0,0.6), rgb(0,1,0,0.6), rgb(0,0,1,0.6), rgb(0,1,1,0.6)), pch = 16)

dev.off()

### make final boxplots and print out summary (mean, median)

### IF EQUILIBRIUM SENSITIVE TO STARTING CONDITIONS?
### IF EQUILIBRIUM IS SENSITIVE TO MUTSPEC (cold versus warm fishes)
### how the same derive with analogous functions (Lynch, Valerian, Stepan)


