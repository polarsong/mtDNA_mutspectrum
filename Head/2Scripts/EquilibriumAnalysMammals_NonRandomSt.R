rm(list=ls(all=T))

SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')

gen_len = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', header =TRUE, sep = '\t') ## ../../
gen_len = gen_len[ ,c(2,11)] 
gen_len$Scientific_name = gsub(' ','_',gen_len$Scientific_name)

poly = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE) 


SynNuc = SynNuc[SynNuc$Gene != 'ND6',] # take 12 genes except ND6


SynNuc = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species), FUN = sum)
names(SynNuc) = c('Species','NeutralA','NeutralT','NeutralG','NeutralC')
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)



polygen = merge(poly, gen_len, by.x = 'Species', by.y = 'Scientific_name')

nuc_and_gen = merge(polygen, SynNuc, by = 'Species') ## take mammals that have codons and genlen


quant = quantile(nuc_and_gen$GenerationLength_d, probs = c(0.1, 0.9))

short = nuc_and_gen[nuc_and_gen$GenerationLength_d <= quant[[1]],]
long = nuc_and_gen[nuc_and_gen$GenerationLength_d >= quant[[2]],]


short_mam = apply(as.matrix(short[,2:13]), 2, mean)
long_mam = apply(as.matrix(long[,2:13]), 2, mean)

# HEAVY CHAIN!
names(short_mam) = c('T_A','T_C','T_G','A_T','A_C','A_G','C_T','C_A','C_G','G_T','G_A','G_C')
names(long_mam) = c('T_A','T_C','T_G','A_T','A_C','A_G','C_T','C_A','C_G','G_T','G_A','G_C')


short_mam_vec_fr = (apply(as.matrix(short[,19:22]), 2, mean)) #T A C G annot HEAVY CHAIN!
long_mam_vec_fr = (apply(as.matrix(long[,19:22]), 2, mean))


##### INITIALIZE PARAMETERS:

Final = data.frame()
GenomeLength = 10000  
SimulationLengthNumberOfGenerations = 1000000

##### A: INITIALIZE GENOME

for (MutSpecProb in c('short-matu mam','long-matu mam'))
{
  ### choose initial nucleotide frequencies 
  for (InitGenome in 1:10)
  {  #InitGenome = 2
    
    
    if (MutSpecProb == 'short-matu mam') 
    {
      frA = long_mam_vec_fr[2]; frG = long_mam_vec_fr[4];
      frC = long_mam_vec_fr[3]; frT = long_mam_vec_fr[1];  
    }
    else if (MutSpecProb == 'long-matu mam')
    {
      frA = short_mam_vec_fr[2]; frG = short_mam_vec_fr[4];
      frC = short_mam_vec_fr[3]; frT = short_mam_vec_fr[1];
    }
    
    
    
    ### make a genome
    genome = sample(c(rep('A',round(frA*GenomeLength)),rep('T',round(frT*GenomeLength)),rep('G',round(frG*GenomeLength)),rep('C',round(frC*GenomeLength))))
    length(genome) # should be equal GenomeLength
    
    ##### B: DEFINE MUTATIONAL SPECTRUM
    # Take MutSpec from the two values: short_mam and long_mam
    
    if (MutSpecProb == 'short-matu mam')
    { #short maturated mammals
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
        # Check if we have all 4 subs
        if (ncol(Res) == 7){
          names(Res)= c('A','C','G','T','Gener','InitGenome', 'MutSpecProb')
          Final = rbind(Final,Res)
        }
      }}}
}


colnames(Final)[7] = 'MutSpecProb' 
row.names(Final) = 1:nrow(Final)
Final$FrA = as.numeric(as.character(Final$A))/length(genome)
Final$FrT = as.numeric(as.character(Final$T))/length(genome)
Final$FrG = as.numeric(as.character(Final$G))/length(genome)
Final$FrC = as.numeric(as.character(Final$C))/length(genome)


ColorA = rgb(1,0,0,0.3)
ColorT = rgb(0,1,0,0.3)
ColorG = rgb(0,0,1,0.3)
ColorC = rgb(0,1,1,0.3)

par(mfrow=c(2,1))

for (InitGenome in 1:10)
{
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'short-matu mam',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'short-matu mam',]$FrA, ylim = c(0,0.6), pch = '.', col = ColorA, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'short-matu mam',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'short-matu mam',]$FrT, ylim = c(0,0.6), pch = '.', col = ColorT, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'short-matu mam',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'short-matu mam',]$FrG, ylim = c(0,0.6), pch = '.', col = ColorG, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'short-matu mam',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'short-matu mam',]$FrC, ylim = c(0,0.6), pch = '.', col = ColorC, main = 'From long- to short- maturated species, neutral equilibrium, heavy chain', xlab = '', ylab = ''); par(new=TRUE)
}
legend(1,0.6, legend = c('Ah','Th','Gh','Ch'), col = c(ColorA, ColorT, ColorG, ColorC), pch = 16, cex=0.6)

for (InitGenome in 1:10)
{
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'long-matu mam',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'long-matu mam',]$FrA, ylim = c(0,0.65), pch = '.', col = ColorA, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'long-matu mam',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'long-matu mam',]$FrT, ylim = c(0,0.65), pch = '.', col = ColorT, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'long-matu mam',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'long-matu mam',]$FrG, ylim = c(0,0.65), pch = '.', col = ColorG, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'long-matu mam',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'long-matu mam',]$FrC, ylim = c(0,0.65), pch = '.', col = ColorC, main = 'From short- to long- maturated species, neutral equilibrium, heavy chain', xlab = '', ylab = ''); par(new=TRUE)
}
legend(1,0.6, legend = c('Ah','Th','Gh','Ch'), col = c(ColorA, ColorT, ColorG, ColorC), pch = 16, cex=0.6)


### aditional script to plot figure like 3A
to_check = Final[Final$Gener==1000000,]
rrA = aggregate(to_check$FrA, by=list(to_check$MutSpecProb), FUN= median)
rrC = aggregate(to_check$FrC, by=list(to_check$MutSpecProb), FUN= median)
rrG = aggregate(to_check$FrG, by=list(to_check$MutSpecProb), FUN= median)
rrT = aggregate(to_check$FrT, by=list(to_check$MutSpecProb), FUN= median)
colnames(rrA) = c('TypeMamm', 'FrAmedian')
colnames(rrC) = c('TypeMamm', 'FrCmedian')
colnames(rrG) = c('TypeMamm', 'FrGmedian')
colnames(rrT) = c('TypeMamm', 'FrTmedian')
fr_df_final = cbind(rrA, rrC$FrCmean)
fr_df_final = cbind(fr_df_final, rrG$FrGmean)
fr_df_final = cbind(fr_df_final, rrT$FrTmean)

ggplot(data = expvsobs, aes(x = observed, y = expected, group=type_of_mam, col = type_of_mam))+
  geom_point(size = 3.5)+
  geom_abline(col = 'gray3', linetype="longdash", size = 0.6)+
  theme_bw()+
  scale_color_manual(name="Type of Mammalia", labels = c('Short Maturated Mammalia', 'Long Maturated Mammalia'), values = c('short_mam'='deepskyblue4', 'long_mam' = 'firebrick3'))+
  geom_line(aes(group = mutation), col = 'black', size = 0.7)+
  geom_text(aes(label=mutation),hjust=-0.40, vjust=-0.40)+
  xlim(0, 0.60)+
  ylim(0,0.65)+
  labs(x = 'Observed Nucleotide Content',y = 'Expected Nucleotide Content')

