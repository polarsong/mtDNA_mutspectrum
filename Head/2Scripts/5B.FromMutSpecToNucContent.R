rm(list=ls(all=TRUE))

Final = read.table("../../Body/2Derived/5A.FromMutSpecToNucContent.R.FinalTable.txt", header = TRUE)

###: PLOT
# real data for cold-water fishes: A ~ 0.24, C ~ 0.08, G ~ 0.34, T ~ 0.34 # total = 1
# real data for warm-water fishes: A ~ 0.20, C ~ 0.06, G ~ 0.36, T ~ 0.38 # total = 1
# summary(Final$InitGenome)
pdf("../../Body/4Figures/5B.FromMutSpecToNucContent.R.pdf",  width = 20, height = 10)
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


