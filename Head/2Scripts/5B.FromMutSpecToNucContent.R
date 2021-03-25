rm(list=ls(all=TRUE))

ColorA = rgb(1,0,0,0.3)
ColorT = rgb(0,1,0,0.3)
ColorG = rgb(0,0,1,0.3)
ColorC = rgb(0,1,1,0.3)

Final = read.table("../../Body/2Derived/5A.FromMutSpecToNucContent.R.FinalTable.txt", header = TRUE)

###: PLOT
# real data for cold-water fishes: A ~ 0.24, C ~ 0.08, G ~ 0.34, T ~ 0.34 # total = 1
# real data for warm-water fishes: A ~ 0.20, C ~ 0.06, G ~ 0.36, T ~ 0.38 # total = 1
# summary(Final$InitGenome)
pdf("../../Body/4Figures/5B.FromMutSpecToNucContent.R.pdf",  width = 20, height = 10) # dev.off()
par(mfrow=c(1,2))
for (InitGenome in 1:10)
{
plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$FrA, ylim = c(0,0.6), pch = '.', col = ColorA, main = '', xlab = '', ylab = ''); par(new=TRUE)
plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$FrT, ylim = c(0,0.6), pch = '.', col = ColorT, main = '', xlab = '', ylab = ''); par(new=TRUE)
plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$FrG, ylim = c(0,0.6), pch = '.', col = ColorG, main = '', xlab = '', ylab = ''); par(new=TRUE)
plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$FrC, ylim = c(0,0.6), pch = '.', col = ColorC, main = 'Cold-water, neutral equilibrium, heavy chain', xlab = '', ylab = ''); par(new=TRUE)
}
G = 0.2775772029; T= 0.4598399750; C=0.07858235961 # data from Valerian
A = 1 - G - C - T # 0.1840005
abline(h = A, col = ColorA, lwd = 3, lt = 2)
abline(h = T, col = ColorT, lwd = 3, lt = 2)
abline(h = G, col = ColorG, lwd = 3, lt = 2)
abline(h = C, col = ColorC, lwd = 3, lt = 2)
legend(1,0.6, legend = c('A','T','G','C'), col = c(ColorA, ColorT, ColorG, ColorC), pch = 16)

for (InitGenome in 1:10)
{
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$FrA, ylim = c(0,0.6), pch = '.', col = ColorA, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$FrT, ylim = c(0,0.6), pch = '.', col = ColorT, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$FrG, ylim = c(0,0.6), pch = '.', col = ColorG, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'warm-water fish',]$FrC, ylim = c(0,0.6), pch = '.', col = ColorC, main = 'Warm-water, neutral equilibrium, heavy chain', xlab = '', ylab = ''); par(new=TRUE)
}
A=0.1423421338; G=0.2897844548; T=0.5225696239; C=0.04530378870; # data from Valerian:
abline(h = A, col = ColorA, lwd = 3, lt = 2)
abline(h = T, col = ColorT, lwd = 3, lt = 2)
abline(h = G, col = ColorG, lwd = 3, lt = 2)
abline(h = C, col = ColorC, lwd = 3, lt = 2)
legend(1,0.6, legend = c('A','T','G','C'), col = c(ColorA, ColorT, ColorG, ColorC), pch = 16)

dev.off()

### make final boxplots and print out summary (mean, median)

### IF EQUILIBRIUM SENSITIVE TO STARTING CONDITIONS?
### IF EQUILIBRIUM IS SENSITIVE TO MUTSPEC (cold versus warm fishes)
### how the same derive with analogous functions (Lynch, Valerian, Stepan)


