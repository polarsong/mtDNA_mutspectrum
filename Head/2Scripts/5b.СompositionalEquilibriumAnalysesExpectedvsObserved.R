rm(list=ls(all=T))

ColorA = rgb(1,0,0,0.3)
ColorT = rgb(0,1,0,0.3)
ColorG = rgb(0,0,1,0.3)
ColorC = rgb(0,1,1,0.3)

Final = read.table("../../Body/2Derived/5A.FromMutSpecToNucContent.R.FinalTable.txt", header = TRUE)

Final = Final[Final$Gener <= 600000,]
###: PLOT

pdf("../../Body/4Figures/5B.FromMutSpecToNucContent.R.pdf",  width = 20, height = 10) # dev.off()
par(mfrow=c(2,2))

for (InitGenome in 1:10)
{
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$FrA, ylim = c(0,0.6), pch = '.', col = ColorA, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$FrT, ylim = c(0,0.6), pch = '.', col = ColorT, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$FrG, ylim = c(0,0.6), pch = '.', col = ColorG, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$Gener,Final[Final$InitGenome == InitGenome & Final$MutSpecProb == 'cold-water fish',]$FrC, ylim = c(0,0.6), pch = '.', col = ColorC, main = 'Cold-water, neutral equilibrium, heavy chain', xlab = '', ylab = ''); par(new=TRUE)
}

G = 0.2378578262; T= 0.5271436896; C=0.07049385836 # data from Valerian
A =  0.1645046270
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
A = 0.1518187869; G=0.3236389788; T = 0.4862288867; C = 0.03831334710 # data from Valerian:
abline(h = A, col = ColorA, lwd = 3, lt = 2)
abline(h = T, col = ColorT, lwd = 3, lt = 2)
abline(h = G, col = ColorG, lwd = 3, lt = 2)
abline(h = C, col = ColorC, lwd = 3, lt = 2)
legend(1,0.6, legend = c('A','T','G','C'), col = c(ColorA, ColorT, ColorG, ColorC), pch = 16)

dev.off()
