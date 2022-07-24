
rm(list=ls(all=TRUE))

pdf("../../Body/4Figures/DuplexSeqDataAnalyses.Downsampling.Plotting.R.pdf", width = 12, height = 14)
par(mfrow=c(3,1))
breaks = seq(-10, 700,2)

LMSlopes = read.table("../../Body/3Results/DuplexSeqDataAnalyses.Downsampling.R.txt", header = TRUE)
names(LMSlopes)
LMSlopes = LMSlopes*10^9

summary(LMSlopes)

hist(LMSlopes$MusYoungAhGhDownSampledSlopes, xlim = c(0,600), ylim = c(0,10), col = 'red', border = 'red', breaks = breaks, main = '', xlab = '', ylab = '')
abline(v=median(LMSlopes$MusYoungAhGhDownSampledSlopes), col = 'dark red', lwd = 4, lt = 2); par(new=TRUE)
hist(LMSlopes$MusYoungChThDownSampledSlopes, xlim = c(0,600), ylim = c(0,10), col = 'dark grey', border = 'dark grey',  breaks = breaks, main = '', xlab = '', ylab = '')
abline(v=median(LMSlopes$MusYoungChThDownSampledSlopes), col = 'black', lwd = 4, lt = 2)

hist(LMSlopes$MusOldAhGhDownSampledSlopes, xlim = c(0,600), ylim = c(0,10), col = 'red', border = 'red', breaks = breaks, main = '', xlab = '', ylab = '')
abline(v=median(LMSlopes$MusOldAhGhDownSampledSlopes), col = 'dark red', lwd = 4, lt = 2); par(new=TRUE)
hist(LMSlopes$MusOldChThDownSampledSlopes, xlim = c(0,600), ylim = c(0,10), col = 'dark grey', border = 'dark grey', breaks = breaks, main = '', xlab = '', ylab = '')
abline(v=median(LMSlopes$MusOldChThDownSampledSlopes), col = 'black', lwd = 4, lt = 2)

hist(LMSlopes$HumanAhGhDownSampledSlopes, xlim = c(0,600), ylim = c(0,3), col = 'red', border = 'red', breaks = breaks, main = '', xlab = '', ylab = '')
abline(v=median(LMSlopes$HumanAhGhDownSampledSlopes), col = 'dark red', lwd = 4, lt = 2); par(new=TRUE)
hist(LMSlopes$HumanChThDownSampledSlopes, xlim = c(0,600), ylim = c(0,3), col = 'dark grey', border = 'dark grey', breaks = breaks, main = '', xlab = '', ylab = '')
abline(v=median(LMSlopes$HumanChThDownSampledSlopes), col = 'black', lwd = 4, lt = 2)

par(mfrow=c(2,2))

median(LMSlopes$MusYoungAhGhDownSampledSlopes) 
median(LMSlopes$MusYoungChThDownSampledSlopes)

median(LMSlopes$MusOldAhGhDownSampledSlopes)
median(LMSlopes$MusOldChThDownSampledSlopes)

median(LMSlopes$HumanAhGhDownSampledSlopes)
median(LMSlopes$HumanChThDownSampledSlopes)

# AhGh increase with age:
median(LMSlopes$MusOldAhGhDownSampledSlopes)/median(LMSlopes$MusYoungAhGhDownSampledSlopes) #  (8.3)
median(LMSlopes$HumanAhGhDownSampledSlopes)/median(LMSlopes$MusOldAhGhDownSampledSlopes)    #  (11.7)
median(LMSlopes$HumanAhGhDownSampledSlopes)/median(LMSlopes$MusYoungAhGhDownSampledSlopes)  #  (96.8)

# ChTh increase with age:
median(LMSlopes$MusOldChThDownSampledSlopes)/median(LMSlopes$MusYoungChThDownSampledSlopes) #  (1.9)
median(LMSlopes$HumanChThDownSampledSlopes)/median(LMSlopes$MusOldChThDownSampledSlopes)    #  (7.7)
median(LMSlopes$HumanChThDownSampledSlopes)/median(LMSlopes$MusYoungChThDownSampledSlopes)  #  (14.5)

######### how to compare datasets => how to prove that 8.3 > 1.9 etc..
#### 1: simple approach: divide vectors normally (each element one by one) and compare results
# 8.3 > 1.9
OldToYoungMiceAhGh = LMSlopes$MusOldAhGhDownSampledSlopes/LMSlopes$MusYoungAhGhDownSampledSlopes # 
OldToYoungMiceChTh = LMSlopes$MusOldChThDownSampledSlopes/LMSlopes$MusYoungChThDownSampledSlopes
boxplot(OldToYoungMiceAhGh,OldToYoungMiceChTh, notch = TRUE, outline = FALSE, col = c('red','grey'))
wilcox.test(OldToYoungMiceAhG,OldToYoungMiceChT)

# 11.7 > 7.7
HumanToOldMiceAhGh = LMSlopes$HumanAhGhDownSampledSlopes/LMSlopes$MusOldAhGhDownSampledSlopes  
HumanToOldMiceChTh = LMSlopes$HumanChThDownSampledSlopes/LMSlopes$MusOldChThDownSampledSlopes
boxplot(HumanToOldMiceAhGh,HumanToOldMiceChTh, notch = TRUE, outline = FALSE, col = c('red','grey'))
wilcox.test(HumanToOldMiceAhGh,HumanToOldMiceChTh)

# 96.8 > 14.5
HumanToYoungMiceAhGh = LMSlopes$HumanAhGhDownSampledSlopes/LMSlopes$MusYoungAhGhDownSampledSlopes  
HumanToYoungMiceChTh = LMSlopes$HumanChThDownSampledSlopes/LMSlopes$MusYoungChThDownSampledSlopes
boxplot(HumanToYoungMiceAhGh,HumanToYoungMiceChTh, notch = TRUE, outline = FALSE, col = c('red','grey'))
wilcox.test(HumanToYoungMiceAhGh,HumanToYoungMiceChTh)

#### 2: my approach: divide vectors in all possible combinations (each element of the first vector divide by all elements of the second) and compare results

# 8.3 > 1.9
OldToYoungMiceAhGh = data.frame(merge(LMSlopes$MusOldAhGhDownSampledSlopes,LMSlopes$MusYoungAhGhDownSampledSlopes))
names(OldToYoungMiceAhGh)=c('MusOldAhGhDownSampledSlopes','MusYoungAhGhDownSampledSlopes')
OldToYoungMiceAhGh$OldToYoungMiceAhGh = OldToYoungMiceAhGh$MusOldAhGhDownSampledSlopes/OldToYoungMiceAhGh$MusYoungAhGhDownSampledSlopes

OldToYoungMiceChTh = data.frame(merge(LMSlopes$MusOldChThDownSampledSlopes,LMSlopes$MusYoungChThDownSampledSlopes))
names(OldToYoungMiceChTh)=c('MusOldChThDownSampledSlopes','MusYoungChThDownSampledSlopes')
OldToYoungMiceChTh$OldToYoungMiceChTh = OldToYoungMiceChTh$MusOldChThDownSampledSlopes/OldToYoungMiceChTh$MusYoungChThDownSampledSlopes

boxplot(OldToYoungMiceAhGh$OldToYoungMiceAhGh,OldToYoungMiceChTh$OldToYoungMiceChTh, notch = TRUE, outline = FALSE, col = c('red','grey'))
wilcox.test(OldToYoungMiceAhGh$OldToYoungMiceAhGh,OldToYoungMiceChTh$OldToYoungMiceChTh)

# 11.7 > 7.7

# 96.8 > 14.5


dev.off()
