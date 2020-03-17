########################
################ ATGC FRACTIONS ALONG GENOMES
########################
## plot segments of ATGC along genomes
## for ND6 we count complementary nucleotides (T instead of A and so on) so that we analyze only light chain content.
## find species with unsusual nucleotide content:

rm(list=ls(all=TRUE))

############ Syn mut

unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip", exdir = "../../Body/3Results/")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")

###### work with Mammals only:
table(SynNuc$Class)
SynNuc = SynNuc[SynNuc$Class == 'Mammalia',]

  
###### make a sequence of synonymous fourfold sies ranked according to the time being signle stranded

Result = c()
ResWholeMatrix = c()

VecOfSpecies = as.character(unique(SynNuc$Species))

for (sp in 1:length(VecOfSpecies))
{ # sp = 1
    TEMP = SynNuc[SynNuc$Species == VecOfSpecies[sp],] # the file is sorted according to the start of each gene
    
    ## concatenate all genes of a given species in correct order:
    Nucl = c(
    as.character(TEMP[TEMP$Gene == 'COX1',]$CodonsNoOverlap),
    as.character(TEMP[TEMP$Gene == 'COX2',]$CodonsNoOverlap),
    as.character(TEMP[TEMP$Gene == 'ATP8',]$CodonsNoOverlap), # Faith et al 2003 deleted it (too short)
    as.character(TEMP[TEMP$Gene == 'ATP6',]$CodonsNoOverlap),
    as.character(TEMP[TEMP$Gene == 'COX3',]$CodonsNoOverlap),
    as.character(TEMP[TEMP$Gene == 'ND3',]$CodonsNoOverlap),  # Faith et al 2003 deleted it (too short)
    as.character(TEMP[TEMP$Gene == 'ND4L',]$CodonsNoOverlap), # Faith et al 2003 deleted it (too short)
    as.character(TEMP[TEMP$Gene == 'ND4',]$CodonsNoOverlap),
    as.character(TEMP[TEMP$Gene == 'ND1',]$CodonsNoOverlap), # methods in Faith et al 2003 
    as.character(TEMP[TEMP$Gene == 'ND5',]$CodonsNoOverlap),
    as.character(TEMP[TEMP$Gene == 'ND2',]$CodonsNoOverlap), # methods in Faith et al 2003
    as.character(TEMP[TEMP$Gene == 'CYTB',]$CodonsNoOverlap))
    Nucl  = paste(Nucl,collapse = '')
    
    ## nucleotides => codons
    Nucl = unlist(strsplit(Nucl,''))
    StartNuc = 1
    Codons = c()
    # if (length(CodonsVec)/3 == integer
    for (i in 1:(length(Nucl)/3))
    {
      Codons = c(Codons,paste(Nucl[StartNuc : (StartNuc+2)],collapse = ''))
      StartNuc = StartNuc+3
    }
    
    ## codons => synon fourfold degener codons
    SynFourfoldDegenerCodons = c(
    'CTA','GTA','TCA','CCA','ACA','GCA','CGA','GGA',
    'CTT','GTT','TCT','CCT','ACT','GCT','CGT','GGT',
    'CTG','GTG','TCG','CCG','ACG','GCG','CGG','GGG',
    'CTC','GTC','TCC','CCC','ACC','GCC','CGC','GGC')
    CodonsNeutral  = data.frame(Codons[Codons %in% SynFourfoldDegenerCodons])
    names(CodonsNeutral) = c('codons')
    
    ## synon fourfond degener codons => third positions only
    LastNuc <- function(x)
    {
      unlist(strsplit(as.character(x),''))[3]  
    }
    CodonsNeutral$third = apply(as.matrix(CodonsNeutral$codons),1,FUN = LastNuc)
    
    ResLine = c()
    total = nrow(CodonsNeutral) # 1627
    steps = round(total/25 - 0.5) # # 1627/25
    for (perc in 1:steps)
    {# perc = 65
      Start = perc*25 - 24
      End = perc*25
      temp = CodonsNeutral[Start:End,] # CodonsNeutral,total - total*(perc/50))
      FrA = nrow(temp[temp$third == 'A',])/nrow(temp)
      FrT = nrow(temp[temp$third == 'T',])/nrow(temp)
      FrG = nrow(temp[temp$third == 'G',])/nrow(temp)
      FrC = nrow(temp[temp$third == 'C',])/nrow(temp)
      AtSkew = (nrow(temp[temp$third == 'A',]) - nrow(temp[temp$third == 'T',]))/(nrow(temp[temp$third == 'A',]) + nrow(temp[temp$third == 'T',]))
      GcSkew = (nrow(temp[temp$third == 'G',]) - nrow(temp[temp$third == 'C',]))/(nrow(temp[temp$third == 'G',]) + nrow(temp[temp$third == 'C',]))
      ResLine = rbind(ResLine,c(perc,FrA,FrT,FrG,FrC,AtSkew,GcSkew))
    }
      
    ResLine = as.data.frame(ResLine)
    names(ResLine)=c('perc','FrA','FrT','FrG','FrC','AtSkew','GcSkew')

    A=cor.test(ResLine$perc,ResLine$FrA, method = 'spearman')
    T=cor.test(ResLine$perc,ResLine$FrT, method = 'spearman')
    G=cor.test(ResLine$perc,ResLine$FrG, method = 'spearman')
    C=cor.test(ResLine$perc,ResLine$FrC, method = 'spearman')
    At=cor.test(ResLine$perc,ResLine$AtSkew, method = 'spearman')
    Gc=cor.test(ResLine$perc,ResLine$GcSkew, method = 'spearman')
    FrA.P = as.numeric(A[3]); FrA.Rho = as.numeric(A[4])
    FrT.P = as.numeric(T[3]); FrT.Rho = as.numeric(T[4])
    FrG.P = as.numeric(G[3]); FrG.Rho = as.numeric(G[4])
    FrC.P = as.numeric(C[3]); FrC.Rho = as.numeric(C[4])
    At.P = as.numeric(At[3]); At.Rho = as.numeric(At[4])
    Gc.P = as.numeric(Gc[3]); Gc.Rho = as.numeric(Gc[4])

    Result = rbind(Result,c(VecOfSpecies[sp],FrA.Rho,FrT.Rho,FrG.Rho,FrC.Rho,At.Rho,Gc.Rho,FrA.P,FrT.P,FrG.P,FrC.P,At.P,Gc.P,as.character(TEMP$Class[1]),as.character(TEMP$Taxonomy[1])))
    
    ResLine$Species = VecOfSpecies[sp]
    ResWholeMatrix = rbind(ResWholeMatrix,ResLine)
}

############### ANALYSIS OF "Result": 
Result = as.data.frame(Result)
names(Result) = c('Species','FrA.Rho','FrT.Rho','FrG.Rho','FrC.Rho','At.Rho','Gc.Rho','FrA.P','FrT.P','FrG.P','FrC.P','At.P','Gc.P','Class','Taxonomy')

Result$FrA.Rho = as.numeric(as.character(Result$FrA.Rho))
Result$FrT.Rho = as.numeric(as.character(Result$FrT.Rho))
Result$FrG.Rho = as.numeric(as.character(Result$FrG.Rho))
Result$FrC.Rho = as.numeric(as.character(Result$FrC.Rho))
Result$At.Rho = as.numeric(as.character(Result$At.Rho))
Result$Gc.Rho = as.numeric(as.character(Result$Gc.Rho))
Result$FrA.P = as.numeric(as.character(Result$FrA.P))
Result$FrT.P = as.numeric(as.character(Result$FrT.P))
Result$FrG.P = as.numeric(as.character(Result$FrG.P))
Result$FrC.P = as.numeric(as.character(Result$FrC.P))
Result$At.P = as.numeric(as.character(Result$At.P))
Result$Gc.P = as.numeric(as.character(Result$Gc.P))

ColG = rgb(0.1,0.1,0.1,0.5)
ColT = rgb(0.1,0.1,1,0.5)
ColC = rgb(0.1,1,0.1,0.5)
ColA = rgb(1,0.1,0.1,0.5)

pdf("../../Body/4Figures/WholeGenomeAnalyses.AtgcAlongGenomes.NoOverlapWindowsNucleotideFreq.Intercept&Slope.R01.pdf",height = 30, width = 30)

par(mfrow=c(2,3), cex = 3)
plot(Result$FrA.Rho,-log10(Result$FrA.P), ylim = c(0,10), col = ColA)
plot(Result$FrT.Rho,-log10(Result$FrT.P), ylim = c(0,10), col = ColT)
plot(Result$At.Rho,-log10(Result$At.P), ylim = c(0,10), col = 'gray')
plot(Result$FrG.Rho,-log10(Result$FrG.P), ylim = c(0,10), col = ColG)
plot(Result$FrC.Rho,-log10(Result$FrC.P), ylim = c(0,10), col = ColC)
plot(Result$Gc.Rho,-log10(Result$Gc.P), ylim = c(0,10), col = 'gray')

Breaks = seq(-0.8,0.8,0.02)
par(mfrow=c(2,3), cex = 3)
hist(Result$FrA.Rho, breaks = Breaks, xlim = c(-0.8,0.8), col = ColA, ylim = c(0,180), xlab = 'Rho', main = '', border = ColA); abline(v=0, col='black', lwd = 3); par(new = TRUE); 
hist(Result[Result$FrA.P < 0.01,]$FrA.Rho, breaks = Breaks, xlim = c(-0.8,0.8), col = ColA, ylim = c(0,180), xlab = 'Rho', main = 'A', border = 'black');
hist(Result$FrT.Rho, breaks = Breaks, xlim = c(-0.8,0.8), col = ColT, ylim = c(0,180), xlab = 'Rho', main = '', border = ColT); abline(v=0, col='black', lwd = 3); par(new = TRUE); 
hist(Result[Result$FrT.P < 0.01,]$FrT.Rho, breaks = Breaks, xlim = c(-0.8,0.8), col = ColT, ylim = c(0,180), xlab = 'Rho', main = 'T', border = 'black');
hist(Result$At.Rho, breaks = Breaks, xlim = c(-0.8,0.8), col = 'gray', ylim = c(0,180), xlab = 'Rho', main = '', border = 'gray'); abline(v=0, col='black', lwd = 3); par(new = TRUE); 
hist(Result[Result$At.P < 0.01,]$At.Rho, breaks = Breaks, xlim = c(-0.8,0.8), col = 'gray', ylim = c(0,180), xlab = 'Rho', main = 'AT skew', border = 'black');

hist(Result$FrG.Rho, breaks = Breaks, xlim = c(-0.8,0.8), col = ColG, ylim = c(0,180), xlab = 'Rho', main = '', border = ColG); abline(v=0, col='black', lwd = 3); par(new = TRUE); 
hist(Result[Result$FrG.P < 0.01,]$FrG.Rho, breaks = Breaks, xlim = c(-0.8,0.8), col = ColG, ylim = c(0,180), xlab = 'Rho', main = 'G', border = 'black');
hist(Result$FrC.Rho, breaks = Breaks, xlim = c(-0.8,0.8), col = ColC, ylim = c(0,180), xlab = 'Rho', main = '', border = ColC); abline(v=0, col='black', lwd = 3); par(new = TRUE); 
hist(Result[Result$FrC.P < 0.01,]$FrC.Rho, breaks = Breaks, xlim = c(-0.8,0.8), col = ColC, ylim = c(0,180), xlab = 'Rho', main = 'C', border = 'black');
hist(Result$Gc.Rho, breaks = Breaks, xlim = c(-0.8,0.8), col = 'gray', ylim = c(0,180), xlab = 'Rho', main = '', border = 'gray'); abline(v=0, col='black', lwd = 3); par(new = TRUE); 
hist(Result[Result$Gc.P < 0.01,]$Gc.Rho, breaks = Breaks, xlim = c(-0.8,0.8), col = 'gray', ylim = c(0,180), xlab = 'Rho', main = 'GCskew', border = 'black');


par(mfrow=c(5,4), cex = 1)
hist(Result[Result$Class == 'Mammalia',]$FrA.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColA); abline(v=0, col='red', lwd = 3)
hist(Result[Result$Class == 'Mammalia',]$FrT.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColT); abline(v=0, col='red', lwd = 3)
hist(Result[Result$Class == 'Mammalia',]$FrG.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColG); abline(v=0, col='red', lwd = 3)
hist(Result[Result$Class == 'Mammalia',]$FrC.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColC); abline(v=0, col='red', lwd = 3)

############# GENERATION LENGTH FOR ALL MAMMALS

GT = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

Mammals = merge(Result,GT, by = 'Species'); nrow(Mammals)

par(mfrow=c(2,2), cex = 3)
plot(log2(Mammals$GenerationLength_d), Mammals$FrA.Rho, col = ColA, ylim = c(-0.6,0.6)); abline(h = 0, lwd = 2, col = 'red', lt = 1)
plot(log2(Mammals$GenerationLength_d), Mammals$FrT.Rho, col = ColT, ylim = c(-0.6,0.6)); abline(h = 0, lwd = 2, col = 'red', lt = 1)
plot(log2(Mammals$GenerationLength_d), Mammals$FrG.Rho, col = ColG, ylim = c(-0.6,0.6)); abline(h = 0, lwd = 2, col = 'red', lt = 1)
plot(log2(Mammals$GenerationLength_d), Mammals$FrC.Rho, col = ColC, ylim = c(-0.6,0.6)); abline(h = 0, lwd = 2, col = 'red', lt = 1)

cor.test(Mammals$FrA.Rho,log2(Mammals$GenerationLength_d), method = 'spearman') #  negative - A is symmetrical, don't include it
cor.test(Mammals$FrT.Rho,log2(Mammals$GenerationLength_d), method = 'spearman') #  positive
cor.test(Mammals$FrG.Rho,log2(Mammals$GenerationLength_d), method = 'spearman') #  negative
cor.test(Mammals$FrC.Rho,log2(Mammals$GenerationLength_d), method = 'spearman') #  positive
A <- lm(log2(Mammals$GenerationLength_d) ~ Mammals$FrC.Rho + Mammals$FrG.Rho + Mammals$FrT.Rho)
summary(A)

### fisher test: signif. non signif, long - short lived
# mammals  - c
A = nrow(Mammals[Mammals$FrC.Rho > 0 & Mammals$FrC.P <0.01 & Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),])
B = nrow(Mammals[Mammals$FrC.Rho > 0 & Mammals$FrC.P <0.01 & Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),])
C = nrow(Mammals[Mammals$FrC.Rho > 0 & Mammals$FrC.P >= 0.01 & Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),])
D = nrow(Mammals[Mammals$FrC.Rho > 0 & Mammals$FrC.P >= 0.01 & Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),])
A/(A+B) # 57% of longlived mammals among significant
C/(C+D) # 45%% of longlived mammals among non significant 
A # 142
B # 105
C # 180
D # 215
# A/B = 1.35; C/D = 0.84; 1.35/0.84 =  1.6
X = cbind(c(A,B),c(C,D))
fisher.test(X) # 1.61, p = 0.003521
mosaicplot(X)

# mammals  - T
A = nrow(Mammals[Mammals$FrT.Rho < 0 & Mammals$FrT.P <0.01 & Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),])
B = nrow(Mammals[Mammals$FrT.Rho < 0 & Mammals$FrT.P <0.01 & Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),])
C = nrow(Mammals[Mammals$FrT.Rho < 0 & Mammals$FrT.P >= 0.01 & Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),])
D = nrow(Mammals[Mammals$FrT.Rho < 0 & Mammals$FrT.P >= 0.01 & Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),])
X = cbind(c(A,B),c(C,D))
fisher.test(X) # 0.86, p = 0.435

# mammals  - G
A = nrow(Mammals[Mammals$FrG.Rho < 0 & Mammals$FrG.P <0.01 & Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),])
B = nrow(Mammals[Mammals$FrG.Rho < 0 & Mammals$FrG.P <0.01 & Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),])
C = nrow(Mammals[Mammals$FrG.Rho < 0 & Mammals$FrG.P >= 0.01 & Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),])
D = nrow(Mammals[Mammals$FrG.Rho < 0 & Mammals$FrG.P >= 0.01 & Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),])
X = cbind(c(A,B),c(C,D))
fisher.test(X) # 0.1144

# mammals  - AtSkew
A = nrow(Mammals[Mammals$At.Rho > 0 & Mammals$At.P <0.01 & Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),])
B = nrow(Mammals[Mammals$At.Rho > 0 & Mammals$At.P <0.01 & Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),])
C = nrow(Mammals[Mammals$At.Rho > 0 & Mammals$At.P >= 0.01 & Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),])
D = nrow(Mammals[Mammals$At.Rho > 0 & Mammals$At.P >= 0.01 & Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),])
X = cbind(c(A,B),c(C,D))
fisher.test(X) # 0.313

# mammals  - GcSkew
A = nrow(Mammals[Mammals$Gc.Rho < 0 & Mammals$Gc.P <0.01 & Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),])
B = nrow(Mammals[Mammals$Gc.Rho < 0 & Mammals$Gc.P <0.01 & Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),])
C = nrow(Mammals[Mammals$Gc.Rho < 0 & Mammals$Gc.P >= 0.01 & Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),])
D = nrow(Mammals[Mammals$Gc.Rho < 0 & Mammals$Gc.P >= 0.01 & Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),])
X = cbind(c(A,B),c(C,D))
A # 92
B # 58
C # 218
D # 241
fisher.test(X) # odds = 1.751909; p = 0.003534
mosaicplot(X)

par(mfrow=c(5,4), cex = 1)
hist(Mammals$FrA.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColA); abline(v=0, col='red', lwd = 3)
hist(Mammals$FrG.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColG); abline(v=0, col='red', lwd = 3)
hist(Mammals$FrT.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColT); abline(v=0, col='red', lwd = 3)
hist(Mammals$FrC.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColC); abline(v=0, col='red', lwd = 3)

hist(Mammals[Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),]$FrA.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColA); abline(v=0, col='red', lwd = 3)
hist(Mammals[Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),]$FrG.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColG); abline(v=0, col='red', lwd = 3)
hist(Mammals[Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),]$FrT.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColT); abline(v=0, col='red', lwd = 3)
hist(Mammals[Mammals$GenerationLength_d > quantile(Mammals$GenerationLength_d,0.5),]$FrC.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColC); abline(v=0, col='red', lwd = 3)

hist(Mammals[Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),]$FrA.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColA); abline(v=0, col='red', lwd = 3)
hist(Mammals[Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),]$FrG.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColG); abline(v=0, col='red', lwd = 3)
hist(Mammals[Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),]$FrT.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColT); abline(v=0, col='red', lwd = 3)
hist(Mammals[Mammals$GenerationLength_d <= quantile(Mammals$GenerationLength_d,0.5),]$FrC.Rho, breaks = 50, xlim = c(-0.8,0.8), col = ColC); abline(v=0, col='red', lwd = 3)


################## ANALYSIS OF ResWholeMatrix
## nucleotide content as a function of GT and time being single stranded (perc)
ResWholeMatrix = merge(ResWholeMatrix,GT, by ='Species')

a<-lm(ResWholeMatrix$FrT ~ ResWholeMatrix$perc + log2(ResWholeMatrix$GenerationLength_d)); summary(a)
a<-lm(ResWholeMatrix$FrT ~ ResWholeMatrix$perc*log2(ResWholeMatrix$GenerationLength_d)); summary(a)
a<-lm(ResWholeMatrix$FrT ~ scale(ResWholeMatrix$perc) + scale(log2(ResWholeMatrix$GenerationLength_d))); summary(a)
a<-lm(ResWholeMatrix$FrT ~ scale(ResWholeMatrix$perc)*scale(log2(ResWholeMatrix$GenerationLength_d))); summary(a) # negative, negative, positive
# R^2 = 0.1274
#  (Intercept)                                                                0.1912182  0.0004448 429.879  < 2e-16 ***
#  scale(ResWholeMatrix$perc)                                                -0.0311871  0.0004450 -70.078  < 2e-16 ***
#  scale(log2(ResWholeMatrix$GenerationLength_d))                            -0.0146407  0.0004448 -32.912  < 2e-16 ***
#  scale(ResWholeMatrix$perc):scale(log2(ResWholeMatrix$GenerationLength_d))  0.0022713  0.0004468   5.083 3.73e-07 ***
# Residual standard error: 0.09114 on 42005 degrees of freedom
# Multiple R-squared:  0.1275,	Adjusted R-squared:  0.1274 
# F-statistic:  2045 on 3 and 42005 DF,  p-value: < 2.2e-16

a<-lm(ResWholeMatrix$FrC ~ ResWholeMatrix$perc + log2(ResWholeMatrix$GenerationLength_d)); summary(a)
a<-lm(ResWholeMatrix$FrC ~ ResWholeMatrix$perc*log2(ResWholeMatrix$GenerationLength_d)); summary(a)
a<-lm(ResWholeMatrix$FrC ~ scale(ResWholeMatrix$perc) + scale(log2(ResWholeMatrix$GenerationLength_d))); summary(a)
a<-lm(ResWholeMatrix$FrC ~ scale(ResWholeMatrix$perc)*scale(log2(ResWholeMatrix$GenerationLength_d))); summary(a) # positive, positive, positive 
# R^2 = 0.1168
#  (Intercept)                                                               0.2824121  0.0004751 594.367  < 2e-16 ***
#  scale(ResWholeMatrix$perc)                                                0.0270797  0.0004754  56.965  < 2e-16 ***
#  scale(log2(ResWholeMatrix$GenerationLength_d))                            0.0219669  0.0004752  46.229  < 2e-16 ***
#  scale(ResWholeMatrix$perc):scale(log2(ResWholeMatrix$GenerationLength_d)) 0.0018501  0.0004773   3.876 0.000106 ***

a<-lm(ResWholeMatrix$FrA ~ ResWholeMatrix$perc + log2(ResWholeMatrix$GenerationLength_d)); summary(a)
a<-lm(ResWholeMatrix$FrA ~ ResWholeMatrix$perc*log2(ResWholeMatrix$GenerationLength_d)); summary(a)
a<-lm(ResWholeMatrix$FrA ~ scale(ResWholeMatrix$perc) + scale(log2(ResWholeMatrix$GenerationLength_d))); summary(a) # perc increases A, GT decreases A
a<-lm(ResWholeMatrix$FrA ~ scale(ResWholeMatrix$perc)*scale(log2(ResWholeMatrix$GenerationLength_d))); summary(a)
# R^2 = 0.02223

a<-lm(ResWholeMatrix$FrG ~ ResWholeMatrix$perc + log2(ResWholeMatrix$GenerationLength_d)); summary(a)
a<-lm(ResWholeMatrix$FrG ~ ResWholeMatrix$perc*log2(ResWholeMatrix$GenerationLength_d)); summary(a)
a<-lm(ResWholeMatrix$FrG ~ scale(ResWholeMatrix$perc) + scale(log2(ResWholeMatrix$GenerationLength_d))); summary(a) # perc decreases G, GT increases G
a<-lm(ResWholeMatrix$FrG ~ scale(ResWholeMatrix$perc)*scale(log2(ResWholeMatrix$GenerationLength_d))); summary(a)
# R^2 = 0.04164

## Try to analyze fraction of T from T and C
ResWholeMatrix$FrT1 =  ResWholeMatrix$FrT / (ResWholeMatrix$FrT + ResWholeMatrix$FrC)
a<-lm(ResWholeMatrix$FrT1 ~ ResWholeMatrix$perc + log2(ResWholeMatrix$GenerationLength_d)); summary(a)
a<-lm(ResWholeMatrix$FrT1 ~ ResWholeMatrix$perc*log2(ResWholeMatrix$GenerationLength_d)); summary(a) # negative, negative, non significant
a<-lm(ResWholeMatrix$FrT1 ~ scale(ResWholeMatrix$perc) + scale(log2(ResWholeMatrix$GenerationLength_d))); summary(a)
a<-lm(ResWholeMatrix$FrT1 ~ scale(ResWholeMatrix$perc)*scale(log2(ResWholeMatrix$GenerationLength_d))); summary(a) # negative, negative, non significant effect
# (Intercept)                                                                0.4011864  0.0008016 500.461   <2e-16 ***
# scale(ResWholeMatrix$perc)                                                -0.0626613  0.0008020 -78.129   <2e-16 ***
# scale(log2(ResWholeMatrix$GenerationLength_d))                            -0.0357323  0.0008017 -44.572   <2e-16 ***
# scale(ResWholeMatrix$perc):scale(log2(ResWholeMatrix$GenerationLength_d))  0.0012187  0.0008053   1.513     0.13   
# Adjusted R-squared:  0.1647 

dev.off()

