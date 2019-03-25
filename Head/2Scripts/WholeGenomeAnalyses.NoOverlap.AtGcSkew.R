################################
################################

rm(list=ls(all=TRUE))

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip", exdir = "../../Body/3Results/")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")

# SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE)
names(SynNuc)

### make ND6 complementary:
NotND6 = SynNuc[SynNuc$Gene != 'ND6',]
ND6 = SynNuc[SynNuc$Gene == 'ND6',]
A = ND6$NeutralT
T = ND6$NeutralA
G = ND6$NeutralC
C = ND6$NeutralG
ND6$NeutralA = A
ND6$NeutralT = T
ND6$NeutralG = G
ND6$NeutralC = C
SynNuc = rbind(NotND6,ND6)

VecOfTaxa = unique(SynNuc$Class)
table(SynNuc$Class)/13

AGG = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species,SynNuc$Class), FUN = sum)
names(AGG) = c('Species','Class','NeutralA','NeutralT','NeutralG','NeutralC')
AGG$ATSkew = (AGG$NeutralA - AGG$NeutralT)/(AGG$NeutralA + AGG$NeutralT)
AGG$GCSkew = (AGG$NeutralG - AGG$NeutralC)/(AGG$NeutralG + AGG$NeutralC)

pdf("../../Body/4Figures/WholeGenomeAnalyses.NoOverlap.AtGcSkew.R.01.pdf", height = 20, width = 40)
par(mfcol=c(2,3), cex = 2)
boxplot(AGG[AGG$Class == 'Actinopterygii',]$ATSkew,AGG[AGG$Class == 'Amphibia',]$ATSkew,AGG[AGG$Class == 'Reptilia',]$ATSkew,AGG[AGG$Class == 'Mammalia',]$ATSkew,AGG[AGG$Class == 'Aves',]$ATSkew, notch = TRUE, outline = FALSE, names = c('Actinopterygii','Amphibia','Reptilia','Mammalia','Aves'), main = 'AT skew', las = 1)
boxplot(AGG[AGG$Class == 'Actinopterygii',]$GCSkew,AGG[AGG$Class == 'Amphibia',]$GCSkew,AGG[AGG$Class == 'Reptilia',]$GCSkew,AGG[AGG$Class == 'Mammalia',]$GCSkew,AGG[AGG$Class == 'Aves',]$GCSkew, notch = TRUE, outline = FALSE, names = c('Actinopterygii','Amphibia','Reptilia','Mammalia','Aves'), main = 'GC skew', las = 1)

############# GENERATION LENGTH FOR MAMMALS
GT = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

Mam = merge(AGG,GT, by ='Species')
plot(log2(Mam$GenerationLength_d),Mam$ATSkew, main = 'Mammalia', xlab = 'log2(Generation Length)', ylab = 'AT skew', ylim=c(-0.12,0.7))
a<-lm(Mam$ATSkew ~ log2(Mam$GenerationLength_d)); summary(a)
abline(a, col = 'red', lwd = 4)
plot(log2(Mam$GenerationLength_d),Mam$GCSkew, main = 'Mammalia', xlab = 'log2(Generation Length)', ylab = 'GC skew', ylim=c(-1,-0.2))
a<-lm(Mam$GCSkew ~ log2(Mam$GenerationLength_d)); summary(a)
abline(a, col = 'red', lwd = 4)

cor.test(log2(Mam$GenerationLength_d),Mam$ATSkew, method = 'spearman') # 0.1183936, 0.0025, 
cor.test(log2(Mam$GenerationLength_d),Mam$GCSkew, method = 'spearman') # -0.05774663, 0.1414
median(log2(Mam$GenerationLength_d)) # 11.09672
cor.test(log2(Mam[log2(Mam$GenerationLength_d) > 11.09672,]$GenerationLength_d),Mam[log2(Mam$GenerationLength_d) > 11.09672,]$GCSkew, method = 'spearman') # -0.31, p = 8.205e-09
cor.test(log2(Mam[log2(Mam$GenerationLength_d) < 11.09672,]$GenerationLength_d),Mam[log2(Mam$GenerationLength_d) < 11.09672,]$GCSkew, method = 'spearman') #  0.2871864, p = 1.256e-07

#### PICS, ALINA
# cor.test(log2(Mam$GenerationLength_d),Mam$ATSkew, method = 'spearman')
# cor.test(log2(Mam$GenerationLength_d),Mam$GCSkew, method = 'spearman')

















############ AnAge
AA = read.table("../../Body/1Raw/anage_data.txt", header = TRUE, sep = '\t')
AA$Species = paste(AA$Genus,AA$Species,sep = '_')

### Mammalia - so, so
Mammalia = AGG[AGG$Class == 'Mammalia',]
Mammalia = merge(Mammalia,AA, by = 'Species')
### Maximum longevity years are two times better known than femaly maturity days
nrow(Mammalia[!is.na(Mammalia$Maximum.longevity..yrs.),]) # 387
nrow(Mammalia[!is.na(Mammalia$Female.maturity..days.),])  # 359
cor.test(Mammalia$Maximum.longevity..yrs.,Mammalia$ATSkew, method = 'spearman') # 0.02185936, 0.6682
cor.test(Mammalia$Maximum.longevity..yrs.,Mammalia$GCSkew, method = 'spearman') # -0.1683312, 0.0008858
median(Mammalia$Maximum.longevity..yrs.)
cor.test(Mammalia$Maximum.longevity..yrs.,Mammalia$GCSkew, method = 'spearman') # -0.1683312, 0.0008858

### FISHES - yes!!!
FISHES = AGG[AGG$Class == 'Actinopterygii',]
FISHES = merge(FISHES,AA, by = 'Species')
### Maximum longevity years are two times better known than femaly maturity days
nrow(FISHES[!is.na(FISHES$Maximum.longevity..yrs.),]) # 206
nrow(FISHES[!is.na(FISHES$Female.maturity..days.),])  # 91
cor.test(FISHES$Maximum.longevity..yrs.,FISHES$ATSkew, method = 'spearman') # 0.3206948, 2.607e-06
cor.test(FISHES$Maximum.longevity..yrs.,FISHES$GCSkew, method = 'spearman') # -0.1975776, 0.004418
plot(log2(FISHES$Maximum.longevity..yrs.),FISHES$ATSkew, main = 'Actinopterygii', xlab = 'log2(Maximal Lifespan)', ylab = 'AT skew', ylim=c(-0.12,0.7))
plot(log2(FISHES$Maximum.longevity..yrs.),FISHES$GCSkew, main = 'Actinopterygii', xlab = 'log2(Maximal Lifespan)', ylab = 'GC skew', ylim=c(-1,-0.2))

### Reptilia - nothing
Reptilia = AGG[AGG$Class == 'Reptilia',]
Reptilia = merge(Reptilia,AA, by = 'Species')
### Maximum longevity years are two times better known than femaly maturity days
nrow(Reptilia[!is.na(Reptilia$Maximum.longevity..yrs.),]) # 86
nrow(Reptilia[!is.na(Reptilia$Female.maturity..days.),])  # 18
cor.test(Reptilia$Maximum.longevity..yrs.,Reptilia$ATSkew, method = 'spearman') #
cor.test(Reptilia$Maximum.longevity..yrs.,Reptilia$GCSkew, method = 'spearman') #

### Amphibia - nothing
Amphibia = AGG[AGG$Class == 'Amphibia',]
Amphibia = merge(Amphibia,AA, by = 'Species')
### Maximum longevity years are two times better known than femaly maturity days
nrow(Amphibia[!is.na(Amphibia$Maximum.longevity..yrs.),]) # 24
nrow(Amphibia[!is.na(Amphibia$Female.maturity..days.),])  # 13
cor.test(Amphibia$Maximum.longevity..yrs.,Amphibia$ATSkew, method = 'spearman') #
cor.test(Amphibia$Maximum.longevity..yrs.,Amphibia$GCSkew, method = 'spearman') #

### Aves - nothing
Aves = AGG[AGG$Class == 'Aves',]
Aves = merge(Aves,AA, by = 'Species')
### Maximum longevity years are two times better known than femaly maturity days
nrow(Aves[!is.na(Aves$Maximum.longevity..yrs.),]) # 147
nrow(Aves[!is.na(Aves$Female.maturity..days.),])  # 111
cor.test(Aves$Maximum.longevity..yrs.,Aves$ATSkew, method = 'spearman') # 
cor.test(Aves$Maximum.longevity..yrs.,Aves$GCSkew, method = 'spearman') # 

dev.off()

