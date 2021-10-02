rm(list=ls(all=TRUE))

if (!require(ggpubr)) install.packages("ggpubr")
if (!require(vioplot)) install.packages("vioplot")

library("ggpubr")
library(vioplot)
###########Taxonomy###################################################################
Taxa = read.table("../../Body/1Raw/TaxaFromKostya.Names.stat", sep = '\t',header = FALSE) 
Taxa$Species = gsub(";.*",'',Taxa$V1); 
for (i in (1:nrow(Taxa)))  {Taxa$Species[i] = paste(unlist(strsplit(as.character(Taxa$Species[i]),split = ' '))[1],unlist(strsplit(as.character(Taxa$Species[i]),split = ' '))[2], sep = '_')}
Taxa$Class = gsub(";Chordata;.*",'',Taxa$V1); Taxa$Class = gsub(".*;",'',Taxa$Class); table(Taxa$Class)
Taxa$Class = gsub('Actinopteri','Actinopterygii',Taxa$Class)
Taxa$Class = gsub("Testudines|Squamata|Crocodylia|Sphenodontia",'Reptilia',Taxa$Class)
length(unique(Taxa$Species)) # 1708
table(Taxa$Class)
Taxa = Taxa[,-1]

TaxaMore = read.table("../../Body/1Raw/TaxaFromKostya.2NeedTaxa.tax.txt", sep = '\t',header = FALSE) 
TaxaMore$Species = ''
for (i in (1:nrow(TaxaMore)))  
{TaxaMore$Species[i] = paste(unlist(strsplit(as.character(TaxaMore$V1[i]),split = ' '))[1],unlist(strsplit(as.character(TaxaMore$V1[i]),split = ' '))[2], sep = '_')}
TaxaMore$Class = gsub("; Chordata;.*",'',TaxaMore$V2); 
TaxaMore$Class = gsub(".*; ",'',TaxaMore$Class); 
TaxaMore$Class = gsub('Actinopteri','Actinopterygii',TaxaMore$Class)
TaxaMore$Class = gsub("Testudines|Squamata|Crocodylia",'Reptilia',TaxaMore$Class)
table(TaxaMore$Class)
TaxaMore = TaxaMore[,-c(1,2)]

Taxa = rbind(Taxa,TaxaMore); Taxa = unique(Taxa)
#########################################################################################

MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)
MUTFROMK = merge(MUT, Taxa)
table(MUTFROMK$Class)
MUTFROMK = merge(MUT, Taxa)
MUTFROMK = MUTFROMK[MUTFROMK$T_C > 0,]
MUTFROMK = MUTFROMK[MUTFROMK$A_G > 0,]
MUTFROMK$TCdivAG = MUTFROMK$T_C / MUTFROMK$A_G

#pdf("../../Body/4Figures/VertebratePolymorphisms.MutSpecComparisons.Ecology.AtoG.DiffClasses.MutSpectrum.pdf")
#boxplot(MUTFROMK[MUTFROMK$Class == "Actinopterygii",]$T_C, MUTFROMK[MUTFROMK$Class == "Amphibia",]$T_C, MUTFROMK[MUTFROMK$Class == "Reptilia",]$T_C, MUTFROMK[MUTFROMK$Class == "Mammalia",]$T_C, MUTFROMK[MUTFROMK$Class == "Aves",]$T_C,
#        names = c(  "Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), ylab = "AH>GH", notch = TRUE)
#boxplot(MUTFROMK[MUTFROMK$Class == "Actinopterygii",]$TCdivAG, MUTFROMK[MUTFROMK$Class == "Amphibia",]$TCdivAG, MUTFROMK[MUTFROMK$Class == "Reptilia",]$TCdivAG, MUTFROMK[MUTFROMK$Class == "Mammalia",]$TCdivAG, MUTFROMK[MUTFROMK$Class == "Aves",]$TCdivAG, 
#        names = c(  "Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), ylab = "AH>GH div TH>CH", outline = FALSE, notch = TRUE)
#dev.off()
#HG = MUTFROMK[MUTFROMK$Species == "Heterocephalus_glaber",]
#MM = MUTFROMK[MUTFROMK$Species == "Mus_musculus",]
#CH = MUTFROMK[MUTFROMK$Species == "Cryptomys_hottentotus",]
#ALL = rbind(HG, CH)
#ALL = rbind(ALL, MM)
#ALL$TCdivAG = ALL$T_C / (ALL$A_G + ALL$T_C)

########AnAge data reading
AA = read.table("../../Body/1Raw/anage_data.txt", header = TRUE, sep = '\t')
AA$Species = paste(AA$Genus,AA$Species,sep = '_')
AA$TemperatureC = AA$Temperature..K. - 273.15
summary(AA$TemperatureC)
table(AA[AA$Class == "Aves",]$Temperature..K.)
######mergw with temperatures
TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)
class(TEMPE$Temperature)
class(TEMPE$Species)
TEMPE = aggregate(Temperature ~ ., median, data = TEMPE)
TEMPE$Class = "Actinopterygii"
Alltemp = data.frame(AA$Species, AA$Class, AA$TemperatureC); names(Alltemp) = c("Species", "Class", "Temperature")
Alltemp = rbind(Alltemp, TEMPE)

####U tests
colder <- c("Actinopterygii", "Amphibia", "Reptilia")
warmer <- c("Mammalia","Aves")
wilcox.test(MUTFROMK[MUTFROMK$Class %in% colder,]$T_C, MUTFROMK[MUTFROMK$Class %in% warmer,]$T_C, paired=F) #p-value = 2.574e-13
wilcox.test(MUTFROMK[MUTFROMK$Class %in% colder,]$A_G, MUTFROMK[MUTFROMK$Class %in% warmer,]$A_G, paired=F)  #p-value < 2.2e-16
wilcox.test(MUTFROMK[MUTFROMK$Class %in% colder,]$TCdivAG, MUTFROMK[MUTFROMK$Class %in% warmer,]$TCdivAG, paired=F) #p-value < 2.2e-16
wilcox.test(Alltemp[Alltemp$Class %in% colder,]$Temperature, Alltemp[Alltemp$Class %in% warmer,]$Temperature, paired=F) #p-value < 2.2e-16

pdf("../../Body/4Figures/5d.VertebratePolymorphisms.MutSpecComparisons.AtoG.DiffClasses.ViolinPlots.pdf", width = 10, height = 5.3)
MUTFROMK = MUTFROMK[MUTFROMK$T_C < 0.45,]
ggviolin(MUTFROMK, x = "Class", y = "T_C", select = c("Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), ylab = "AH>GH",
         order=c("Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), add = "boxplot", fill="Class", palette=c("#6760db", "#7849bf", "#9145c4", "#c73a69", "#c2464c"))
MUTFROMK = MUTFROMK[MUTFROMK$A_G < 0.15,]
ggviolin(MUTFROMK, x = "Class", y = "A_G", select = c("Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), ylab = "TH>CH",
         order=c("Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), add = "boxplot", fill="Class", palette=c("#6760db", "#7849bf", "#9145c4", "#c73a69", "#c2464c"))
MUTFROMK = MUTFROMK[MUTFROMK$TCdivAG < 20,]
ggviolin(MUTFROMK, x = "Class", y = "TCdivAG", select = c("Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), ylab = "AH>GH div TH>CH",
         order=c("Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), add = "boxplot", fill="Class", palette=c("#6760db", "#7849bf", "#9145c4", "#c73a69", "#c2464c"))
Alltemp=Alltemp[!is.na(Alltemp$Temperature),]
table(Alltemp[Alltemp$Class == "Aves",]$Temperature)
ggviolin(Alltemp, x = "Class", y = "Temperature", select = c("Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), ylab = "Body temperature, °C",
         order=c("Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), add = "boxplot", fill="Class", palette=c("#6760db", "#7849bf", "#9145c4", "#c73a69", "#c2464c"),  ylim = c(0, 50))
dev.off()

ac <- nrow(MUTFROMK[MUTFROMK$Class == "Actinopterygii"  & !is.na(MUTFROMK$T_C),]) #460
a <- nrow(MUTFROMK[MUTFROMK$Class == "Amphibia" & !is.na(MUTFROMK$T_C),]) #70
r <- nrow(MUTFROMK[MUTFROMK$Class == "Reptilia" & !is.na(MUTFROMK$T_C),]) #183
m <- nrow(MUTFROMK[MUTFROMK$Class == "Mammalia" & !is.na(MUTFROMK$T_C),]) #467
av <- nrow(MUTFROMK[MUTFROMK$Class == "Aves" & !is.na(MUTFROMK$T_C),]) #170
ourdataset <- ac + a + r + m + av #1350
nrow(MUTFROMK)

ac <- nrow(Alltemp[Alltemp$Class == "Actinopterygii" & !is.na(Alltemp$Temperature),]) #1590
a <- nrow(Alltemp[Alltemp$Class == "Amphibia" & !is.na(Alltemp$Temperature),]) #18
r <- nrow(Alltemp[Alltemp$Class == "Reptilia" & !is.na(Alltemp$Temperature),]) #16
m <- nrow(Alltemp[Alltemp$Class == "Mammalia" & !is.na(Alltemp$Temperature),]) #457
av <- nrow(Alltemp[Alltemp$Class == "Aves" & !is.na(Alltemp$Temperature),]) #3
ourdataset <- ac + a + r + m + av #2084
