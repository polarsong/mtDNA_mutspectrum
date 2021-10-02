###################################
rm(list=ls(all=TRUE))
library(ggpubr)
library(caper)
library(geiger)


### reading whole genomes database
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

SynNuc = SynNuc[SynNuc$Gene != 'ND6',]

####### obtaining neutral nucleotide fractions in whole genomes
SynNuc = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species), FUN = sum)
names(SynNuc) = c('Species','NeutralA','NeutralT','NeutralG','NeutralC')
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 

### merge whole genomes with temperature
TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)
TEMPE = aggregate(Temperature ~ ., median, data = TEMPE); summary(TEMPE$Temperature)
SynNuc = merge(TEMPE,SynNuc, by = 'Species'); summary(SynNuc$Temperature)

##############Rank corr 
#Suppl. mat. 2.a
cor.test(log2(SynNuc$Temperature),SynNuc$FrA, method = "spearman")
cor.test(log2(SynNuc$Temperature),SynNuc$FrT, method = "spearman")
cor.test(log2(SynNuc$Temperature),SynNuc$FrG, method = "spearman")
cor.test(log2(SynNuc$Temperature),SynNuc$FrC, method = "spearman")


###### merge whole genomes and temperature with time of maturation
MATUTM = read.table('../../Body/1Raw/FishBaseMaturity_Tm.txt',  header = TRUE)
MATUTM = aggregate(Tm ~ ., median, data = MATUTM) 
SynNuc = merge(SynNuc, MATUTM, by = 'Species'); nrow(SynNuc)

SynNuc$CtoTSkew = (SynNuc$FrC-SynNuc$FrT)/(SynNuc$FrC+SynNuc$FrT); summary(SynNuc$CtoTSkew) 
SynNuc$GtoASkew = (SynNuc$FrG-SynNuc$FrA)/(SynNuc$FrG+SynNuc$FrA); summary(SynNuc$GtoASkew)
SynNuc$AtoGSkew = (SynNuc$FrA-SynNuc$FrG)/(SynNuc$FrA+SynNuc$FrG); summary(SynNuc$AtoGSkew)
SynNuc$TtoCSkew = (SynNuc$FrT-SynNuc$FrC)/(SynNuc$FrT+SynNuc$FrC); summary(SynNuc$TtoCSkew)

SynNuc$TG = SynNuc$FrT+SynNuc$FrG
SynNuc$AC = SynNuc$FrA+SynNuc$FrC
SynNuc$TG_ACSkew = (SynNuc$TG-SynNuc$AC)/(SynNuc$TG+SynNuc$AC); summary(SynNuc$TG_ACSkew)
SynNuc$AC_TGSkew = -(SynNuc$TG-SynNuc$AC)/(SynNuc$TG+SynNuc$AC); summary(SynNuc$AC_TGSkew) #OUR



### ANALYSES:
summary(SynNuc$Temperature)
summary(SynNuc$Tm)
#Suppl. Mat. 2d
summary(lm(FrT ~ scale(Temperature)+scale(Tm), data = SynNuc))
summary(lm(FrT ~ log2(Temperature + 2)*log2(Tm), data = SynNuc))  # keep it for presentation!!!
summary(lm(FrT ~ log2(Temperature + 2)+log2(Tm), data = SynNuc))
summary(lm(FrG ~ log2(Temperature + 2)+log2(Tm), data = SynNuc)) # strong
summary(lm(FrA ~ log2(Temperature + 2)+log2(Tm), data = SynNuc)) # strong

summary(lm(GtoASkew ~ log2(Temperature + 2)+log2(Tm), data = SynNuc)) # the highest R^2 = 0.17 
summary(lm(CtoTSkew ~ log2(Temperature + 2)+log2(Tm), data = SynNuc))
summary(lm(TtoCSkew ~ log2(Temperature + 2)+log2(Tm), data = SynNuc))
summary(lm(TG_ACSkew ~ log2(Temperature + 2)+log2(Tm), data = SynNuc))
summary(lm(TG_ACSkew ~ scale(Temperature + 2)+scale(Tm), data = SynNuc))
#Suppl. Mat. 2b
summary(lm(AC_TGSkew ~ log2(Temperature + 2)+log2(Tm), data = SynNuc))
summary(lm(AC_TGSkew ~ scale(Temperature + 2)+scale(Tm), data = SynNuc))# ###PICS







###################################################### phylogenetic inertia analysis

tree = read.tree('../1Raw/mtalign.aln.treefile.rooted')

row.names(SynNuc) = SynNuc$Species

tree_pruned = treedata(tree, SynNuc, sort=T, warnings=T)$phy 
# 5104 sp in SynNuc, 3608 sp in tree_pruned

data<-as.data.frame(treedata(tree_pruned, SynNuc, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)

data$AC_TGSkew = as.numeric(as.character(data$AC_TGSkew))
data$Temperature = as.numeric(as.character(data$Temperature))
data$Tm = as.numeric(as.character(data$Tm))

data_comp <- comparative.data(tree_pruned, data[, c('Species', 'AC_TGSkew',
                                                    'Temperature', 'Tm')], Species, vcv=TRUE)

model = pgls(AC_TGSkew ~ scale(Temperature+2) + scale(Tm), data_comp, lambda="ML")
summary(model)

# lambda [ ML]  : 0.992
# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)            0.4182931  0.1214322  3.4447 0.0007736 ***
#   scale(Temperature + 2) 0.0125032  0.0063876  1.9574 0.0524728 .  
# scale(Tm)              0.0050481  0.0051593  0.9784 0.3296996    

model2 = pgls(AC_TGSkew ~ scale(Temperature+2), data_comp, lambda="ML")
summary(model2)

# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)            0.4231533  0.1223940  3.4573 0.0007395 ***
#   scale(Temperature + 2) 0.0101526  0.0060137  1.6882 0.0937819 .  


model3 = pgls(AC_TGSkew ~ log2(Temperature + 2) + log2(Tm), data_comp, lambda="ML")
summary(model3)

# lambda [ ML]  : 0.990
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)           0.3185648  0.1248595  2.5514  0.01191 *
#   log2(Temperature + 2) 0.0186083  0.0086042  2.1627  0.03242 *
#   log2(Tm)              0.0099184  0.0048863  2.0298  0.04445 *

model4 = pgls(AC_TGSkew ~ Temperature, data_comp, lambda="ML")
summary(model4)



#xlim = c(8, 14.5), ylim = c(-0.75, -0.25)


















######### plotting scatter with temperature and two groups
medianTm = median(SynNuc[!is.na(SynNuc$Tm),]$Tm)
SynNuc$Lifespan = "Na"

for (i in 1:nrow(SynNuc)){
  if (SynNuc$Tm[i] < medianTm & !is.na(SynNuc$Tm[i])){
    SynNuc$Lifespan[i] = "ShortMaturated"
  }
  if (SynNuc$Tm[i] > medianTm & !is.na(SynNuc$Tm[i])){
    SynNuc$Lifespan[i] = "LongMaturated"
  }
}
temp = SynNuc[!SynNuc$Lifespan == "Na",]
summary(lm(-TG_ACSkew ~ Temperature + Lifespan, data = temp))

plot(temp$Temperature, temp$Tm, col="#42cbf5", xlab="Temperature", ylab="Time of maturation, years", pch = 1, cex = (temp$TG_ACSkew*-7), ylim =c(0, 30))


pdf("../../Body/4Figures/WholeGenomeAnalyses.NucContent.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1D.pdf", width = 7, height = 8.5)
plot(temp[temp$Lifespan == "ShortMaturated",]$Temperature, temp[temp$Lifespan == "ShortMaturated",]$TG_ACSkew*-1, col="#4da36c", xlab="Mean annual water temperature, C?", ylab="STG-SAC skew", ylim=c(0.1, 0.65), xlim=c(0,30), pch = 16)
abline((0.331911-0.049196), 0.006172, col="#4da36c", lwd = 2)
par(new=TRUE)
plot(temp[temp$Lifespan == "LongMaturated",]$Temperature, temp[temp$Lifespan == "LongMaturated",]$TG_ACSkew*-1, col="#42cbf5", xlab="", ylab="", ylim=c(0.1, 0.65), xlim=c(0,30), pch = 16)
abline(0.331911, 0.006172, col="#42cbf5", lwd = 2)
legend("bottomright", legend=c( "Short time of maturation","Long time of maturation"), col=c("#4da36c","#42cbf5"), pch = 16)
dev.off()

pdf("../../Body/4Figures/WholeGenomeAnalyses.NucContent.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1CC.pdf", width = 7, height = 8.5)
plot(SynNuc$Temperature, SynNuc$AC_TGSkew, ylim=c(0.1, 0.65), xlim=c(0,30), pch = 16, col="#c99bc9", xlab="Mean annual water temperature, C?", ylab="STG-SAC skew" )
abline(lm(Temperature ~ AC_TGSkew, data = SynNuc))

#### Rank corr
cor.test(log2(SynNuc$Temperature),SynNuc$FrA, method = "spearman")
cor.test(log2(SynNuc$Temperature),SynNuc$FrT, method = "spearman")
cor.test(log2(SynNuc$Temperature),SynNuc$FrG, method = "spearman")
cor.test(log2(SynNuc$Temperature),SynNuc$FrC, method = "spearman")

tempnumper = SynNuc[!is.na(SynNuc$Temperature),]
tempnumper = tempnumper[!is.na(tempnumper$FrT),]

samplesize = paste("N==", as.character(nrow(tempnumper)), sep="")

pdf("../../Body/4Figures/WholeGenomeAnalyses.NucContent.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1?.pdf")
ggscatter(SynNuc, x = "Temperature", y = "FrT",
          color = "#e61a0b", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Mean annual water temperature, ?C", ylab="Whole genome neutral fraction of A")+ stat_cor(
            aes(label = paste(..rr.label.., ..p.label.., samplesize, sep = "~`,`~")), 
            label.x = 3
          )

ggscatter(SynNuc, x = "Temperature", y = "FrG",
          color = "#009414", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Mean annual water temperature, ?C", ylab="Whole genome neutral fraction of C")+ stat_cor(
            aes(label = paste(..rr.label.., ..p.label.., samplesize, sep = "~`,`~")), 
            label.x = 3
          )
ggscatter(SynNuc, x = "Temperature", y = "FrC",
          color = "#5c5c5c", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Mean annual water temperature, ?C", ylab="Whole genome neutral fraction of G")+ stat_cor(
            aes(label = paste(..rr.label.., ..p.label.., samplesize, sep = "~`,`~")), 
            label.x = 3
          )
ggscatter(SynNuc, x = "Temperature", y = "FrA",
          color = "#0918e6", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Mean annual water temperature, ?C", ylab="Whole genome neutral fraction of T")+ stat_cor(
            aes(label = paste(..rr.label.., ..p.label.., samplesize, sep = "~`,`~")), 
            label.x = 3
          )
dev.off()


pdf("../../Body/4Figures/WholeGenomeAnalyses.NucContent.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1E.pdf", width = 7, height = 8.5)
ggscatter(SynNuc, x = "Temperature", y = "AC_TGSkew",
          color = "#c99bc9", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Mean annual water temperature, ?C", ylab="STG-SAC skew", ylim=c(0.1, 0.65), xlim=c(0,30))+stat_cor(aes(label = paste(..rr.label.., ..p.label.., samplesize, sep = "~`,`~")))
dev.off()


ggscatter(SynNuc, x = "Temperature", y = "TtoCSkew",
          color = "#009414", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Temperature, C", ylab="AGskew")+ stat_cor(
            aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
            label.x = 3
          )
ggscatter(SynNuc, x = "Temperature", y = "AtoGSkew",
          color = "#009414", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Temperature, C", ylab="TCskew")+ stat_cor(
            aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
            label.x = 3
          )

ggscatter(SynNuc, x = "Temperature", y = "GtoASkew",
          color = "#009414", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Temperature, C", ylab="GAskew")+ stat_cor(
            aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
            label.x = 3
          )

pdf("../../Body/4Figures/WholeGenomeAnalyses.NucContent.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1D.pdf")
temp = SynNuc[!SynNuc$Lifespan == "Na",]
ggscatter(temp, x = "Temperature", y = "TG_ACSkew",
          color = "Lifespan", shape = "Lifespan",
          palette = c("#4da36c", "#c9a157"),
          ellipse = TRUE, mean.point = TRUE, add = "reg.line", xlab="Mean annual water temperature, ?C", ylab="AC_TGSkew")
dev.off()

#################################metabolic rate approximation
SynNuc$MR=(SynNuc$Tm+1)^0.75
SynNuc$TemperatureK = 273.15 + SynNuc$Temperature
SynNuc$B=SynNuc$MR * exp(-1.2/((8.617*10^-5)*SynNuc$TemperatureK))
cor.test(SynNuc$B, SynNuc$TtoCSkew, method="spearman") #rho  -0.3155981 
cor.test(SynNuc$Temperature, SynNuc$TtoCSkew, method = "spearman")


form=SynNuc[!is.na(SynNuc$Tm),]
form=form[!is.na(form$TG_ACSkew),]
form=form[!is.na(form$TemperatureK),]
form$TG_ACSkew=form$TG_ACSkew * (-1)

form$x = (1/log(form$Tm))*(log(form$TG_ACSkew)+0.62/((8.617*10^-5)*form$TemperatureK))





####################################
####################################
###Full ecology for Kuptsov
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




FISHsystematix2 = merge(TEMPE, SynNuc, all = T)
FISHsystematix2= merge(FISHsystematix2, MATULM, all = T)
FISHsystematix2= merge(FISHsystematix2, MATUTM, all = T)
FISHsystematix2$TemperatureC=FISHsystematix2$Temperature
FISHsystematix2 = FISHsystematix2[,-2]
str(FISHsystematix2)

FISHsystematix2= merge(FISHsystematix2, Taxa, by="Species", all = T)

table(FISHsystematix2$Class)
my = c("Ceratodontiformes", "Chondrichthyes", "Actinopterygii", "Cladistia")
FISHsystematix2 = FISHsystematix2[FISHsystematix2$Class %in% my,]
FISHsystematix2 = FISHsystematix2[!is.na(FISHsystematix2$FrA),]
FISHsystematix2 = FISHsystematix2[order(FISHsystematix2$FrA),]

table(FISHsystematix2$Class)
write.csv(FISHsystematix2, file = '../../Body/3Results/WholeGenomeAnalyses.AllActinopterygiiSystematix.csv', quote = FALSE)

