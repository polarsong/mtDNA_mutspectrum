###################################
rm(list=ls(all=TRUE))
library(ggpubr)

setwd("../../Body/3Results")
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

# SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE)
names(SynNuc)
SynNuc = SynNuc[SynNuc$Gene != 'ND6',]

SynNuc = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species), FUN = sum)
names(SynNuc) = c('Species','NeutralA','NeutralT','NeutralG','NeutralC')
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 

### merge with temperature
TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)
TEMPE = aggregate(Temperature ~ ., median, data = TEMPE); summary(TEMPE$Temperature)
SynNuc = merge(TEMPE,SynNuc, by = 'Species', all = TRUE); summary(SynNuc$Temperature)

###### merge with fishbase maturation

MATUTM = read.table('../../Body/1Raw/FishBaseMaturity_Tm.txt',  header = TRUE)
MATUTM = aggregate(Tm ~ ., median, data = MATUTM) 
SynNuc = merge(MATUTM,SynNuc, by = 'Species', all = TRUE); nrow(SynNuc)

SynNuc$CtoTSkew = (SynNuc$FrC-SynNuc$FrT)/(SynNuc$FrC+SynNuc$FrT); summary(SynNuc$CtoTSkew) 
SynNuc$GtoASkew = (SynNuc$FrG-SynNuc$FrA)/(SynNuc$FrG+SynNuc$FrA); summary(SynNuc$GtoASkew)
SynNuc$AtoGSkew = (SynNuc$FrA-SynNuc$FrG)/(SynNuc$FrA+SynNuc$FrG); summary(SynNuc$AtoGSkew)
SynNuc$TtoCSkew = (SynNuc$FrT-SynNuc$FrC)/(SynNuc$FrT+SynNuc$FrC); summary(SynNuc$TtoCSkew)
SynNuc$TG = SynNuc$FrT+SynNuc$FrG
SynNuc$AC = SynNuc$FrA+SynNuc$FrC
SynNuc$TG_ACSkew = (SynNuc$TG-SynNuc$AC)/(SynNuc$TG+SynNuc$AC); summary(SynNuc$TG_ACSkew)
summary(SynNuc$TG_ACSkew *-1)


### ANALYSES:
summary(SynNuc$Temperature)
summary(SynNuc$Tm)
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
summary(lm(-TG_ACSkew ~ Temperature + Tm, data = SynNuc))

summary(lm(TG_ACSkew ~ Temperature * Tm, data = SynNuc))


#xlim = c(8, 14.5), ylim = c(-0.75, -0.25)


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


plot(temp$Temperature, temp$Tm, col="#42cbf5", xlab="Temperature", ylab="Time of maturetion, years", pch = 1, cex = (temp$TG_ACSkew*-7), ylim =c(0, 30))


pdf("../../Body/4Figures/WholeGenomeAnalyses.NucContent.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1D.pdf", width = 5.5, height = 5.5)
plot(temp[temp$Lifespan == "ShortMaturated",]$Temperature, temp[temp$Lifespan == "ShortMaturated",]$TG_ACSkew, col="#4da36c", xlab="Temperature", ylab="Fraction of AC_TCSkew", xlim = c(0, 30), ylim = c(-0.65, -0.1))
abline((-0.272616+-0.006709 ), -0.006184 , col="#4da36c", lwd = 2)
par(new=TRUE)
plot(temp[temp$Lifespan == "LongMaturated",]$Temperature, temp[temp$Lifespan == "LongMaturated",]$TG_ACSkew, col="#42cbf5", xlab="Temperature", ylab="Fraction of AC_TCSkew", xlim = c(0, 30), ylim = c(-0.65, -0.1))
abline(-0.272616, -0.006184 , col="#42cbf5", lwd = 2)
dev.off()


#### Rank corr
cor.test(log2(SynNuc$Temperature),SynNuc$FrA, method = "spearman")
cor.test(log2(SynNuc$Temperature),SynNuc$FrT, method = "spearman")
cor.test(log2(SynNuc$Temperature),SynNuc$FrG, method = "spearman")
cor.test(log2(SynNuc$Temperature),SynNuc$FrC, method = "spearman")

tempnumper = SynNuc[!is.na(SynNuc$Temperature),]
tempnumper = tempnumper[!is.na(tempnumper$FrT),]

samplesize = paste("N==", as.character(nrow(tempnumper)), sep="")

pdf("../../Body/4Figures/WholeGenomeAnalyses.NucContent.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1Ñ.pdf")
ggscatter(SynNuc, x = "Temperature", y = "FrT",
          color = "#e61a0b", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Mean annual water temperature, °C", ylab="Whole genome neutral fraction of A")+ stat_cor(
            aes(label = paste(..rr.label.., ..p.label.., samplesize, sep = "~`,`~")), 
            label.x = 3
          )

ggscatter(SynNuc, x = "Temperature", y = "FrG",
          color = "#009414", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Mean annual water temperature, °C", ylab="Whole genome neutral fraction of C")+ stat_cor(
            aes(label = paste(..rr.label.., ..p.label.., samplesize, sep = "~`,`~")), 
            label.x = 3
          )
ggscatter(SynNuc, x = "Temperature", y = "FrC",
          color = "#5c5c5c", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Mean annual water temperature, °C", ylab="Whole genome neutral fraction of G")+ stat_cor(
            aes(label = paste(..rr.label.., ..p.label.., samplesize, sep = "~`,`~")), 
            label.x = 3
          )
ggscatter(SynNuc, x = "Temperature", y = "FrA",
          color = "#0918e6", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Mean annual water temperature, °C", ylab="Whole genome neutral fraction of T")+ stat_cor(
            aes(label = paste(..rr.label.., ..p.label.., samplesize, sep = "~`,`~")), 
            label.x = 3
          )
dev.off()


pdf("../../Body/4Figures/WholeGenomeAnalyses.NucContent.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1E.pdf")
ggscatter(SynNuc, x = "Temperature", y = "TG_ACSkew",
          color = "#c99bc9", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Mean annual water temperature, °C", ylab="ACH-TGH skew")+ stat_cor(
            aes(label = paste(..rr.label.., ..p.label.., samplesize, sep = "~`,`~")), 
            label.x = 3
          )
dev.off()


#################################metabolic rate approximation
SynNuc$MR=(SynNuc$Tm+1)^0.75
SynNuc$TemperatureK = 273.15 + SynNuc$Temperature
SynNuc$B=SynNuc$MR * exp(-1.2/((8.617*10^-5)*SynNuc$TemperatureK))
cor.test(SynNuc$B, SynNuc$TtoCSkew, method="spearman") #rho  -0.3155981 
cor.test(SynNuc$Temperature, SynNuc$TtoCSkew, method = "spearman")

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
          ellipse = TRUE, mean.point = TRUE, add = "reg.line", xlab="Mean annual water temperature, °C", ylab="AC_TGSkew")
dev.off()

form=SynNuc[!is.na(SynNuc$Tm),]
form=form[!is.na(form$TG_ACSkew),]
form=form[!is.na(form$TemperatureK),]
form$TG_ACSkew=form$TG_ACSkew * (-1)

form$x = (1/log(form$Tm))*(log(form$TG_ACSkew)+0.62/((8.617*10^-5)*form$TemperatureK))
