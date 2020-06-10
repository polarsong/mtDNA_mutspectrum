rm(list=ls(all=TRUE))
library("ggpubr")

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

###########################MUT spectrum in Actinoptery##################################
MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)
MUTACTINOPTERITAXAFROMK = merge(MUT, Taxa)
MUTACTINOPTERITAXAFROMK = MUTACTINOPTERITAXAFROMK[MUTACTINOPTERITAXAFROMK$Class == "Actinopterygii",]

MUTACTINOPTERITAXAFROMK=MUTACTINOPTERITAXAFROMK[,-1]
MUTACTINOPTERITAXAFROMK=MUTACTINOPTERITAXAFROMK[,-13]
temp = data.frame(t(MUTACTINOPTERITAXAFROMK))
a=c()

for(i in 1:ncol(temp)){
  b = data.frame(temp[i], rownames(temp))
  colnames(b)=c("Freq", "Sub")
  a = rbind(a, b)
}

pdf("../../Body/4Figures/VertebratePolymorphisms.MutSpecComparisons.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1A.pdf")
ggbarplot(a, x = "Sub", y = "Freq", fill = "Sub", color = "Sub",
          palette = c("#bdbdbd", "#876363", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#1b7c98", "#b83e3e", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2a836f", "#bdbdbd"), 
          xlab="Substitution types", ylab="Normalised frequencies", add = "mean_se")  
dev.off()
##########################################################################################


##########TEMPERATURE 
TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)
class(TEMPE$Temperature)
class(TEMPE$Species)
TEMPE = aggregate(Temperature ~ ., median, data = TEMPE)
TemperMut = merge(MUT, TEMPE) 

cor.test(TemperMut$A_T,TemperMut$Temperature, method = 'spearman')   #rho  
cor.test(TemperMut$A_G,TemperMut$Temperature, method = 'spearman')   #rho     -0.3581037 p-value = 3.321e-05
cor.test(TemperMut$A_C,TemperMut$Temperature, method = 'spearman')   #rho   
cor.test(TemperMut$T_A,TemperMut$Temperature, method = 'spearman')   #rho   
cor.test(TemperMut$T_G,TemperMut$Temperature, method = 'spearman')   #rho   
cor.test(TemperMut$T_C,TemperMut$Temperature, method = 'spearman')   #rho     0.2648037 p-value = 0.002522
cor.test(TemperMut$G_A,TemperMut$Temperature, method = 'spearman')   #rho   
cor.test(TemperMut$G_T,TemperMut$Temperature, method = 'spearman')   #rho  
cor.test(TemperMut$G_C,TemperMut$Temperature, method = 'spearman')   #rho  
cor.test(TemperMut$C_A,TemperMut$Temperature, method = 'spearman')   #rho   
cor.test(TemperMut$C_T,TemperMut$Temperature, method = 'spearman')   #rho   
cor.test(TemperMut$C_G,TemperMut$Temperature, method = 'spearman')   #rho   

pdf("../../Body/4Figures/VertebratePolymorphisms.MutSpecComparisons.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1B.pdf")
ggscatter(TemperMut, x = "Temperature", y = "T_C",
          color = "#2a836f", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n"), xlab="Temperature, C", ylab="AH>GH")

ggscatter(TemperMut, x = "Temperature", y = "A_G",
          color = "#876363", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n"), xlab="Temperature, C", ylab="TH>CH")
dev.off()

##########MATURITY Lm (mean length at first maturity in  )  and Tm (Mean or median age at first maturity)
MATULM = read.table('../../Body/1Raw/FishBaseMaturity_Lm.txt',  header = TRUE, stringsAsFactors = FALSE)
MATUTM = read.table('../../Body/1Raw/FishBaseMaturity_Tm.txt',  header = TRUE)
class(MATULM$Lm)
class(MATUTM$Tm)
MATULM$Lm = as.numeric(MATULM$Lm)
table(MATULM$Lm)
MATULM = MATULM[!is.na(MATULM$Lm),]
MATUTM = aggregate(Tm ~ ., median, data = MATUTM) 
MATULM = aggregate(Lm ~ ., median, data = MATULM)
MATULmmut = merge(MUT,MATULM) 
MATUTmmut = merge(MUT,MATUTM)

cor.test(MATULmmut$A_T,MATULmmut$Lm, method = 'spearman')   #rho  
cor.test(MATULmmut$A_G,MATULmmut$Lm, method = 'spearman')   #rho            
cor.test(MATULmmut$A_C,MATULmmut$Lm, method = 'spearman')   #rho   
cor.test(MATULmmut$T_A,MATULmmut$Lm, method = 'spearman')   #rho   
cor.test(MATULmmut$T_G,MATULmmut$Lm, method = 'spearman')   #rho   
cor.test(MATULmmut$T_C,MATULmmut$Lm, method = 'spearman')   #rho      
cor.test(MATULmmut$G_A,MATULmmut$Lm, method = 'spearman')   #rho   
cor.test(MATULmmut$G_T,MATULmmut$Lm, method = 'spearman')   #rho   
cor.test(MATULmmut$G_C,MATULmmut$Lm, method = 'spearman')   #rho  -0.266859         p-value = 0.008235
cor.test(MATULmmut$C_A,MATULmmut$Lm, method = 'spearman')   #rho     
cor.test(MATULmmut$C_T,MATULmmut$Lm, method = 'spearman')   #rho    
cor.test(MATULmmut$C_G,MATULmmut$Lm, method = 'spearman')   #rho   

cor.test(MATUTmmut$A_T,MATUTmmut$Tm, method = 'spearman')   #rho  
cor.test(MATUTmmut$A_G,MATUTmmut$Tm, method = 'spearman')   #rho            
cor.test(MATUTmmut$A_C,MATUTmmut$Tm, method = 'spearman')   #rho   
cor.test(MATUTmmut$T_A,MATUTmmut$Tm, method = 'spearman')   #rho   
cor.test(MATUTmmut$T_G,MATUTmmut$Tm, method = 'spearman')   #rho   
cor.test(MATUTmmut$T_C,MATUTmmut$Tm, method = 'spearman')   #rho       
cor.test(MATUTmmut$G_A,MATUTmmut$Tm, method = 'spearman')   #rho   
cor.test(MATUTmmut$G_T,MATUTmmut$Tm, method = 'spearman')   #rho   
cor.test(MATUTmmut$G_C,MATUTmmut$Tm, method = 'spearman')   #rho  -0.2219755     p-value = 0.0222
cor.test(MATUTmmut$C_A,MATUTmmut$Tm, method = 'spearman')   #rho     
cor.test(MATUTmmut$C_T,MATUTmmut$Tm, method = 'spearman')   #rho    
cor.test(MATUTmmut$C_G,MATUTmmut$Tm, method = 'spearman')   #rho   


############Multiple models
allparameters=TemperMut#128 species

ltest = lm(formula = Temperature ~ scale(T_C) + scale(A_G), data = allparameters)
summary(ltest)    

allparameters=merge(TemperMut, MATUTM)#65 species

ltest = lm(formula = Temperature ~ scale(T_C) + scale(A_G), data = allparameters)
summary(ltest)   

ltest = lm(formula = T_C ~ scale(Temperature) + scale(Tm), data = allparameters)
summary(ltest)

ltest = lm(formula = A_G ~ scale(Temperature) + scale(Tm), data = allparameters)
summary(ltest)

ltest = lm(formula = T_C ~ scale(Temperature), data = allparameters)
summary(ltest)

ltest = lm(formula = A_G ~ scale(Temperature), data = allparameters)
summary(ltest)


### TC divided by AG rank corr and log2 models
allparameters=allparameters[allparameters$A_G != 0,]
allparameters=allparameters[allparameters$T_C != 0,]
allparameters$TCdivAG=allparameters$T_C/allparameters$A_G
cor.test(allparameters$TCdivAG,allparameters$Temperature, method = 'spearman')  

ltest = lm(formula = log2(TCdivAG) ~ scale(Temperature) + scale(Tm), data = allparameters)
summary(ltest)

ltest = lm(formula = Temperature ~ scale(TCdivAG), data = allparameters)
summary(ltest) 

pdf("../../Body/4Figures/VertebratePolymorphisms.MutSpecComparisons.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1C.pdf")
ggscatter(allparameters, x = "Temperature", y = "TCdivAG",
          color = c("#814194"), # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n"), yscale = "log2", xlab="Temperature, C", ylab="log2 TCdivAG")
dev.off()

