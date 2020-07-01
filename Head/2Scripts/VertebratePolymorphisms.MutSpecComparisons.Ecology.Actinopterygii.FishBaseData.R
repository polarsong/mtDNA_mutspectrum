rm(list=ls(all=TRUE))

if (!require(caper)) install.packages("caper")
if (!require(geiger)) install.packages("geiger")

library(caper)
library(geiger)
library("ggpubr")

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
          color = "#036a5b", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n"), xlab="Temperature, C", ylab="AH>GH")

ggscatter(TemperMut, x = "Temperature", y = "A_G",
          color = "#73514f", # Points color, shape and size
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

summary(lm(formula = Temperature ~ scale(T_C) + scale(A_G), data = allparameters))

allparameters=merge(TemperMut, MATUTM)#65 species

summary(lm(formula = Temperature ~ scale(T_C) + scale(A_G), data = allparameters))

summary(lm(formula = T_C ~ scale(Temperature) + scale(Tm), data = allparameters))

summary(lm(formula = A_G ~ scale(Temperature) + scale(Tm), data = allparameters))

summary(lm(formula = T_C ~ scale(Temperature), data = allparameters))

summary(lm(formula = A_G ~ scale(Temperature), data = allparameters))


### TC divided by AG rank corr and log2 models
allparameters=TemperMut
allparameters=allparameters[allparameters$A_G != 0,]
allparameters=allparameters[allparameters$T_C != 0,]
allparameters$TCdivAG=allparameters$T_C/allparameters$A_G
#cor.test(allparameters$TCdivAG,allparameters$Temperature, method = 'spearman')  

summary(lm(formula = log2(TCdivAG) ~ scale(Temperature), data = allparameters))
summary(lm(formula = Temperature ~ scale(TCdivAG), data = allparameters))

allparameters=merge(TemperMut, MATUTM)
allparameters=allparameters[allparameters$A_G != 0,]
allparameters=allparameters[allparameters$T_C != 0,]
allparameters$TCdivAG=allparameters$T_C/allparameters$A_G
summary(lm(formula = Temperature ~ scale(TCdivAG), data = allparameters))
summary(lm(formula = log2(TCdivAG) ~ scale(Temperature) + scale(Tm), data = allparameters))

pdf("../../Body/4Figures/VertebratePolymorphisms.MutSpecComparisons.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1C.pdf")
ggscatter(allparameters, x = "Temperature", y = "TCdivAG",
          color = c("#814194"), # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n"), yscale = "log2", xlab="Temperature, C", ylab="log2 TCdivAG")
dev.off()

#####################################################################################
### PICs

tree = read.tree('../../Body/1Raw/mtalign.aln.treefile.rooted')

row.names(allparameters) = allparameters$Species

tree_pruned = treedata(tree, allparameters, sort=T, warnings=T)$phy

#   Not found in the tree and were dropped from the dataframe:
# Boops_boops
# Helicolenus_dactylopterus
# Lepidorhombus_whiffiagonis
# Lethrinus_olivaceus
# Lithognathus_mormyrus
# Lutjanus_vitta
# Pagellus_acarne
# Plagioscion_squamosissimus
# Sarda_sarda
# Scomberomorus_commerson
# Solea_solea
# Squalius_cephalus
# Squalius_pyrenaicus
# Trisopterus_minutus

data<-as.data.frame(treedata(tree_pruned, allparameters, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)

data$TCdivAG = as.numeric(as.character(data$TCdivAG))
data$Temperature = as.numeric(as.character(data$Temperature))
data$Tm = as.numeric(as.character(data$Tm))

data_comp <- comparative.data(tree_pruned, data, Species, vcv=TRUE)

model = pgls(TCdivAG ~ scale(Temperature) + scale(Tm), data_comp, lambda="ML")
summary(model)

# lambda [ ML]  : 0.731
# lower bound : 0.000, p = 0.03528
# upper bound : 1.000, p = 1.2987e-08
# 95.0% CI   : (0.217, 0.900)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)         5.76508    1.93865  2.9738 0.004715 **
#   scale(Temperature)  0.55015    0.53135  1.0354 0.306023   
# scale(Tm)          -0.14488    0.57032 -0.2540 0.800630   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6.819 on 45 degrees of freedom
# Multiple R-squared: 0.03013,	Adjusted R-squared: -0.01297 
# F-statistic: 0.6991 on 2 and 45 DF,  p-value: 0.5024 

model = pgls(log2(TCdivAG) ~ scale(Temperature) + scale(Tm), data_comp, lambda="ML")
summary(model)

# lambda [ ML]  : 0.535
# lower bound : 0.000, p = 0.14424
# upper bound : 1.000, p = 3.3975e-12
# 95.0% CI   : (NA, 0.830)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)         1.69761    0.63241  2.6843  0.01014 *
#   scale(Temperature)  0.47952    0.21438  2.2368  0.03030 *
#   scale(Tm)          -0.08112    0.22464 -0.3611  0.71971  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.436 on 45 degrees of freedom
# Multiple R-squared: 0.1189,	Adjusted R-squared: 0.07972 
# F-statistic: 3.036 on 2 and 45 DF,  p-value: 0.05798 

summary(lm(pic(data$TCdivAG, tree_pruned) ~ pic(data$Temperature, tree_pruned) +
             pic(data$Tm, tree_pruned)))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                        -0.961202   2.371449  -0.405    0.687
# pic(data$Temperature, tree_pruned) -0.005479   0.091820  -0.060    0.953
# pic(data$Tm, tree_pruned)           0.110200   0.127524   0.864    0.392
