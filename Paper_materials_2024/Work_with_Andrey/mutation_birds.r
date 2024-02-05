rm(list = ls(all=TRUE))
BMR <- read.csv2("GlobalBMRbase.csv", header=T)
DNA <- read.csv2("Birds_dataset_paper.csv", header=T)
DNA <- DNA[c('Species', 'ghahSkew', 'chthSkew')]
DNA <- subset(DNA, Species!="Nothoprocta perdicaria") 

#BMR <- BMR[BMR$Order == "Passeriformes",]
#BMR <- BMR[BMR$Order != "Passeriformes",]
#BMR <- BMR[BMR$Trait == "BMR",]
#BMR <- BMR[BMR$Bad1Good2 == 2,]
#BMR <- BMR[BMR$Climate == "Temperate",]
#BMR <- BMR[BMR$Climate == "Tropical",]
#BMR <- BMR[BMR$Climate == "Sub",]
#BMR <- BMR[BMR$Method == "flow",]

BMR <- BMR[c('Species', 'Mass', 'TraitValue')]


DNA <- aggregate (DNA[,-1], by=list(DNA$Species), mean)
BMR <- aggregate (BMR[,-1], by=list(BMR$Species), mean)
BMR$log_Mass <- log10(BMR$Mass)
BMR$log_BMR <- log10(BMR$TraitValue)
BMR <- merge(BMR, DNA, by="Group.1", all.x=T)
BMR <- na.omit(BMR)

colnames(BMR)[1] <- "Species"
#write.csv2(BMR, file="MeanBMR.csv")

model_1 <- lm(ghahSkew ~ log_Mass + log_BMR, data=BMR)
model_2 <- lm(log10(ghahSkew) ~ log_Mass + log_BMR, data=BMR)

model_3 <- lm(chthSkew ~ log_Mass + log_BMR, data=BMR)
model_4 <- lm(log10(chthSkew) ~ log_Mass + log_BMR, data=BMR)

library(bbmle)
AICctab (model_1, model_2, weights=T) 
AICctab (model_3, model_4, weights=T) 

model <- model_1
#model <- model_4

#################################################################### Diagnostics

library(olsrr) # To perform Residual Diagnostics

ols_plot_resid_qq(model) # Graph for detecting violation of normality assumption
ols_test_normality(model) # Test for detecting violation of normality assumption
ols_test_correlation(model) # Correlation between observed residuals and expected residuals under normality
ols_plot_resid_hist(model) # Histogram of residuals for detecting violation of normality assumption
ols_plot_resid_fit(model) # Residual vs Fitted Values Plot (scatter plot of residuals on the y axis and fitted values on the x axis to detect non-linearity, unequal error variances, and outliers)
# Characteristics of a well behaved residual vs fitted plot:
# 1) The residuals spread randomly around the 0 line indicating that the relationship is linear.
# 2) The residuals form an approximate horizontal band around the 0 line indicating homogeneity of error variance.
# 3) No one residual is visibly away from the random pattern of the residuals indicating that there are no outliers.

plot(model)

library(nortest)
library(see)
library (performance)
check_normality(model) # same as shapiro.test(rstandard(model)) for LM
plot(check_normality(model)) # plot residuals and Gaussian distribution

shapiro.test(residuals(model))
shapiro.test(rstandard(model))
shapiro.test(rstudent(model))

library(gvlma)
summary(gvlma(model)) #Tests the main assumptions of regression analysis

plot(gvlma(model))

# Evaluate homoscedasticity
library(car)
ncvTest(model) # non-constant error variance test (Breusch–Pagan test). Null-hypothesis (homoskedasticity) is rejected when p < 0.05 (i.e. presence of Heteroscedasticity)
spreadLevelPlot(model) # plot studentized residuals vs. fitted values

#####################################################


summary (model_1)


plot(residuals(lm(log_BMR~log_Mass, data=BMR)) ~ log10(chthSkew), data=BMR)
abline(lm(residuals(lm(log_BMR~log_Mass, data=BMR)) ~ log10(chthSkew), data=BMR))

###########################################################
library(ape); library(phytools);  library(geiger)
feathertree <- read.tree("anc_kg.tre")
feathertree$node.label <- NULL # Remove internal node labels (if any)
is.ultrametric(feathertree)
is.binary(feathertree)
is.rooted(feathertree)

feathertree <- root(feathertree, outgroup = "Struthio_camelus", resolve.root = TRUE)

mycalibration <- makeChronosCalib(feathertree, node="root", age.max=94)
feathertree <- chronos(feathertree, lambda=0, model = 'correlated', calibration = mycalibration)

is.ultrametric(feathertree)
is.binary(feathertree)
is.rooted(feathertree)

writeNexus(feathertree, file="Ultrametric_feathertree.nex")

BMR$Species <- gsub(" ", "_", fixed = TRUE, BMR$Species) # Replaces spaces with _
listBMR <- BMR$Species
listTree <- feathertree$tip.label
#PGLS  for binar metrics
df_flight = read.csv('../../Paper_materials_2024/flight_and_gene.csv')
df_flight = df_flight[,c(2,3,4,5,6)]
df_flight$flight_num = 0
df_flight[df_flight$ability_to_fly != 'Flying',]$flight_num = 1
df_flight = df_flight[df_flight$ability_to_fly !='Sphenisciformes',]
rownames(df_flight) = df_flight$species_name
listSkew = df_flight$species_name
SpeciesToDrop <- setdiff(listTree, listSkew)
drop.tip(feathertree, SpeciesToDrop) -> Fly_skew_tree
rownames(df_flight) <- df_flight[,1] 
df_flight <- df_flight[match(Fly_skew_tree$tip.label,rownames(df_flight)),]
attach(df_flight)
names(GhAhSkew) = rownames(df_flight)
names(flight_num) = rownames(df_flight)
name.check(Fly_skew_tree, df_flight)
library(nlme)
m1 <- gls(GhAhSkew~flight_num, data=df_flight, correlation=corPagel(value = 1, Fly_skew_tree, form = ~species_name, fixed=FALSE), method="ML") # ML estimation
summary(m1)
library(caper)
CompBMR <- comparative.data(Fly_skew_tree, df_flight, 'species_name', na.omit=FALSE, vcv=TRUE, vcv.dim=3) #vcv.dim=2 is default
m2 <- pgls(log_BMR~log_Mass, data=CompBMR, lambda='ML') #The branch length transformations can be optimised between bounds using maximum likelihood by setting the value for a transformation to 'ML'
summary(m2)

model_1phy <- pgls(ghahSkew ~ log_Mass + log_BMR, data=CompBMR, lambda='ML')
model_2phy <- pgls(log10(ghahSkew) ~ log_Mass + log_BMR, data=CompBMR, lambda='ML')

AICctab (model_1phy, model_2phy, weights=T)

summary(model_1phy)









SpeciesToDrop <- setdiff(listTree, listBMR) #Вычисляет, какие виды убрать из общего дерева
drop.tip(feathertree, SpeciesToDrop) -> BMRTree #Убирает лишние виды из дерева
plotTree(BMRTree, fsize=0.5)

rownames(BMR) <- BMR[,1] 
BMR <- BMR[match(BMRTree$tip.label,rownames(BMR)),] # rearranging the rows of the data frame to match the order of species names in the tree object.
attach(BMR)
names(log_Mass) <- rownames(BMR)
names(log_BMR) <- rownames(BMR)
names(chthSkew) <- rownames(BMR)
names(ghahSkew) <- rownames(BMR)
name.check(BMRTree, BMR)
BMRTree
feathertree
library(nlme)
m1 <- gls(log_BMR~log_Mass, data=BMR, correlation=corPagel(value = 1, BMRTree, form = ~Species, fixed=FALSE), method="ML") # ML estimation
summary(m1)

library(caper)
CompBMR <- comparative.data(BMRTree, BMR, 'Species', na.omit=FALSE, vcv=TRUE, vcv.dim=3) #vcv.dim=2 is default
str(BMRTree)
m2 <- pgls(log_BMR~log_Mass, data=CompBMR, lambda='ML') #The branch length transformations can be optimised between bounds using maximum likelihood by setting the value for a transformation to 'ML'
summary(m2)

model_1phy <- pgls(ghahSkew ~ log_Mass + log_BMR, data=CompBMR, lambda='ML')
model_2phy <- pgls(log10(ghahSkew) ~ log_Mass + log_BMR, data=CompBMR, lambda='ML')

AICctab (model_1phy, model_2phy, weights=T)

summary(model_1phy)

#############################################################################################################################3



#Phylogenetic signal in Skews using full tree
graphics.off()
rm(list=ls())
library(ape); library(geiger)
DNA1 <- read.csv2("Birds_dataset_paper.csv", header=T)
DNA1 <- DNA1[c('Species', 'ghahSkew', 'chthSkew')]
DNA1 <- aggregate (DNA1[,-1], by=list(DNA1$Species), mean)
colnames(DNA1)[1] <- "Species"
DNA1$Species <- gsub(" ", "_", fixed = TRUE, DNA1$Species) # Replaces spaces with _

FullTree <- read.nexus("Ultrametric_feathertree.nex")
rownames(DNA1) <- DNA1[,1] 
DNA1 <- DNA1[match(FullTree$tip.label, rownames(DNA1)),] # rearranging the rows of the data frame to match the order of species names in the tree object.
attach(DNA1)
names(chthSkew) <- rownames(DNA1)
names(ghahSkew) <- rownames(DNA1)
name.check(FullTree, DNA1)

# 1) Estimating Pagel's lambda - 1st method
library(caper)
CompData <- comparative.data(FullTree, DNA1, 'Species', na.omit=FALSE, vcv=TRUE) #vcv.dim=2 is default)
model_1phy <- pgls(ghahSkew ~ 1, data=CompData, lambda='ML')
summary(model_1phy)

# 2) Estimating Pagel's lambda - 2nd method
library(nlme)
model_2phy <- gls(ghahSkew ~ 1, data=DNA1, correlation=corPagel(value = 1, FullTree, form = ~Species, fixed=FALSE), method="ML")
summary(model_2phy)

# 3) Estimating Pagel's lambda - 3rd method
library(phytools)
phylosig(FullTree, ghahSkew, method="lambda", test=T) #difference of Pagel's lambda from 0

# Testing difference of Pagel's lambda from 1 "by hands"
mFullTree <- FullTree
mFullTree$mapped.edge <- matrix(FullTree$edge.length, nrow(FullTree$edge), 1)
colnames(mFullTree$mapped.edge) <- "1"
resSkew <- brownie.lite(mFullTree, ghahSkew)
resLambda <- phylosig(FullTree, ghahSkew, method="lambda")
LR <- 2*(resLambda$logL - resSkew$logL1)
pchisq(LR, df=1, lower.tail=F)

# 4) Estimating Pagel's lambda - 4th method
library(phylolm)
model_3phy <- phylolm(ghahSkew ~ 1, data=DNA1, phy=FullTree, model="lambda")
summary(model_3phy)

# 5) Estimating Pagel's lambda - 5th method
library(geiger)
fitContinuous(FullTree, ghahSkew, model="lambda")

