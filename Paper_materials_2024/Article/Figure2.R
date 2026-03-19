rm(list = ls(all=TRUE))
library(ggplot2)
df_mtdna = read.csv('../Work_with_Andrey/Birds_dataset_paper.csv', header = TRUE, sep = ';')
df_mtdna$Mass = gsub(',', '.', df_mtdna$Mass)
df_mtdna$ghahSkew = gsub(',', '.', df_mtdna$ghahSkew)
df_mtdna$chthSkew = gsub(',', '.', df_mtdna$chthSkew)
df_mtdna$Mass = as.numeric(as.character(df_mtdna$Mass))
df_mtdna$ghahSkew = as.numeric(as.character(df_mtdna$ghahSkew))
df_mtdna$chthSkew = as.numeric(as.character(df_mtdna$chthSkew))
names_v = unique(df_mtdna$Species)
df_short = data.frame()
for (i in names_v)
{
  df1 = df_mtdna[df_mtdna$Species == i,]
  a = sum(df1$ghahSkew)/12
  b = sum(df1$chthSkew)/12
  v = sum(df1$Mass)/12
  ab = c(i, a, b, v)
  df_short = rbind(df_short, ab)
}
names(df_short) = c('Species', 'GhAhSkew', 'ThChSkew', 'Mass')
df_short$Mass = as.numeric(df_short$Mass)
df_short$GhAhSkew = as.numeric(df_short$GhAhSkew)
df_short$ThChSkew = as.numeric(df_short$ThChSkew)
df_short$log_mass = log10(df_short$Mass)

ggplot(df_short, aes(x = log_mass, y = GhAhSkew))+
  geom_point()+
  annotate('text', x = 4, y = 0.1, label = 'N = 766')+
  xlab('Decimal logarithm of mass')+
  geom_smooth(method = lm)

#Clutch
df_par = read.csv('../Work_with_Andrey/Species_life-histories.csv')
df_mtdna_par = merge(df_short, df_par, by = 'Species')
df_clutch = df_mtdna_par[,c(1,2,3,15)]
df_clutch = na.omit(df_clutch)
df_clutch$logclutch = log10(df_clutch$Clutch)

ggplot(df_clutch, aes(x = Clutch, y = GhAhSkew))+
  geom_point()+
  annotate('text', x = 4, y = 0.1, label = 'N = 203')+
  geom_smooth(method = lm)

ggplot(df_clutch, aes(x = logclutch, y = GhAhSkew))+
  geom_point()+
  annotate('text', x = 1, y = 0.1, label = 'N = 203')+
  geom_smooth(method = lm)

#BMR
df_bmr = read.csv('../Work_with_Andrey/GlobalBMRbase.csv', sep = ';')
df_bmr_e = df_bmr[df_bmr$Trait == 'BMR',]
names_v = unique(df_bmr_e$Species)
df_short_1 = data.frame()
df_bmr_e$TraitValue = gsub(',', '.', df_bmr_e$TraitValue)
df_bmr_e$TraitValue = suppressWarnings(as.numeric(df_bmr_e$TraitValue))
for (i in names_v)
{
  df1 = df_bmr_e[df_bmr_e$Species == i,]
  a1 = sum(df1$TraitValue)
  a2 = nrow(df1)
  a = a1/a2
  b = 'BMR'
  ab = c(i, b, a)
  df_short_1 = rbind(df_short_1, ab)
}
names(df_short_1) = c('Species', 'Trait', 'BMR_value')
df_mtdna_bmr = merge(df_short, df_short_1)
df_mtdna_bmr$BMR_value = as.numeric(df_mtdna_bmr$BMR_value)
df_mtdna_bmr$logBMR = log10(df_mtdna_bmr$BMR_value)
ggplot(df_mtdna_bmr, aes(x = BMR_value, y = GhAhSkew))+
  geom_point()+
  annotate('text', x = 2500, y = 0.1, label = 'N = 186')+
  geom_smooth(method = lm)
ggplot(df_mtdna_bmr, aes(x = logBMR, y = GhAhSkew))+
  geom_point()+
  annotate('text', x = 2, y = 0.1, label = 'N = 186')+
  geom_smooth(method = lm)

#longevity
df_long = read.csv('../Work_with_Andrey/AVES_longevity.csv')
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
df_long$scinam = firstup(df_long$scinam)
names(df_long) = c('Species', 'Longevity', 'Origin', 'Data')
df_long_correct = data.frame()
long_birds = unique(df_long$Species)
for (i in long_birds)
{
  bird = df_long[df_long$Species == i,]
  a = sum(bird$Longevity)/nrow(bird)
  df_long_correct = rbind(df_long_correct, c(i,a))
}
names(df_long_correct) = c('Species', 'Longevity')
df_short$Species = gsub(' ','_', df_short$Species)
df_long_correct$Species = gsub(' ','_',df_long_correct$Species)
df_long_mtdna = merge(df_long_correct, df_short)
df_long_mtdna$Longevity = as.numeric(as.character(df_long_mtdna$Longevity))
ggplot(df_long_mtdna, aes(x = Longevity, y = GhAhSkew))+
  geom_point()+
  annotate('text', x = 10, y = 0.1, label = 'N = 264')+
  geom_smooth(method = lm)

#normal stats
#mass
masslog_lm= lm(GhAhSkew ~ log_mass, data = df_short)
summary(masslog_lm)
mass_lm= lm(GhAhSkew ~ Mass, data = df_short)
summary(mass_lm)
#longevity
long_lm = lm(GhAhSkew ~ Longevity, data = df_long_mtdna)
summary(long_lm)
#BMR
bmr_lm = lm(GhAhSkew ~ logBMR, data = df_mtdna_bmr)
summary(bmr_lm)
#Clutch
clutch_lm = lm(GhAhSkew ~ Clutch, data = df_clutch)
summary(clutch_lm)

#pgls
library(ape) 
library(phytools) 
library(geiger)
library(nlme)
feathertree <- read.nexus("../Work_with_Andrey/Ultrametric_feathertree.nex")
feathertree$node.label <- NULL # Remove internal node labels (if any)
is.ultrametric(feathertree)
is.binary(feathertree)
is.rooted(feathertree)
#mass
row.names(df_short) = df_short$Species
name.check(feathertree, df_short)
df_short[df_short$Species == "Agapornis_pullarius" | df_short$Species == "Mergus_squamatus" | df_short$Species == "Vestiaria_coccinea",] = NA
df_short = na.omit(df_short)
name.check(feathertree, df_short)
spp_mass = rownames(df_short)
corLambda_mass = corPagel(value = 1, phy = feathertree, form=~spp_mass)
pgls_mass = gls(GhAhSkew~log_mass,
                data=df_short, correlation=corLambda_mass)
summary(pgls_mass)

mass = setNames(df_short[,"log_mass"], rownames(df_short))
Gh_mass = setNames(df_short[,"GhAhSkew"], rownames(df_short))
masspic = pic(mass, feathertree)
Gh_massPIC = pic(Gh_mass, feathertree)
fit_pic_mass = lm(masspic~Gh_massPIC+0)
fit_pic_mass
summary(fit_pic_mass)
plot(Gh_massPIC~masspic)

fit_pic_mass1 = lm(Gh_massPIC~masspic+0)
fit_pic_mass1
summary(fit_pic_mass1)
plot(masspic~Gh_massPIC)
#longevity
row.names(df_long_mtdna) = df_long_mtdna$Species
name.check(feathertree, df_long_mtdna)
df_long_mtdna[df_long_mtdna$Species == "Agapornis_pullarius" | df_long_mtdna$Species == "Vestiaria_coccinea",] = NA
df_long_mtdna = na.omit(df_long_mtdna)
listSkew_long = df_long_mtdna$Species
listTree_long <- feathertree$tip.label
SpeciesToDrop <- setdiff(listTree_long, listSkew_long)
drop.tip(feathertree, SpeciesToDrop) -> long_tree
spp_long = rownames(df_long_mtdna)
corLambda_long = corPagel(value = 1, phy = long_tree, form=~spp_long)
pgls_long = gls(GhAhSkew~Longevity,
                data=df_long_mtdna, correlation=corLambda_long)
summary(pgls_long)

long = setNames(df_long_mtdna[,"Longevity"], rownames(df_long_mtdna))
Gh_long = setNames(df_long_mtdna[,"GhAhSkew"], rownames(df_long_mtdna))
longpic = pic(long, long_tree)
Gh_longPIC = pic(Gh_long, long_tree)
fit_pic_long = lm(longpic~Gh_longPIC+0)
fit_pic_long
summary(fit_pic_long)
plot(Gh_longPIC~longpic)

#bmr
df_mtdna_bmr$Species = gsub(' ', '_', df_mtdna_bmr$Species)
row.names(df_mtdna_bmr) = df_mtdna_bmr$Species
name.check(feathertree, df_mtdna_bmr)
listSkew_bmr = df_mtdna_bmr$Species
listTree_bmr <- feathertree$tip.label
SpeciesToDrop_bmr <- setdiff(listTree_bmr, listSkew_bmr)
drop.tip(feathertree, SpeciesToDrop_bmr) -> bmr_tree
spp_bmr = rownames(df_mtdna_bmr)
corLambda_bmr = corPagel(value = 1, phy = bmr_tree, form=~spp_bmr)
pgls_bmr = gls(GhAhSkew~logBMR,
                data=df_mtdna_bmr, correlation=corLambda_bmr)
summary(pgls_bmr)

bmr = setNames(df_mtdna_bmr[,"logBMR"], rownames(df_mtdna_bmr))
Gh_bmr = setNames(df_mtdna_bmr[,"GhAhSkew"], rownames(df_mtdna_bmr))
bmrpic = pic(bmr, bmr_tree)
Gh_bmrPIC = pic(Gh_bmr, bmr_tree)
fit_pic_bmr = lm(bmrpic~Gh_bmrPIC+0)
fit_pic_bmr
summary(fit_pic_bmr)
plot(Gh_bmrPIC~bmrpic)

#clutch

df_clutch$Species = gsub(' ', '_', df_clutch$Species)
row.names(df_clutch) = df_clutch$Species
name.check(feathertree, df_clutch)
listSkew_clutch = df_clutch$Species
listTree_clutch <- feathertree$tip.label
SpeciesToDrop_clutch <- setdiff(listTree_clutch, listSkew_clutch)
drop.tip(feathertree, SpeciesToDrop_clutch) -> clutch_tree
spp_clutch = rownames(df_clutch)
corLambda_clutch = corPagel(value = 1, phy = clutch_tree, form=~spp_clutch)
pgls_clutch = gls(GhAhSkew~Clutch,
               data=df_clutch, correlation=corLambda_clutch)
summary(pgls_clutch)

clutch = setNames(df_clutch[,"Clutch"], rownames(df_clutch))
Gh_clutch = setNames(df_clutch[,"GhAhSkew"], rownames(df_clutch))
clutchpic = pic(clutch, clutch_tree)
Gh_clutchPIC = pic(Gh_clutch, clutch_tree)
fit_pic_clutch = lm(clutchpic~Gh_clutchPIC+0)
fit_pic_clutch
summary(fit_pic_clutch)
plot(Gh_clutchPIC~clutchpic)
fit_pic_clutch_r = lm(Gh_clutchPIC~clutchpic+0)
fit_pic_clutch_r
summary(fit_pic_clutch_r)
plot(clutchpic~Gh_clutchPIC)


