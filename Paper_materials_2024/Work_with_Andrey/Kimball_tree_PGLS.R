rm(list = ls(all=TRUE))
library(ape)
library(geiger)
df_kim = read.tree('Kimball.tre')
df_kim1 = df_kim[[1]]
a = df_kim1$tip.label
a = gsub("[^_]*-(.*)", "\\1",a)
a = gsub("[^_]*_(.*)", "\\1",a)
a = gsub("[^_]*_(.*)", "\\1",a)
df_kim1$tip.label = a
df_kim1$tip.label
df_cytb = read.csv('Midori2_birds_cytb_ghahskew_better.csv')
df_cytb = df_cytb[,c(2,3,4)]
df_cytb$species_name = gsub(' ', '_', df_cytb$species_name)
avo_data = read.csv('../../Body/1Raw/Avonet_data.csv')
avo_data1 = avo_data[,c('Species3', 'Migration')]
names(avo_data1) = c('species_name', 'migration')
avo_data1 = na.omit(avo_data1)
avo_data1$migration_try1 = 0
avo_data1[avo_data1$migration == 2 | avo_data1$migration == 3,]$migration_try1 = 1
avo_data1$migration_try2 = 0
avo_data1[avo_data1$migration == 3,]$migration_try2 = 1
avo_data1$species_name = gsub(' ', '_', avo_data1$species_name)
df_try = merge(df_cytb, avo_data1)

listSkew = df_try$species_name
listTree <- df_kim1$tip.label
SpeciesToDrop <- setdiff(listTree, listSkew)
drop.tip(df_kim1, SpeciesToDrop) -> Avo_cytb_tree
rownames(df_try) <- df_try[,1] 
df_try <- df_try[match(Avo_cytb_tree$tip.label,rownames(df_try)),]
attach(df_try)
names(GhAhSkew_start_codon) = rownames(df_try)
names(GhAhSkew_seq_begining) = rownames(df_try)
names(migration_try1) = rownames(df_try)
names(migration_try2) = rownames(df_try)
name.check(Avo_cytb_tree, df_try)
library(nlme)
m <- gls(GhAhSkew_start_codon~migration_try1, data=df_try, correlation=corPagel(value = 1, Avo_cytb_tree, form = ~species_name, fixed=FALSE), method="ML") # ML estimation
#summary(m)
library(caper)
CompBMR <- comparative.data(Avo_cytb_tree, df_try, 'species_name', na.omit=FALSE, vcv=TRUE, vcv.dim=3) #vcv.dim=2 is default
m_extra <- pgls(GhAhSkew_seq_begining~migration_try1, data=CompBMR, lambda='ML') #The branch length transformations can be optimised between bounds using maximum likelihood by setting the value for a transformation to 'ML'
summary(m_extra)