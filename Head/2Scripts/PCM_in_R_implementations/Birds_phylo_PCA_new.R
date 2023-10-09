rm(list = ls(all=TRUE))
library(ape)
library(phangorn)
library(phytools)
library(geiger)
#birds_tree<-read.tree(file="../../2Scripts/PCM_in_R_implementations/anc_kg.treefile")
df_mtdna = read.csv('Birds_dataset_paper.csv')
df_realms = unique(df_mtdna[,c('species_name', 'realm')])
df_realms$species_name = gsub(' ', '_', df_realms$species_name)
#df_temp = read.csv('Birds_temperature_table.csv')
#df_char1= unique(df_mtdna[,c('species_name',"Beak_length_Culmen", "Beak_length_Nares", "Beak_width", "Beak_depth", "Tarsus_length", "Wing_length", "Kipps_distance", "Hand_wing_index", "Tail_length", "Mass" )])
#df_char2 = df_temp[,c('Species.name', "Latitude", "AnnualTemp", "TempRange", "AnnualPrecip", "PrecipRange")]
#names(df_char2)= c('species_name', "Latitude", "AnnualTemp", "TempRange", "AnnualPrecip", "PrecipRange")
#df_char1 = df_char1[df_char1$species_name != "Agapornis pullarius",]
#df_char1 = df_char1[df_char1$species_name != "Mergus squamatus",]
#df_char3 = merge(df_char1, df_char2)
#df_char1$species_name = gsub(' ', '_', df_char1$species_name)
#df_char3$species_name = gsub(' ', '_', df_char3$species_name)
#df_char3 = na.omit(df_char3)
#need_species = setdiff(df_char1$species_name, df_char3$species_name)
#correct_need_species = setdiff(df_char1$species_name, need_species)
#temp_tree<-keep.tip(birds_tree,correct_need_species)
#temp_tree
#write.tree(temp_tree ,file="temperature_birds_tree.tre")
#write.csv(df_char3, 'birds_metrics.csv')
df_temp_fly = read.csv('birds_metrics.csv')
rownames(df_temp_fly) = df_temp_fly$species_name
df_flight_names = read.csv('flight_and_gene.csv')
df_flight_names = df_flight_names[,c(2:3)]
df_flight_names = merge(df_flight_names, df_realms)
rownames(df_flight_names) = df_flight_names$species_name
df_temp_fly = df_temp_fly[,c(2:17)]
df_temp_fly = merge(df_temp_fly, df_flight_names)
rownames(df_temp_fly) = df_temp_fly$species_name
#df_temp_fly = df_temp_fly[,c(2:17)]


temp_tree = read.tree('flight_and_temp.tre')
df_temp_fly$Latitude = gsub(',', '.', df_temp_fly$Latitude)
df_temp_fly$AnnualTemp = gsub(',', '.', df_temp_fly$AnnualTemp)
df_temp_fly$TempRange = gsub(',', '.', df_temp_fly$TempRange)
df_temp_fly$AnnualPrecip = gsub(',', '.', df_temp_fly$AnnualPrecip)
df_temp_fly$PrecipRange= gsub(',', '.', df_temp_fly$PrecipRange)

df_temp_fly$Beak_length_Culmen = as.numeric(as.character(df_temp_fly$Beak_length_Culmen))
df_temp_fly$Beak_length_Nares = as.numeric(as.character(df_temp_fly$Beak_length_Nares))
df_temp_fly$Beak_width = as.numeric(as.character(df_temp_fly$Beak_width))
df_temp_fly$Beak_depth = as.numeric(as.character(df_temp_fly$Beak_depth))
df_temp_fly$Tarsus_length = as.numeric(as.character(df_temp_fly$Tarsus_length))
df_temp_fly$Wing_length = as.numeric(as.character(df_temp_fly$Wing_length))
df_temp_fly$Kipps_distance = as.numeric(as.character(df_temp_fly$Kipps_distance))
df_temp_fly$Hand_wing_index = as.numeric(as.character(df_temp_fly$Hand_wing_index))
df_temp_fly$Tail_length = as.numeric(as.character(df_temp_fly$Tail_length))
df_temp_fly$Mass = as.numeric(as.character(df_temp_fly$Mass))
df_temp_fly$Latitude = as.numeric(as.character(df_temp_fly$Latitude))
df_temp_fly$AnnualTemp = as.numeric(as.character(df_temp_fly$AnnualTemp))
df_temp_fly$TempRange = as.numeric(as.character(df_temp_fly$TempRange))
df_temp_fly$AnnualPrecip = as.numeric(as.character(df_temp_fly$AnnualPrecip))
df_temp_fly$PrecipRange = as.numeric(as.character(df_temp_fly$PrecipRange))

df_temp_fly$Beak_length_Culmen = log10(df_temp_fly$Beak_length_Culmen)
df_temp_fly$Beak_length_Nares = log10(df_temp_fly$Beak_length_Nares)
df_temp_fly$Beak_width = log10(df_temp_fly$Beak_width)
df_temp_fly$Beak_depth = log10(df_temp_fly$Beak_depth)
df_temp_fly$Tarsus_length = log10(df_temp_fly$Tarsus_length)
df_temp_fly$Wing_length = log10(df_temp_fly$Wing_length)
df_temp_fly$Kipps_distance = log10(df_temp_fly$Kipps_distance)
df_temp_fly$Hand_wing_index = log10(df_temp_fly$Hand_wing_index)
df_temp_fly$Tail_length = log10(df_temp_fly$Tail_length)
df_temp_fly$Mass = log10(df_temp_fly$Mass)


df_temp_fly$TempRange = log10(df_temp_fly$TempRange)
df_temp_fly$AnnualPrecip = log10(df_temp_fly$AnnualPrecip)
df_temp_fly$PrecipRange = log10(df_temp_fly$PrecipRange)

#adding everything to PCA
df_flight = read.csv('flight_and_gene.csv')
df_flight = df_flight[,c(2,4,5,6)]
df_temp_fly = merge(df_temp_fly, df_flight)

df_mutspec = read.table('C:/Users/User/Desktop/mutspec12.tsv', header = TRUE, fill = TRUE)
df_ff = df_mutspec[df_mutspec$Label == 'ff',]
df_syn = df_mutspec[df_mutspec$Label == 'syn',]
df_ff = df_ff[!grepl('Node', df_ff$AltNode),]
df_ff = df_ff[,c(1,5,7)]
df_ff1 = reshape(data = df_ff, idvar = 'AltNode',
                 timevar = 'Mut',
                 direction = 'wide')
names(df_ff1) = c('species_name', 'Mutation_TG', 'Mutation_TC', 'Mutation_TA', 'Mutation_GT', 'Mutation_GC', 
                  'Mutation_GA', 'Mutation_CT', 'Mutation_CG', 'Mutation_CA', 'Mutation_AT', 'Mutation_AG', 'Mutation_AC')
df_syn = df_syn[!grepl('Node', df_syn$AltNode),]
df_syn = df_syn[,c(1,5,7)]
df_syn1 = reshape(data = df_syn, idvar = 'AltNode',
                  timevar = 'Mut',
                  direction = 'wide')
names(df_syn1) = c('species_name', 'Mutation_TG_syn', 'Mutation_TC_syn', 'Mutation_TA_syn', 'Mutation_GT_syn', 'Mutation_GC_syn', 
                   'Mutation_GA_syn', 'Mutation_CT_syn', 'Mutation_CG_syn', 'Mutation_CA_syn', 'Mutation_AT_syn', 'Mutation_AG_syn', 'Mutation_AC_syn')

df_temp_fly1 = merge(df_temp_fly, df_ff1)
df_temp_fly1 = merge(df_temp_fly1, df_syn1)

need_species = setdiff(df_temp_fly$species_name, df_temp_fly1$species_name)
correct_need_species = setdiff(df_temp_fly$species_name, need_species)
big_data_tree1<-keep.tip(temp_tree,correct_need_species)

df_temp_fly1$GhAhSkew = as.numeric(as.character(df_temp_fly1$GhAhSkew))
df_temp_fly1$fAn = as.numeric(as.character(df_temp_fly1$fAn))
df_temp_fly1$fGn = as.numeric(as.character(df_temp_fly1$fGn))
df_temp_fly1$Mutation_TG = as.numeric(as.character(df_temp_fly1$Mutation_TG))
df_temp_fly1$Mutation_TC = as.numeric(as.character(df_temp_fly1$Mutation_TC))
df_temp_fly1$Mutation_TA = as.numeric(as.character(df_temp_fly1$Mutation_TA))
df_temp_fly1$Mutation_GT = as.numeric(as.character(df_temp_fly1$Mutation_GT))
df_temp_fly1$Mutation_GA = as.numeric(as.character(df_temp_fly1$Mutation_GA))
df_temp_fly1$Mutation_GC = as.numeric(as.character(df_temp_fly1$Mutation_GC))
df_temp_fly1$Mutation_CT = as.numeric(as.character(df_temp_fly1$Mutation_CT))
df_temp_fly1$Mutation_CG = as.numeric(as.character(df_temp_fly1$Mutation_CG))
df_temp_fly1$Mutation_CA = as.numeric(as.character(df_temp_fly1$Mutation_CA))
df_temp_fly1$Mutation_AT = as.numeric(as.character(df_temp_fly1$Mutation_AT))
df_temp_fly1$Mutation_AG = as.numeric(as.character(df_temp_fly1$Mutation_AG))
df_temp_fly1$Mutation_AC = as.numeric(as.character(df_temp_fly1$Mutation_AC))

df_temp_fly1$Mutation_TG_syn = as.numeric(as.character(df_temp_fly1$Mutation_TG_syn))
df_temp_fly1$Mutation_TC_syn = as.numeric(as.character(df_temp_fly1$Mutation_TC_syn))
df_temp_fly1$Mutation_TA_syn = as.numeric(as.character(df_temp_fly1$Mutation_TA_syn))
df_temp_fly1$Mutation_GT_syn = as.numeric(as.character(df_temp_fly1$Mutation_GT_syn))
df_temp_fly1$Mutation_GA_syn = as.numeric(as.character(df_temp_fly1$Mutation_GA_syn))
df_temp_fly1$Mutation_GC_syn = as.numeric(as.character(df_temp_fly1$Mutation_GC_syn))
df_temp_fly1$Mutation_CT_syn = as.numeric(as.character(df_temp_fly1$Mutation_CT_syn))
df_temp_fly1$Mutation_CG_syn = as.numeric(as.character(df_temp_fly1$Mutation_CG_syn))
df_temp_fly1$Mutation_CA_syn = as.numeric(as.character(df_temp_fly1$Mutation_CA_syn))
df_temp_fly1$Mutation_AT_syn = as.numeric(as.character(df_temp_fly1$Mutation_AT_syn))
df_temp_fly1$Mutation_AG_syn = as.numeric(as.character(df_temp_fly1$Mutation_AG_syn))
df_temp_fly1$Mutation_AC_syn = as.numeric(as.character(df_temp_fly1$Mutation_AC_syn))
row.names(df_temp_fly1) = df_temp_fly1$species_name
#SET ALMOST EVERYTHING TO LOG10
temp_birds_pca<-phyl.pca(big_data_tree1,df_temp_fly1[,c(2:11, 13:16, 19:45)])
par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(temp_birds_pca,main="")
#flying_birds_pca$Evec[,1]<--flying_birds_pca$Evec[,1]
#flying_birds_pca$L[,1]<--flying_birds_pca$L[,1]
#flying_birds_pca$S<-scores(flying_birds_pca,
#newdata=df_flight_pca)

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))

phylomorphospace(big_data_tree1,
                 scores(temp_birds_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1",
                 ylab="PC2")
#color legend!!!! do later
df_temp_fly1$ability_to_fly = as.factor(df_temp_fly1$ability_to_fly)
df_temp_fly1$species_name = as.factor(df_temp_fly1$species_name)
df_temp_fly1$realm = as.factor(df_temp_fly1$realm)
eco<-setNames(df_temp_fly1[,17],rownames(df_temp_fly1))

ECO<-to.matrix(eco,levels(eco))
tiplabels(pie=ECO[big_data_tree1$tip.label,],cex=0.3)
legend(x="topleft",legend=levels(eco),cex=0.5,pch=21,
       pt.bg=rainbow(n=length(levels(eco))),pt.cex=1.5)

a = as.data.frame(temp_birds_pca$S)
a$PC1 = a$PC1 - min(a$PC1) + 1
a$PC2 = a$PC2 - min(a$PC2) + 1
a$PC3 = a$PC3 - min(a$PC3) + 1
a$PC4 = a$PC4 - min(a$PC4) + 1
a$PC5 = a$PC5 - min(a$PC5) + 1
a$PC6 = a$PC6 - min(a$PC6) + 1
a$PC7 = a$PC7 - min(a$PC7) + 1
a$PC8 = a$PC8 - min(a$PC8) + 1
a$PC9 = a$PC9 - min(a$PC9) + 1
a$PC10 = a$PC10 - min(a$PC10) + 1

a$PC1 = log10(a$PC1)
a$PC2 = log10(a$PC2)
a$PC3 = log10(a$PC3)
a$PC4 = log10(a$PC4)
a$PC5 = log10(a$PC5)
a$PC6 = log10(a$PC6)
a$PC7 = log10(a$PC7)
a$PC8 = log10(a$PC8)
a$PC9 = log10(a$PC9)
a$PC10 = log10(a$PC10)

phylomorphospace(temp_tree,
                 a[,9:10],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC9",
                 ylab="PC10")
#color legend!!!! do later
df_temp_fly$ability_to_fly = as.factor(df_temp_fly$ability_to_fly)
df_temp_fly$species_name = as.factor(df_temp_fly$species_name)
eco<-setNames(df_temp_fly[,17],rownames(df_temp_fly))

ECO<-to.matrix(eco,levels(eco))
tiplabels(pie=ECO[temp_tree$tip.label,],cex=0.3)
legend(x="topleft",legend=levels(eco),cex=0.6,pch=21,
       pt.bg=rainbow(n=length(levels(eco))),pt.cex=1.5)

#only non-flying
df_nf = df_temp_fly[df_temp_fly$ability_to_fly != 'Flying',]
#not_fly_species = setdiff(df_temp_fly$species_name, df_nf$species_name)
#correct_not_fly_species = setdiff(df_temp_fly$species_name, not_fly_species)
#not_fly_tree<-keep.tip(temp_tree,correct_not_fly_species)
#not_fly_tree
write.tree(not_fly_tree ,file="non_flying_birds_tree.tre")

nf_tree = read.tree('non_flying_birds_tree.tre')
nf_birds_pca<-phyl.pca(nf_tree,df_nf[,c(2:16)])
par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(nf_birds_pca,main="")
#flying_birds_pca$Evec[,1]<--flying_birds_pca$Evec[,1]
#flying_birds_pca$L[,1]<--flying_birds_pca$L[,1]
#flying_birds_pca$S<-scores(flying_birds_pca,
#newdata=df_flight_pca)

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))

phylomorphospace(nf_tree,
                 scores(nf_birds_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1",
                 ylab="PC2")
#color legend!!!! do later
df_nf$ability_to_fly = as.factor(df_nf$ability_to_fly)
df_nf$species_name = as.factor(df_nf$species_name)
df_nf$realm = as.factor(df_nf$realm)
eco<-setNames(df_nf[,17],rownames(df_nf))

ECO<-to.matrix(eco,levels(eco))
tiplabels(pie=ECO[nf_tree$tip.label,],cex=0.3)
legend(x="bottomright",legend=levels(eco),cex=0.6,pch=21,
       pt.bg=rainbow(n=length(levels(eco))),pt.cex=1.5)
