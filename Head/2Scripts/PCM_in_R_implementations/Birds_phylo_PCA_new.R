rm(list = ls(all=TRUE))
library(ape)
library(phangorn)
library(phytools)
library(geiger)
#birds_tree<-read.tree(file="../../2Scripts/PCM_in_R_implementations/anc_kg.treefile")
#df_mtdna = read.csv('Birds_dataset_paper.csv')
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
df_temp_fly = df_temp_fly[,c(3:17)]
temp_tree = read.tree('temperature_birds_tree.tre')
df_temp_fly$Latitude = gsub(',', '.', df_temp_fly$Latitude)
df_temp_fly$AnnualTemp = gsub(',', '.', df_temp_fly$AnnualTemp)
df_temp_fly$TempRange = gsub(',', '.', df_temp_fly$TempRange)
df_temp_fly$AnnualPrecip = gsub(',', '.', df_temp_fly$AnnualPrecip)
df_temp_fly$PrecipRange= gsub(',', '.', df_temp_fly$PrecipRange)

df_temp_fly$Beak_length_Culmen = as.numeric(as.character(df_temp_fly$Beak_length_Culmen))
df_temp_fly$Beak_length_Naresn = as.numeric(as.character(df_temp_fly$Beak_length_Nares))
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



#SET ALMOST EVERYTHING TO LOG10
temp_birds_pca<-phyl.pca(temp_tree,df_temp_fly)
temp_birds_pca
par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(flying_birds_pca,main="")
#flying_birds_pca$Evec[,1]<--flying_birds_pca$Evec[,1]
#flying_birds_pca$L[,1]<--flying_birds_pca$L[,1]
#flying_birds_pca$S<-scores(flying_birds_pca,
#newdata=df_flight_pca)

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))

phylomorphospace(temp_tree,
                 scores(temp_birds_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1",
                 ylab="PC2")
#color legend!!!! do later

df_flight_names[df_flight_names$Flight == '0',]$Flight = 'Flying'
df_flight_names$Flight = as.factor(df_flight_names$Flight)
eco<-setNames(df_flight_names[,3],rownames(df_flight_names))

ECO<-to.matrix(eco,levels(eco))
tiplabels(pie=ECO[flight_tree$tip.label,],cex=0.3)
legend(x="center",legend=levels(eco),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(levels(eco))),pt.cex=1.5)
