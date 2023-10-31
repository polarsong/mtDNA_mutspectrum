# devtools::install_github('r-lib/systemfonts')
# BiocManager::install("svglite", dependencies = TRUE)


rm(list=ls(all=T))

library(ggtree)
library(ggplot2)
library(treeio)
library(dplyr)
library(ape)
library(geiger)
peng_data_tree = read.tree('peng_and_fly_all_data.tre')
birds_tree = read.tree('anc_kg.treefile')
df_mtdna = read.csv('Birds_dataset_paper.csv')

df_need = data.frame()
for (i in unique(df_mtdna$species_name))
{
  a = df_mtdna[df_mtdna$species_name == i,]
  b = sum(a$chthSkew)/12
  ab = c(i,b)
  df_need = rbind(df_need, ab)
}

names(df_need) = c('species_name', 'ChThSkew')
df_need$species_name = gsub(' ', '_', df_need$species_name)

df_flight = read.csv('flight_and_gene.csv')
df_flight = df_flight[,c(2,3,4,5,6)]
df_flight$flight_num = 0
df_flight[df_flight$ability_to_fly != 'Flying',]$flight_num = 1
rownames(df_flight) = df_flight$species_name
flying_tree = read.tree('flying_birds_tree.tre')
name.check(flying_tree, df_flight)
df_flight$ability_to_fly = as.character(df_flight$ability_to_fly)
df_flight$GhAhSkew = as.numeric(as.character(df_flight$GhAhSkew))
df_flight$fAn = as.numeric(as.character(df_flight$fAn))
df_flight$fGn = as.numeric(as.character(df_flight$fGn))
df_flight$flight_num = as.numeric(as.character(df_flight$flight_num))

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
df_flight1 = merge(df_flight, df_ff1)
df_flight1 = merge(df_flight1, df_syn1)

need_species = setdiff(df_flight$species_name, df_flight1$species_name)
correct_need_species = setdiff(df_flight$species_name, need_species)
birds_ms_tree<-keep.tip(flying_tree,correct_need_species)

df_flight1$GhAhSkew = as.numeric(df_flight1$GhAhSkew)
df_flight1$Mutation_AG = as.numeric(df_flight1$Mutation_AG)
df_flight1$Mutation_AG_syn = as.numeric(df_flight1$Mutation_AG_syn)
df_flight1$Mutation_CT = as.numeric(df_flight1$Mutation_CT)
df_flight1$Mutation_CT_syn = as.numeric(df_flight1$Mutation_CT_syn)
rownames(df_flight1)= df_flight1$species_name
df_flight1$double_fly = 0
df_flight1[df_flight1$ability_to_fly == 'Flying',]$double_fly = 1
df_temp_fly = read.csv('birds_metrics.csv')
df_temp_fly = df_temp_fly[,c(2:17)]
df_flight2 = merge(df_flight1, df_temp_fly)
need_species = setdiff(df_flight1$species_name, df_flight2$species_name)
correct_need_species = setdiff(df_flight1$species_name, need_species)
birds_ms_and_temp_tree<-keep.tip(birds_ms_tree,correct_need_species)
row.names(df_flight2) = df_flight2$species_name
df_flight2$Latitude = gsub(',', '.', df_flight2$Latitude)
df_flight2$TempRange = gsub(',', '.', df_flight2$TempRange)
df_flight2$AnnualPrecip = gsub(',', '.', df_flight2$AnnualPrecip)
df_flight2$PrecipRange= gsub(',', '.', df_flight2$PrecipRange)
df_flight2$AnnualTemp= gsub(',', '.', df_flight2$AnnualTemp)

df_flight2$AnnualTemp = as.numeric(as.character(df_flight2$AnnualTemp))
df_flight2$Beak_length_Culmen = as.numeric(as.character(df_flight2$Beak_length_Culmen))
df_flight2$Beak_length_Nares = as.numeric(as.character(df_flight2$Beak_length_Nares))
df_flight2$Beak_width = as.numeric(as.character(df_flight2$Beak_width))
df_flight2$Beak_depth = as.numeric(as.character(df_flight2$Beak_depth))
df_flight2$Tarsus_length = as.numeric(as.character(df_flight2$Tarsus_length))
df_flight2$Wing_length = as.numeric(as.character(df_flight2$Wing_length))
df_flight2$Kipps_distance = as.numeric(as.character(df_flight2$Kipps_distance))
df_flight2$Hand_wing_index = as.numeric(as.character(df_flight2$Hand_wing_index))
df_flight2$Tail_length = as.numeric(as.character(df_flight2$Tail_length))
df_flight2$Mass = as.numeric(as.character(df_flight2$Mass))
df_flight2$Latitude = as.numeric(as.character(df_flight2$Latitude))
df_flight2$AnnualTemp = as.numeric(as.character(df_flight2$AnnualTemp))
df_flight2$TempRange = as.numeric(as.character(df_flight2$TempRange))
df_flight2$AnnualPrecip = as.numeric(as.character(df_flight2$AnnualPrecip))
df_flight2$PrecipRange = as.numeric(as.character(df_flight2$PrecipRange))

df_flight2$Mass = log10(df_flight2$Mass)
df_flight2$Beak_length_Culmen = log10(df_flight2$Beak_length_Culmen)
df_flight2$Beak_length_Nares = log10(df_flight2$Beak_length_Nares)
df_flight2$Beak_width = log10(df_flight2$Beak_width)
df_flight2$Beak_depth = log10(df_flight2$Beak_depth)
df_flight2$Tarsus_length = log10(df_flight2$Tarsus_length)
df_flight2$Wing_length = log10(df_flight2$Wing_length)
df_flight2$Kipps_distance = log10(df_flight2$Kipps_distance)
df_flight2$Hand_wing_index = log10(df_flight2$Hand_wing_index)
df_flight2$Tail_length = log10(df_flight2$Tail_length)

rm(df_for_color1)
df_for_color = df_mutspec[df_mutspec$Label == 'syn',]
df_for_color = df_for_color[,c(1,5,7)]
df_for_color1 = reshape(data = df_for_color, idvar = 'AltNode',
                  timevar = 'Mut',
                  direction = 'wide')
df_for_color1$label = 'A>G'
df_for_color1 = df_for_color1[,c(1,14,12)]
names(df_for_color1) = c('species_name', 'label','Mutation_AG_syn')
df_for_color1$Mutation_AG_syn = as.numeric(df_for_color1$Mutation_AG_syn)
typeof(df_for_color1$Mutation_AG_syn)
df_for_color1$Mutation_AG_syn = df_for_color1$Mutation_AG_syn*10000000
df_for_color1$Mutation_AG_syn = round(df_for_color1$Mutation_AG_syn)
# function_name <- function(arg_1, arg_2, ...) {
#   Function body 
# }
typeof(df_for_color1$Mutation_AG_syn)
#my_tree = read.tree("./data/interim/iqtree_runs/drun1/anc.treefile")

#plot_sbs_on_tree <- function(label, indir='./data/processed/sbs_on_tree/', outdir='./figures/') {
#  mutspec_file <- paste(indir, label, '_ff.tsv', sep="")
#  print(mutspec_file)
#  grad_val = read.table(mutspec_file, sep="\t", header = T)
#  names(grad_val) = c('RefNode', 'label', 'MutSpec')  # this label is not the label in input args
row.names(df_for_color1) = df_for_color1$species_name
test = full_join(as_tibble(birds_tree), tibble(df_for_color1[,c(2,3)]), by = 'label')
der = as.treedata(test)
?tibble
#der <- groupOTU(der, c("Halicephalobus_mephisto", "Caenorhabditis_elegans"))
ggtree(der, aes(color = der@data$species_name)) + geom_tiplab(size = 3, colour = "darkgray") +
  scale_color_continuous(low="green", high="red") +
  theme(legend.position="bottom") +
  labs(col=label) +
  geom_tippoint(aes(alpha = group), col = "red") +
  scale_color_manual(values = c(0,1), aesthetics = "alpha")
  

ggsave(out_image_file,  width = 2000, height = 2000, units = 'px', scale=2.2)
#}

label <- "T>G"
plot_sbs_on_tree(label)


for (n1 in c('A', 'C', 'G', 'T')) {
  for (n2 in c('A', 'C', 'G', 'T')) {
    if (n1 == n2) {
      next
    }
    label <- paste(n1, n2, sep = ">")
    print(label)
    plot_sbs_on_tree(label)
  }
}

# A>T / T>A
label = 'A>T_T>A_ratio'
plot_sbs_on_tree(label)




# T on tree
label <- 'T share'
mutspec_file <- './data/processed/percent_of_T.tsv'
# out_image_file <- './figures/T_ff_tree.svg'
out_image_file <- './figures/T_tree.png'


grad_val = read.table('percent_of_T.tsv', sep="\t", header = T)
names(grad_val) = c('RefNode', 'label', 'MutSpec')

test = full_join(as_tibble(my_tree), tibble(grad_val[,c(2,3)]), by = 'label')
der = as.treedata(test)

# svg(out_image_file, width = 20, height = 20)
png(out_image_file, width = 2000, height = 2000)
der <- groupOTU(der, c("Halicephalobus_mephisto", "Caenorhabditis_elegans"))
ggtree(der, aes(color = der@data$MutSpec)) + geom_tiplab(size = 4, colour = "darkgray") +
  scale_color_continuous(low="green", high="red") +
  theme(legend.position="bottom")+
  labs(col=label)+
  geom_tippoint(aes(alpha = group), col = "purple")+
  scale_color_manual(values = c(0,1), aesthetics = "alpha")

dev.off()