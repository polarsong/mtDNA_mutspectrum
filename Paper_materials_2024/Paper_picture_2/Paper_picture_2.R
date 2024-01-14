rm(list = ls(all=TRUE))
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggbiplot)
library(phytools)
library(nlme)
library(geiger)
library(ggtree)
library(stringr)
library(dplyr)

df_mtdna = read.csv('../Paper_materials_2024/Birds_dataset_paper.csv')

df_nd6 = read.csv('../Paper_materials_2024/Birds_mtDNA_data.csv')
df_nd6$GhAhSkew = (df_nd6$neutral_g - df_nd6$neutral_A)/(df_nd6$neutral_g + df_nd6$neutral_A)
df_nd6$ThChSkew = (df_nd6$neutral_T - df_nd6$neutral_c)/(df_nd6$neutral_T + df_nd6$neutral_c)
df_nd6$fTn = df_nd6$neutral_T/df_nd6$neutral_amount
df_nd6$fAn = df_nd6$neutral_A/df_nd6$neutral_amount
df_nd6$fCn = df_nd6$neutral_c/df_nd6$neutral_amount
df_nd6$fGn = df_nd6$neutral_g/df_nd6$neutral_amount

df_nd6$GhAhSkew = (df_nd6$neutral_c- df_nd6$neutral_T)/(df_nd6$neutral_c + df_nd6$neutral_T)
df_nd6$ThChSkew = (df_nd6$neutral_A - df_nd6$neutral_g)/(df_nd6$neutral_A + df_nd6$neutral_g)
df_nd6$fTn = df_nd6$neutral_A/df_nd6$neutral_amount
df_nd6$fAn = df_nd6$neutral_T/df_nd6$neutral_amount
df_nd6$fCn = df_nd6$neutral_g/df_nd6$neutral_amount
df_nd6$fGn = df_nd6$neutral_c/df_nd6$neutral_amount

names_v = unique(df_mtdna$species_name)
df_short = data.frame()
for (i in names_v)
{
  df1 = df_mtdna[df_mtdna$species_name == i,]
  a = sum(df1$ghahSkew)/12
  b = sum(df1$chthSkew)/12
  v = sum(df1$Mass)/12
  ab = c(i, a, b, v)
  df_short = rbind(df_short, ab)
}
names(df_short) = c('species_name', 'GhAhSkew', 'ThChSkew', 'Mass')
df_short$Mass = as.numeric(df_short$Mass)
df_short$GhAhSkew = as.numeric(df_short$GhAhSkew)
df_short$ThChSkew = as.numeric(df_short$ThChSkew)
df_short$log_mass = log10(df_short$Mass)
mass_skew = ggplot(df_short, aes(x = log_mass, y = GhAhSkew))+
  geom_point()+
  geom_smooth(method=lm)+
  xlab('Decimal logarithm of mass')

skew_eco = ggplot(data = df_mtdna, aes(x = realm, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Birds realms')+
  ylab('GhAhSkew')+
  xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

skew_niche = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic niche')+
  ylab('GhAhSkew')+
  xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

df_pca = df_mtdna[c('species_name','gene_name','fAn', 'fGn', 'fCn', 'fTn', 'Stg_Sac','ghahSkew', 'chthSkew')]
gene_vector = c('fAn', 'fGn', 'fCn', 'fTn', 'Stg_Sac','ghahSkew', 'chthSkew')
gene_stats = data.frame(unique(df_pca$species_name))
for ( char in gene_vector){
  
  stats1 = aggregate(df_pca[,char], by = list(df_pca$species_name), FUN = 'sum')[2]
  stats1 = stats1/12
  gene_stats = cbind(gene_stats, stats1)
  
}
names(gene_stats) = c('species_name', gene_vector)
df_realm = df_mtdna[c('species_name', 'realm', 'Trophic_niche')]
gene_stats = merge(gene_stats, df_realm, by = 'species_name')
gene_stats = unique(gene_stats)
row.names(gene_stats) = gene_stats$species_name
gene_stats = gene_stats[, colSums(is.na(gene_stats)) < nrow(gene_stats)]
stats_pca = prcomp(gene_stats[c(2,3,4,5,6,7,8)], center = TRUE, scale. = TRUE)
summary(stats_pca)
PCA_mtDNA_metrics = autoplot(stats_pca,
                             loadings = TRUE,
                             loadings.label = TRUE)

birds_pca = data.frame(stats_pca$x)
birds_pca = birds_pca[,c(1,2)]
birds_pca$species_name = row.names(birds_pca)
gene_stats$species_name = row.names(gene_stats)
gene_stats = merge(gene_stats, birds_pca, by = 'species_name')
row.names(gene_stats) = gene_stats$species_name
gene_stats = gene_stats[,c(2:12)]
Density_PC1_eco = ggplot(gene_stats, aes(x=PC1, color=realm)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')+
  ylab('Density')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))
Density_PC1_niche = ggplot(gene_stats, aes(x=PC1, color=Trophic_niche)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')+
  ylab('Density')+
  scale_colour_manual(name="Origin", values= c("black", "black", "black", "red", "black", "black", "black", "black", "black", 'black'))

mass_skew
skew_eco
skew_niche
PCA_mtDNA_metrics
Density_PC1_eco
Density_PC1_niche
