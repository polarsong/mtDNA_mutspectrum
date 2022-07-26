rm(list = ls(all=TRUE))
library(ggplot2)
df_mtdna = read.csv('../../Body/3Results/Birds_mtDNA_data.csv')
df_mtdna$ghahSkew = (df_mtdna$neutral_c - df_mtdna$neutral_T)/(df_mtdna$neutral_c + df_mtdna$neutral_T)
cox1 = subset(df_mtdna, df_mtdna$gene_name == 'COX1',)
cox2 = subset(df_mtdna, df_mtdna$gene_name == 'COX2',)
cox3 = subset(df_mtdna, df_mtdna$gene_name == 'COX3',)
atp8 = subset(df_mtdna, df_mtdna$gene_name == 'ATP8',)
atp6 = subset(df_mtdna, df_mtdna$gene_name == 'ATP6',)
nd2 = subset(df_mtdna, df_mtdna$gene_name == 'ND2',)
nd3 = subset(df_mtdna, df_mtdna$gene_name == 'ND3',)
nd4l = subset(df_mtdna, df_mtdna$gene_name == 'ND4L',)
nd4 = subset(df_mtdna, df_mtdna$gene_name == 'ND4',)
nd5 = subset(df_mtdna, df_mtdna$gene_name == 'ND5',)
cytb = subset(df_mtdna, df_mtdna$gene_name == 'CYTB',)
nd6 = subset(df_mtdna, df_mtdna$gene_name == 'ND6',)
nd1 = subset(df_mtdna, df_mtdna$gene_name == 'ND1',)

skew_all = ggplot(data = df_mtdna, aes(x = gene_name, y = ghahSkew))+
  geom_boxplot()
skew_all = skew_all + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
skew_all

skew_eco = ggplot(data = df_mtdna, aes(x = realm, y = ghahSkew))+
  geom_boxplot()
skew_eco = skew_eco + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
skew_eco
