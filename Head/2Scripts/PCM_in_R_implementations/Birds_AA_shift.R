rm(list = ls(all=TRUE))
library(ggplot2)
df_mtdna = read.csv('Birds_dataset_paper.csv')
codon_names = c('ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC', 'GTA', 'GTG', 'ACT', 'ACC', 'ACA', 'ACG',
                'GCT', 'GCC', 'GCA', 'GCG', 'AAT', 'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG',
                'AGT', 'AGC', 'AGA', 'AGG', 'GGT', 'GGC', 'GGA', 'GGG')
#first position
df_mtdna$Ile = df_mtdna$TAA + df_mtdna$TAG
df_mtdna$Met = df_mtdna$TAT + df_mtdna$TAC
df_mtdna$Val_Ile = df_mtdna$CAA + df_mtdna$CAG
df_mtdna$Val_Met = df_mtdna$CAT + df_mtdna$CAC
df_mtdna$Thr = df_mtdna$TGA + df_mtdna$TGG + df_mtdna$TGT + df_mtdna$TGC
df_mtdna$Ala = df_mtdna$CGA + df_mtdna$CGG + df_mtdna$CGT + df_mtdna$CGC
df_mtdna$Asn = df_mtdna$TTA + df_mtdna$TTG
df_mtdna$Lys = df_mtdna$TTT + df_mtdna$TTC
df_mtdna$Asp = df_mtdna$CTA + df_mtdna$CTG
df_mtdna$Glu = df_mtdna$CTT + df_mtdna$CTC
df_mtdna$Ser = df_mtdna$TCA + df_mtdna$TCG
df_mtdna$stop = df_mtdna$TCT + df_mtdna$TCC
df_mtdna$Gly_Ser = df_mtdna$CCA + df_mtdna$CCG
df_mtdna$Gly_stop = df_mtdna$CCT + df_mtdna$CCC

df_aa = df_mtdna[,c('species_name','Ile','Met','Val_Ile','Val_Met','Thr','Ala','Asn','Lys','Asp',
                    'Glu','Ser','stop','Gly_Ser', 'Gly_stop')]

df_need = data.frame()
for (i in unique(df_aa$species_name))
{
  a = df_aa[df_aa$species_name == i,]
  b = a[,c(2:15)]
  c = colSums(b)
  c$species_name = i
  df_need = rbind(df_need, c)
}
df_need1 = df_need[,c(1:14)]
df_need1 = as.data.frame(colSums(df_need1))
df_need1$AA = row.names(df_need1)
names(df_need1) = c('number', 'AA')
ggplot(df_need1, aes(x = AA, y=number))+
  geom_point() +
  xlim('Ile','Val_Ile','Met','Val_Met','Thr','Ala','Asn','Asp','Lys',
           'Glu','Ser','Gly_Ser','stop','Gly_stop')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
