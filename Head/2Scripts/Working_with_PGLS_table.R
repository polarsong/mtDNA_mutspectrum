rm(list = ls(all=TRUE))
library(ggplot2)
library(heatmaply)


#reading table
df = read.csv('../../Head/2Scripts/Birds_all_PGLS.csv')

#adjusting p-values
df$p_adj_bon = p.adjust(df$p_value, method = 'bonferroni')
df$p_adj_holm = p.adjust(df$p_value, method = 'holm')
df$p_adj_hochberg = p.adjust(df$p_value, method = 'hochberg')
df$p_adj_hommel = p.adjust(df$p_value, method = 'hommel')
df$p_adj_bh = p.adjust(df$p_value, method = 'BH')
df$p_adj_fdr = p.adjust(df$p_value, method = 'fdr')
df$p_adj_by = p.adjust(df$p_value, method = 'BY')
df$p_adj_none = p.adjust(df$p_value, method = 'none')

#discarding p-values
df_disc = df[df$p_value <= 0.05,]
df_disc = df_disc[is.na(df_disc$p_value) == FALSE,]
df_disc = df_disc[df_disc$p_value != 0,]

#genetics df
df1 = df_disc[df_disc$Ecology1 =='ghahSkew',]
df2 = df_disc[df_disc$Ecology1 =='Stg_Sac',]
df3 = df_disc[df_disc$Ecology1 =='med_c',]
df4 = df_disc[df_disc$Ecology1 =='med_a',]
df5 = df_disc[df_disc$Ecology1 =='fTn',]
df6 = df_disc[df_disc$Ecology1 =='fCn',]
df7 = df_disc[df_disc$Ecology1 =='fAn',]
df8 = df_disc[df_disc$Ecology1 =='fGn',]
df_genetic = rbind(df1, df2, df3, df4, df5, df6, df7, df8)

#birds measurements df
df1 = df_disc[df_disc$Ecology1 =='Beak_length_Culmen',]
df2 = df_disc[df_disc$Ecology1 =='Beak_length_Nares',]
df3 = df_disc[df_disc$Ecology1 =='Beak_width',]
df4 = df_disc[df_disc$Ecology1 =='Beak_depth',]
df5 = df_disc[df_disc$Ecology1 =='Tarsus_length',]
df6 = df_disc[df_disc$Ecology1 =='Wing_length',]
df7 = df_disc[df_disc$Ecology1 =='Kipps_distance',]
df8 = df_disc[df_disc$Ecology1 =='Hand_wing_index',]
df9 = df_disc[df_disc$Ecology1 =='Tail_length',]
df10 = df_disc[df_disc$Ecology1 =='Mass',]
df_measurements = rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10)

#birds_lifestyles df
df_lifestyles = df_disc[!df_disc$Ecology1 %in% df_genetic$Ecology1,]
df_lifestyles = df_lifestyles[!df_lifestyles$Ecology1 %in% df_measurements$Ecology1,]

#drawing graphs
#genetic
g1 = ggplot(df_genetic, aes(x = Ecology1, y = Ecology2, fill = p_value))+
  geom_tile()
g1 = g1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
g1
g2 = ggplot(df_genetic, aes(x = Ecology1, y = Ecology2, fill = mult_r_squared))+
  geom_tile()
g2 = g2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
g2
g3 = ggplot(df_genetic, aes(x = Ecology1, y = Ecology2, fill = effect_size))+
  geom_tile()
g3 = g3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
g3

#measurements
m1 = ggplot(df_measurements, aes(x = Ecology1, y = Ecology2, fill = p_value))+
  geom_tile()
m1 = m1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
m1
m2 = ggplot(df_measurements, aes(x = Ecology1, y = Ecology2, fill = mult_r_squared))+
  geom_tile()
m2 = m2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
m2
m3 = ggplot(df_measurements, aes(x = Ecology1, y = Ecology2, fill = effect_size))+
  geom_tile()
m3 = m3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
m3

#lifestyles
l1 = ggplot(df_lifestyles, aes(x = Ecology1, y = Ecology2, fill = p_value))+
  geom_tile()
l1 = l1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
l1
l2 = ggplot(df_lifestyles, aes(x = Ecology1, y = Ecology2, fill = mult_r_squared))+
  geom_tile()
l2 = l2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
l2
l3 = ggplot(df_lifestyles, aes(x = Ecology1, y = Ecology2, fill = effect_size))+
  geom_tile()
l3 = l3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
l3







