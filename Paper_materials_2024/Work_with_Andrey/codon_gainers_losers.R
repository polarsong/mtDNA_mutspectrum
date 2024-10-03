rm(list = ls(all=TRUE))
library(ggplot2)
library(ggpubr)
df_mtdna = read.csv('../../Paper_materials_2024/Birds_dataset_paper.csv')
df_codons = df_mtdna[,c(2,32:95)]
df_fly = read.csv('../../Paper_materials_2024/flying_birds.csv')
df_codons1 = aggregate(. ~ species_name, FUN = mean, data = df_codons)

df_fly = df_fly[,c(2,3,4)]
names(df_fly) = c('species_name', 'flightless', 'diving')
df_fly_clean1 = df_fly[df_fly$flightless =='Flightless',]
df_fly_clean= df_fly[df_fly$flightless == 'Almost_flightless',]
df_fly_clean = na.omit(df_fly_clean)
df_fly_clean1 = na.omit(df_fly_clean1)
df_dive = df_fly
df_fly = df_fly[df_fly$flightless != 'Flightless',]
df_fly = df_fly[df_fly$flightless != 'Almost_flightless',]
df_fly_clean$flightless = 'Tinamiformes'
df_fly_clean1$flightless = 'Tinamiformes'
df_fly_big = rbind(df_fly, df_fly_clean, df_fly_clean1)


beasts = df_codons1$species_name
df_data = data.frame()
for (i in beasts)
{ 
  a = data.frame(df_codons1[df_codons1$species_name == i,])
  a$mega_gain_gly = sum(a$CCT + a$CCA + a$CCG + a$CCC) #mega_gain_gly
  a$gain_serg_arg = sum(a$TCT + a$TCA + a$TCG + a$TCC) #gain_serg_arg
  a$gain_asp_glu = sum(a$CTT + a$CTA + a$CTG + a$CTC) #gain_asp_glu
  a$lose_asn_lys = sum(a$TTT + a$TTA + a$TTG + a$TTC) #lose_asn_lys
  
  a$mega_gain_cys_stop_trp = sum(a$ACT + a$ACA + a$ACG + a$ACC) #mega_gain_cys_stop_trp
  a$gain_arg = sum(a$GCT + a$GCA + a$GCG + a$GCC) #gain_arg
  a$gain_tyr_stop = sum(a$ATT + a$ATA + a$ATG + a$ATC) #gain_tyr_stop
  a$lose_his_gln = sum(a$GTT + a$GTA + a$GTG + a$GTC) #lose_gis_gln
  
  a$mega_gain_val = sum(a$CAT + a$CAA + a$CAG + a$CAC) #mega_gain_val
  a$gain_ile_met = sum(a$TAT + a$TAA + a$TAG + a$TAC) #gain_ile_met
  a$gain_ala = sum(a$CGT + a$CGA + a$CGG + a$CGC) #gain ala
  a$lose_thr = sum(a$TGT + a$TGA + a$TGG + a$TGC) #lose_thr
  
  a$mega_gain_leu_phe = sum(a$AAT + a$AAA + a$AAG + a$AAC) #mega_gain_leu_phe
  a$gain_ser = sum(a$AGT + a$AGA + a$AGG + a$AGC) #gain_ser
  a$gain_leu = sum(a$GAT + a$GAA + a$GAG + a$GAC) #gain_leu
  a$lose_pro = sum(a$GGT + a$GGA + a$GGG + a$GGC)#lose pro
  df_data = rbind(df_data, a)
}

df_fly_final = merge(df_fly_big, df_data)
df_fly_final = df_fly_final[df_fly_final$flightless != 'Galliformes',]
df_fly_final[df_fly_final$flightless == '0',]$flightless = 'Flying birds'
a = ggplot(df_fly_final, aes(x = flightless, y = mega_gain_gly))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Flightless birds groups')+
  ylab('Gly')+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
b = ggplot(df_fly_final, aes(x = flightless, y = gain_serg_arg))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Flightless birds groups')+
  ylab('Ser and stop')+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
c = ggplot(df_fly_final, aes(x = flightless, y = gain_asp_glu))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Flightless birds groups')+
  ylab('Asp and glu')+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
d = ggplot(df_fly_final, aes(x = flightless, y = lose_asn_lys))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Flightless birds groups')+
  ylab('Asn and lys')+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
first_graph = ggarrange(d, b, c, a,
                        ncol = 2, nrow = 2)
first_graph

a1 =ggplot(df_fly_final, aes(x = flightless, y = mega_gain_cys_stop_trp))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Flightless birds groups')+
  ylab('Cys and Trp')+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
b1 = ggplot(df_fly_final, aes(x = flightless, y = gain_arg))+
  geom_boxplot()+
  xlab('Flightless birds groups')+
  ylab('Arg')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
c1 = ggplot(df_fly_final, aes(x = flightless, y = gain_tyr_stop))+
  geom_boxplot()+
  xlab('Flightless birds groups')+
  ylab('Tyr and stop')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
d1 = ggplot(df_fly_final, aes(x = flightless, y = lose_his_gln))+
  geom_boxplot()+
  xlab('Flightless birds groups')+
  ylab('His and gln')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
second_graph = ggarrange(c1,a1,d1,b1,
                         ncol = 2, nrow = 2)
second_graph

a2 = ggplot(df_fly_final, aes(x = flightless, y = mega_gain_val))+
  geom_boxplot()+
  xlab('Flightless birds groups')+
  ylab('Val')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
b2 = ggplot(df_fly_final, aes(x = flightless, y = gain_ile_met))+
  geom_boxplot()+
  xlab('Flightless birds groups')+
  ylab('Ile and met')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
c2 = ggplot(df_fly_final, aes(x = flightless, y = gain_ala))+
  geom_boxplot()+
  xlab('Flightless birds groups')+
  ylab('Ala')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
d2 = ggplot(df_fly_final, aes(x = flightless, y = lose_thr))+
  geom_boxplot()+
  xlab('Flightless birds groups')+
  ylab('Thr')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
third_graph = ggarrange(b2, d2, a2, c2,
                         ncol = 2, nrow = 2)
third_graph

a3 = ggplot(df_fly_final, aes(x = flightless, y = mega_gain_leu_phe))+
  geom_boxplot()+
  xlab('Flightless birds groups')+
  ylab('Phe and leu')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
b3 = ggplot(df_fly_final, aes(x = flightless, y = gain_ser))+
  geom_boxplot()+
  xlab('Flightless birds groups')+
  ylab('Ser')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
c3 = ggplot(df_fly_final, aes(x = flightless, y = gain_leu))+
  geom_boxplot()+
  xlab('Flightless birds groups')+
  ylab('Leu')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
d3 = ggplot(df_fly_final, aes(x = flightless, y = lose_pro))+
  geom_boxplot()+
  xlab('Flightless birds groups')+
  ylab('Pro')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")
fourth_graph = ggarrange(a3, b3, c3, d3,
                         nrow = 2, ncol = 2)
fourth_graph

df_dive_final = merge(df_dive, df_data, by = 'species_name')
df_dive_final = df_dive_final[df_dive_final$diving != 'waterbird',]
df_dive_final[df_dive_final$diving == '0',]$diving = 'Non-diving birds'

a = ggplot(df_dive_final, aes(x = diving, y = mega_gain_gly))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Gly')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
b = ggplot(df_dive_final, aes(x = diving, y = gain_serg_arg))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Ser and stop')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
c = ggplot(df_dive_final, aes(x = diving, y = gain_asp_glu))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Asp and glu')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
d = ggplot(df_dive_final, aes(x = diving, y = lose_asn_lys))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Asn and lys')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
fifth_graph = ggarrange(d,b,c,a,
                        nrow = 2, ncol = 2)
fifth_graph

a1 = ggplot(df_dive_final, aes(x = diving, y = mega_gain_cys_stop_trp))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Cys and Trp')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
b1 = ggplot(df_dive_final, aes(x = diving, y = gain_arg))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Arg')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
c1 = ggplot(df_dive_final, aes(x = diving, y = gain_tyr_stop))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Tyr and stop')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
d1 = ggplot(df_dive_final, aes(x = diving, y = lose_his_gln))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('His and gln')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
sixth_graph = ggarrange(c1, a1, d1, b1,
                        nrow = 2, ncol = 2)
sixth_graph


a2 = ggplot(df_dive_final, aes(x = diving, y = mega_gain_val))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Val')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
b2 = ggplot(df_dive_final, aes(x = diving, y = gain_ile_met))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Ile and met')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
c2 = ggplot(df_dive_final, aes(x = diving, y = gain_ala))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Ala')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
d2 = ggplot(df_dive_final, aes(x = diving, y = lose_thr))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Thr')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
seventh_graph = ggarrange(b2,d2,a2,c2,
                          ncol = 2, nrow = 2)
seventh_graph

a3 = ggplot(df_dive_final, aes(x = diving, y = mega_gain_leu_phe))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Leu and phe')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
b3 = ggplot(df_dive_final, aes(x = diving, y = gain_ser))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Ser')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
c3 = ggplot(df_dive_final, aes(x = diving, y = gain_leu))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Leu')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
d3 = ggplot(df_dive_final, aes(x = diving, y = lose_pro))+
  geom_boxplot()+
  xlab('Diving birds groups')+
  ylab('Pro')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim('Non-diving birds', "Anseriformes", "Sphenisciformes", "Podicipediformes", "Gaviiformes", "Suliformes", "Coraciiformes", "Passeriformes", "Gruiformes", "Charadriiformes", "Procellariiformes")
eighth_graph = ggarrange(a3,b3,c3,d3,
                         ncol = 2, nrow = 2)
eighth_graph

df_codons1$mega_gain_gly = sum(df_codons1$CCT + df_codons1$CCA + df_codons1$CCG + df_codons1$CCC)
mega_gain_gly = sum(df_codons1$CCT + df_codons1$CCA + df_codons1$CCG + df_codons1$CCC) #mega_gain_gly
gain_serg_arg = sum(df_codons1$TCT + df_codons1$TCA + df_codons1$TCG + df_codons1$TCC) #gain_serg_arg
gain_asp_glu = sum(df_codons1$CTT + df_codons1$CTA + df_codons1$CTG + df_codons1$CTC) #gain_asp_glu
lose_asn_lys = sum(df_codons1$TTT + df_codons1$TTA + df_codons1$TTG + df_codons1$TTC) #lose_asn_lys

mega_gain_cys_stop_trp = sum(df_codons1$ACT + df_codons1$ACA + df_codons1$ACG + df_codons1$ACC) #mega_gain_cys_stop_trp
gain_arg = sum(df_codons1$GCT + df_codons1$GCA + df_codons1$GCG + df_codons1$GCC) #gain_arg
gain_tyr_stop = sum(df_codons1$ATT + df_codons1$ATA + df_codons1$ATG + df_codons1$ATC) #gain_tyr_stop
lose_his_gln = sum(df_codons1$GTT + df_codons1$GTA + df_codons1$GTG + df_codons1$GTC) #lose_gis_gln

mega_gain_val = sum(df_codons1$CAT + df_codons1$CAA + df_codons1$CAG + df_codons1$CAC) #mega_gain_val
gain_ile_met = sum(df_codons1$TAT + df_codons1$TAA + df_codons1$TAG + df_codons1$TAC) #gain_ile_met
gain_ala = sum(df_codons1$CGT + df_codons1$CGA + df_codons1$CGG + df_codons1$CGC) #gain ala
lose_thr = sum(df_codons1$TGT + df_codons1$TGA + df_codons1$TGG + df_codons1$TGC) #lose_thr

mega_gain_leu_phe = sum(df_codons1$AAT + df_codons1$AAA + df_codons1$AAG + df_codons1$AAC) #mega_gain_leu_phe
gain_ser = sum(df_codons1$AGT + df_codons1$AGA + df_codons1$AGG + df_codons1$AGC) #gain_ser
gain_leu = sum(df_codons1$GAT + df_codons1$GAA + df_codons1$GAG + df_codons1$GAC) #gain_leu
lose_pro = sum(df_codons1$GGT + df_codons1$GGA + df_codons1$GGG + df_codons1$GGC)#lose pro
df_data = merge(df_data, df_fly)

df_fly_final = merge(df_fly_big, df_short)
df_fly_final = df_fly_final[df_fly_final$flightless != 'Galliformes',]
df_fly_final[df_fly_final$flightless == '0',]$flightless = 'Flying birds'
ggplot(df_data, aes(x = flightless, y = mega_gain_gly))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#first_sector = data.frame(mega_gain_cys_stop_trp, gain_arg, gain_tyr_stop , lose_his_gln)
first_sector = data.frame()
first_sector = rbind(first_sector,c("lose_his_gln", lose_his_gln))
first_sector = rbind(first_sector,c("gain_tyr_stop", gain_tyr_stop))
first_sector = rbind(first_sector,c("gain_arg", gain_arg))
first_sector = rbind(first_sector,c("mega_gain_cys_stop_trp", mega_gain_cys_stop_trp))
names(first_sector) = c("aa_type", "value")
first_sector$value = as.numeric(as.character(first_sector$value))
second_sector = data.frame()
second_sector = rbind(second_sector, c("lose_pro", lose_pro))
second_sector = rbind(second_sector, c("gain_leu", gain_leu))
second_sector = rbind(second_sector, c("gain_ser", gain_ser))
second_sector = rbind(second_sector, c("mega_gain_leu_phe", mega_gain_leu_phe))
names(second_sector) = c("aa_type", "value")
second_sector$value = as.numeric(as.character(second_sector$value))
third_sector = data.frame()
third_sector = rbind(third_sector, c("lose_thr", lose_thr))
third_sector = rbind(third_sector, c("gain_ala", gain_ala))  
third_sector = rbind(third_sector, c("gain_ile_met", gain_ile_met))
third_sector = rbind(third_sector, c("mega_gain_val", mega_gain_val))
names(third_sector) = c("aa_type", "value")
third_sector$value = as.numeric(as.character(third_sector$value))
fourth_sector = data.frame()
fourth_sector = rbind(fourth_sector, c("lose_asn_lys", lose_asn_lys))
fourth_sector = rbind(fourth_sector, c("gain_asp_glu", gain_asp_glu))
fourth_sector = rbind(fourth_sector, c("gain_serg_arg", gain_serg_arg))
fourth_sector = rbind(fourth_sector, c("mega_gain_gly", mega_gain_gly))
names(fourth_sector) = c("aa_type", "value")
fourth_sector$value = as.numeric(as.character(fourth_sector$value))
ggplot(first_sector, aes(x = aa_type, y = value))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(second_sector, aes(x = aa_type, y = value))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(third_sector, aes(x = aa_type, y = value))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(fourth_sector, aes(x = aa_type, y = value))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
