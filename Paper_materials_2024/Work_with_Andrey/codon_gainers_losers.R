rm(list = ls(all=TRUE))
library(ggplot2)
df_mtdna = read.csv('../../Paper_materials_2024/Birds_dataset_paper.csv')
df_codons = df_mtdna[,c(2,32:95)]
df_codons1 = aggregate(. ~ species_name, FUN = mean, data = df_codons)
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
